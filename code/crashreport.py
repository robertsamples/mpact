"""
MPACT
Copyright 2022, Robert M. Samples, Sara P. Puckett, and Marcy J. Balunas

Qt-free crash/error reporting. Installs a ``sys.excepthook`` that, on any
otherwise-unhandled exception:

1. formats a full report (traceback + environment: MPACT/Python/platform
   versions, timestamp, optional context such as the tail of the run log),
2. writes it to a timestamped file under a crash-log directory (so there's a
   durable record even if the user dismisses the dialog), and
3. hands the report to a GUI callback that asks the user whether to send it.

The "send" path is deliberately backend-free: it builds a pre-filled GitHub
*new issue* URL (title + body) for the MPACT repo, so reporting is one click
in the browser and nothing leaves the user's machine until they choose to
submit it. That satisfies "prompt the user before sending" without any cloud
egress, DSN, or account.

Why not Sentry (the obvious off-the-shelf option): ``sentry-sdk`` is excellent
for hosted/web services but (a) sends events to a Sentry project by default --
exactly the silent-egress this tool should avoid for a desktop research app,
(b) needs a DSN/account to be provisioned, and (c) still needs a custom
``before_send`` hook + dialog to honour "ask first." For a single-user desktop
tool the local-log + pre-filled-GitHub-issue flow gives the same practical
benefit (a complete traceback in the maintainer's hands) with no infrastructure
and no privacy surprise. If MPACT ever ships to many non-technical users and a
central error feed becomes worth it, Sentry with ``before_send`` gating is the
documented upgrade path.

This module is Qt-free and unit-tested (see ``tests/test_crashreport.py``); the
GUI dialog is injected as a plain callback.
"""

import os
import platform
import sys
import time
import traceback
import urllib.parse

DEFAULT_REPO = 'robertsamples/mpact'
# GitHub rejects extremely long issue URLs; keep the prefilled body well under
# the practical limit so the link always opens (the full report is always in
# the log file regardless).
_MAX_ISSUE_BODY = 6000


def _app_version():
    try:
        from mpactupdate import __version__
        return __version__
    except Exception:
        return 'unknown'


def format_report(exc_type, exc_value, exc_tb, context=None, now=None):
    """Build the human-readable crash report text.

    Args:
        exc_type/exc_value/exc_tb: the ``sys.exc_info()``-style triple.
        context: optional extra text appended under a "Context" heading
            (e.g. the last lines of the run log, the current dataset name).
        now: epoch seconds for the timestamp (injectable for tests).

    Returns:
        A multi-section plain-text report.
    """
    now = time.time() if now is None else now
    stamp = time.strftime('%Y-%m-%d %H:%M:%S', time.localtime(now))
    tb_text = ''.join(traceback.format_exception(exc_type, exc_value, exc_tb))
    lines = [
        'MPACT crash report',
        '==================',
        'Time: ' + stamp,
        'MPACT version: ' + _app_version(),
        'Python: ' + sys.version.split()[0],
        'Platform: ' + platform.platform(),
        '',
        'Traceback:',
        tb_text.rstrip(),
    ]
    if context:
        lines += ['', 'Context:', str(context).rstrip()]
    return '\n'.join(lines) + '\n'


def write_log(report, log_dir, now=None):
    """Write ``report`` to a timestamped file under ``log_dir``.

    Creates ``log_dir`` if needed. Returns the path written, or ``None`` if the
    write failed (reporting must never raise from inside an excepthook).
    """
    now = time.time() if now is None else now
    try:
        os.makedirs(log_dir, exist_ok=True)
        fname = 'mpact_crash_' + time.strftime('%Y%m%d_%H%M%S', time.localtime(now)) + '.log'
        path = os.path.join(log_dir, fname)
        with open(path, 'w', encoding='utf-8', errors='replace') as handle:
            handle.write(report)
        return path
    except Exception:
        return None


def one_line_summary(exc_type, exc_value):
    """A concise ``TypeName: message`` for use as an issue title."""
    name = getattr(exc_type, '__name__', str(exc_type))
    message = str(exc_value).strip().splitlines()[0] if str(exc_value).strip() else ''
    return (name + ': ' + message).strip().rstrip(':').strip() if message else name


def build_issue_url(report, title, repo=DEFAULT_REPO):
    """Build a GitHub 'new issue' URL with a prefilled title and body.

    The body is the report wrapped in a code fence and truncated to
    :data:`_MAX_ISSUE_BODY` so the URL stays openable. The full untruncated
    report always lives in the on-disk log.
    """
    body = report
    if len(body) > _MAX_ISSUE_BODY:
        body = body[:_MAX_ISSUE_BODY] + '\n...\n[truncated -- see attached crash log]'
    body_md = ('**Describe what you were doing when this happened:**\n\n\n'
               '---\n```\n' + body + '\n```\n')
    query = urllib.parse.urlencode({'title': title, 'body': body_md})
    return 'https://github.com/' + repo + '/issues/new?' + query


def make_excepthook(report_handler, log_dir=None, repo=DEFAULT_REPO,
                    context_provider=None, prev_hook=None):
    """Build (but don't install) an excepthook.

    Args:
        report_handler: callable ``handler(report, log_path, issue_url)`` that
            shows the user the report and offers to send it. Exceptions raised
            by the handler are swallowed (an excepthook must not itself raise).
        log_dir: directory for crash logs (skipped if None).
        repo: GitHub repo for the prefilled issue URL.
        context_provider: optional zero-arg callable returning extra context
            text to embed (called defensively; failure is ignored).
        prev_hook: a previous excepthook to chain to (defaults to the standard
            ``sys.__excepthook__`` so the traceback still reaches the console).

    Returns:
        A function with the ``(exc_type, exc_value, exc_tb)`` signature.
    """
    prev_hook = prev_hook if prev_hook is not None else sys.__excepthook__

    def _hook(exc_type, exc_value, exc_tb):
        # Always let the default hook print to stderr first (and never let our
        # own reporting suppress that or raise over it).
        try:
            prev_hook(exc_type, exc_value, exc_tb)
        except Exception:
            pass
        try:
            context = None
            if context_provider is not None:
                try:
                    context = context_provider()
                except Exception:
                    context = None
            report = format_report(exc_type, exc_value, exc_tb, context=context)
            log_path = write_log(report, log_dir) if log_dir else None
            title = 'Crash: ' + one_line_summary(exc_type, exc_value)
            issue_url = build_issue_url(report, title, repo=repo)
            if report_handler is not None:
                report_handler(report, log_path, issue_url)
        except Exception:
            # Reporting failed -- the default hook already printed the real
            # traceback, so just give up quietly rather than masking it.
            pass

    return _hook


def install_excepthook(report_handler, log_dir=None, repo=DEFAULT_REPO,
                       context_provider=None):
    """Install the crash excepthook as ``sys.excepthook``; return the previous
    hook (so callers can restore it)."""
    prev = sys.excepthook
    sys.excepthook = make_excepthook(
        report_handler, log_dir=log_dir, repo=repo,
        context_provider=context_provider, prev_hook=prev)
    return prev
