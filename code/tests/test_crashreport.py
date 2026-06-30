"""Unit tests for the crash/error reporter (``crashreport.py``).

A real exception is manufactured (raise + catch) so there's a genuine
traceback object to format; the GUI dialog is replaced by a recording
callback, and the excepthook is exercised by calling it directly.
"""

import sys
import urllib.parse

import crashreport as cr


def _make_exc_info():
    try:
        raise ValueError('boom happened in feature 0.80_418.1451n')
    except ValueError:
        return sys.exc_info()


# --------------------------------------------------------------------------- #
# format_report
# --------------------------------------------------------------------------- #

def test_format_report_contains_traceback_and_environment():
    et, ev, tb = _make_exc_info()
    report = cr.format_report(et, ev, tb, now=0)
    assert 'MPACT crash report' in report
    assert 'ValueError: boom happened' in report
    assert 'Traceback' in report
    assert 'Python:' in report and 'Platform:' in report


def test_format_report_includes_context_when_given():
    et, ev, tb = _make_exc_info()
    report = cr.format_report(et, ev, tb, context='dataset=PTY087I2; step=heatmap')
    assert 'Context:' in report
    assert 'PTY087I2' in report


# --------------------------------------------------------------------------- #
# one_line_summary
# --------------------------------------------------------------------------- #

def test_one_line_summary_uses_type_and_message():
    et, ev, _ = _make_exc_info()
    assert cr.one_line_summary(et, ev) == 'ValueError: boom happened in feature 0.80_418.1451n'


def test_one_line_summary_handles_empty_message():
    try:
        raise RuntimeError()
    except RuntimeError:
        et, ev, _ = sys.exc_info()
    assert cr.one_line_summary(et, ev) == 'RuntimeError'


# --------------------------------------------------------------------------- #
# write_log
# --------------------------------------------------------------------------- #

def test_write_log_creates_timestamped_file(tmp_path):
    path = cr.write_log('hello report', tmp_path / 'crashlogs', now=0)
    assert path is not None
    with open(path) as f:
        assert f.read() == 'hello report'
    assert 'mpact_crash_' in path and path.endswith('.log')


def test_write_log_returns_none_on_failure(tmp_path):
    # Point log_dir at a path that exists as a *file* so makedirs fails.
    afile = tmp_path / 'not_a_dir'
    afile.write_text('x')
    assert cr.write_log('r', str(afile / 'sub')) is None


# --------------------------------------------------------------------------- #
# build_issue_url
# --------------------------------------------------------------------------- #

def test_build_issue_url_is_wellformed_and_encoded():
    url = cr.build_issue_url('a traceback & stuff', 'Crash: ValueError: boom',
                             repo='robertsamples/mpact')
    assert url.startswith('https://github.com/robertsamples/mpact/issues/new?')
    parsed = urllib.parse.urlparse(url)
    params = urllib.parse.parse_qs(parsed.query)
    assert params['title'] == ['Crash: ValueError: boom']
    assert 'a traceback & stuff' in params['body'][0]


def test_build_issue_url_truncates_huge_body():
    huge = 'x' * 50000
    url = cr.build_issue_url(huge, 'Crash', repo='r/m')
    params = urllib.parse.parse_qs(urllib.parse.urlparse(url).query)
    assert 'truncated' in params['body'][0]
    assert len(params['body'][0]) < 7000


# --------------------------------------------------------------------------- #
# excepthook
# --------------------------------------------------------------------------- #

def test_excepthook_invokes_handler_and_writes_log(tmp_path):
    received = {}

    def handler(report, log_path, issue_url):
        received['report'] = report
        received['log_path'] = log_path
        received['issue_url'] = issue_url

    chained = []
    hook = cr.make_excepthook(handler, log_dir=str(tmp_path / 'logs'),
                              repo='robertsamples/mpact',
                              prev_hook=lambda *a: chained.append(a))
    et, ev, tb = _make_exc_info()
    hook(et, ev, tb)

    assert 'ValueError' in received['report']
    assert received['log_path'] is not None and received['log_path'].endswith('.log')
    assert received['issue_url'].startswith('https://github.com/robertsamples/mpact/issues/new?')
    # The previous hook was still chained (traceback reaches the console).
    assert len(chained) == 1


def test_excepthook_swallows_handler_errors(tmp_path):
    def bad_handler(report, log_path, issue_url):
        raise RuntimeError('dialog blew up')

    hook = cr.make_excepthook(bad_handler, log_dir=str(tmp_path), prev_hook=lambda *a: None)
    et, ev, tb = _make_exc_info()
    # Must not raise, even though the handler does.
    hook(et, ev, tb)


def test_excepthook_uses_context_provider(tmp_path):
    received = {}

    def handler(report, log_path, issue_url):
        received['report'] = report

    hook = cr.make_excepthook(handler, log_dir=str(tmp_path),
                              context_provider=lambda: 'active dataset: foo',
                              prev_hook=lambda *a: None)
    et, ev, tb = _make_exc_info()
    hook(et, ev, tb)
    assert 'active dataset: foo' in received['report']


def test_install_excepthook_restores(tmp_path):
    original = sys.excepthook
    try:
        prev = cr.install_excepthook(lambda *a: None, log_dir=str(tmp_path))
        assert prev is original
        assert sys.excepthook is not original
    finally:
        sys.excepthook = original
