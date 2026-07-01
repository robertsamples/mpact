"""
MPACT
Copyright 2026, Robert M. Samples

Qt-free self-update checker. Queries the GitHub Releases API for the MPACT
repository (Robert Samples' fork, ``robertsamples/mpact``, by default),
compares the latest published release tag against the locally-running version,
and -- if newer -- hands the GUI the information it needs to ask the user
whether to update. The actual update of a git checkout is a best-effort
``git pull --ff-only``.

Why not an off-the-shelf updater framework: the established Python option,
``pyupdater``, targets *frozen* (PyInstaller/cx_Freeze) apps and needs its own
patch-server + signing infrastructure -- heavyweight for a tool that's run
from a git clone (and the portable PyInstaller build is a separate, infrequent
artifact). For a source checkout the meaningful "update" is ``git pull``, and
"is there a newer release" is one GitHub API call + a version compare. That's
what this module does, with no third-party dependency beyond ``packaging``
(already present; a tuple-based fallback covers its absence).

This module is Qt-free and unit-tested (see ``tests/test_mpactupdate.py``); the
network and the git call are both injectable so the tests never touch either.
"""

import json
import subprocess
import urllib.request

#: The version of MPACT currently running. Bump this when cutting a release,
#: and create a matching GitHub release/tag (e.g. ``v1.0.1``) on the fork so
#: this checker can see it. Kept here as the single in-code source of truth;
#: keep it consistent with main.py's ``label_credits`` display string
#: (currently shows ``v1.00.02``).
__version__ = '1.0.02'

DEFAULT_REPO = 'robertsamples/mpact'
_RELEASES_LATEST = 'https://api.github.com/repos/{repo}/releases/latest'


class UpdateInfo:
    """Result of an update check.

    Attributes:
        available: True if the latest release is newer than the running version.
        current: the running version string.
        latest: the latest release's tag (raw, e.g. ``v2.1.0``) or None.
        url: the release's html_url (for "view release" / manual download).
        notes: the release body/changelog text (may be '').
    """
    __slots__ = ('available', 'current', 'latest', 'url', 'notes')

    def __init__(self, available, current, latest=None, url=None, notes=''):
        self.available = available
        self.current = current
        self.latest = latest
        self.url = url
        self.notes = notes

    def __repr__(self):
        return ('UpdateInfo(available=%r, current=%r, latest=%r)'
                % (self.available, self.current, self.latest))


def _normalize(tag):
    """Strip a leading 'v'/'V' and surrounding whitespace from a tag."""
    tag = str(tag).strip()
    if tag[:1] in ('v', 'V'):
        tag = tag[1:]
    return tag


def is_newer(latest_tag, current_version):
    """True if ``latest_tag`` represents a newer version than ``current_version``.

    Uses ``packaging.version`` when available (PEP 440 aware), falling back to a
    dotted-integer tuple comparison so a missing ``packaging`` never breaks the
    check. A tag that can't be parsed at all is treated as "not newer" (fail
    safe -- never nag about an unparseable tag).
    """
    latest = _normalize(latest_tag)
    current = _normalize(current_version)
    try:
        from packaging.version import parse as _parse
        return _parse(latest) > _parse(current)
    except Exception:
        def _tuple(text):
            parts = []
            for chunk in text.replace('-', '.').split('.'):
                if chunk.isdigit():
                    parts.append(int(chunk))
                else:
                    break
            return tuple(parts)
        lt, ct = _tuple(latest), _tuple(current)
        if not lt:
            return False
        return lt > ct


def fetch_latest_release(repo=DEFAULT_REPO, opener=None, timeout=10):
    """Fetch the latest published release from the GitHub API.

    Returns the parsed JSON dict, or ``None`` if the repo has no published
    releases (the API 404s) or the response can't be parsed. ``opener`` is an
    injectable ``opener(request_or_url, timeout=...)`` returning a readable,
    context-managed response (defaults to ``urllib.request.urlopen``).
    """
    opener = opener if opener is not None else urllib.request.urlopen
    url = _RELEASES_LATEST.format(repo=repo)
    # GitHub recommends a User-Agent; the v3 Accept header pins the response shape.
    request = urllib.request.Request(url, headers={
        'Accept': 'application/vnd.github+json',
        'User-Agent': 'MPACT-update-check',
    })
    try:
        try:
            response = opener(request, timeout=timeout)
        except TypeError:
            response = opener(request)
        with response as resp:
            payload = resp.read()
    except Exception:
        # Network error, 404 (no releases yet), DNS failure, offline -- all
        # non-fatal: an update check must never break startup.
        return None
    try:
        data = json.loads(payload.decode('utf-8') if isinstance(payload, bytes) else payload)
    except Exception:
        return None
    if not isinstance(data, dict) or 'tag_name' not in data:
        return None
    return data


def check_for_update(current_version=__version__, repo=DEFAULT_REPO, opener=None,
                     timeout=10):
    """Check whether a newer MPACT release exists.

    Returns an :class:`UpdateInfo`. ``available`` is False (and ``latest`` is
    None) when the check can't reach the API or finds no newer release -- the
    GUI can simply do nothing in that case.
    """
    release = fetch_latest_release(repo=repo, opener=opener, timeout=timeout)
    if release is None:
        return UpdateInfo(available=False, current=current_version)
    latest_tag = release.get('tag_name')
    return UpdateInfo(
        available=is_newer(latest_tag, current_version),
        current=current_version,
        latest=latest_tag,
        url=release.get('html_url'),
        notes=release.get('body') or '',
    )


def apply_git_update(repo_dir, runner=None, remote='origin', branch='main'):
    """Best-effort ``git pull --ff-only`` of a source checkout.

    Only meaningful when MPACT is running from a git clone (not a frozen
    build). ``runner`` is an injectable ``subprocess.run``-compatible callable
    (so tests don't shell out). Returns ``(success, output)`` where ``output``
    is combined stdout/stderr text.
    """
    runner = runner if runner is not None else subprocess.run
    try:
        completed = runner(
            ['git', '-C', str(repo_dir), 'pull', '--ff-only', remote, branch],
            capture_output=True, text=True, timeout=120,
        )
    except Exception as exc:
        return False, 'git pull could not be run: ' + str(exc)
    output = (getattr(completed, 'stdout', '') or '') + (getattr(completed, 'stderr', '') or '')
    return completed.returncode == 0, output.strip()
