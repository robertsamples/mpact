"""Unit tests for the MPACT self-update checker (``mpactupdate.py``).

The GitHub API and the git call are both injected, so no network or
subprocess is touched.
"""

import io
import json

import mpactupdate as mu


class _FakeResponse(io.BytesIO):
    def __enter__(self):
        return self

    def __exit__(self, *exc):
        self.close()
        return False


def _release_opener(tag, *, html_url='http://x/rel', body='notes', record=None):
    payload = json.dumps({'tag_name': tag, 'html_url': html_url, 'body': body}).encode()

    def opener(request, timeout=None):
        if record is not None:
            # request is a urllib Request; record its full URL.
            record.append(getattr(request, 'full_url', request))
        return _FakeResponse(payload)
    return opener


# --------------------------------------------------------------------------- #
# version comparison
# --------------------------------------------------------------------------- #

def test_is_newer_basic():
    assert mu.is_newer('2.1.0', '2.0.0') is True
    assert mu.is_newer('2.0.0', '2.0.0') is False
    assert mu.is_newer('1.9.0', '2.0.0') is False


def test_is_newer_strips_v_prefix():
    assert mu.is_newer('v2.1.0', '2.0.0') is True
    assert mu.is_newer('V2.0.1', 'v2.0.0') is True


def test_is_newer_numeric_not_lexicographic():
    # 2.10.0 > 2.9.0 numerically (lexicographically it would be "<").
    assert mu.is_newer('2.10.0', '2.9.0') is True


def test_is_newer_unparseable_tag_is_not_newer():
    assert mu.is_newer('not-a-version', '2.0.0') is False


# --------------------------------------------------------------------------- #
# release fetch + check
# --------------------------------------------------------------------------- #

def test_check_reports_available_update():
    info = mu.check_for_update(current_version='2.0.0',
                               opener=_release_opener('v2.5.0'))
    assert info.available is True
    assert info.latest == 'v2.5.0'
    assert info.url == 'http://x/rel'
    assert info.notes == 'notes'


def test_check_reports_no_update_when_same_version():
    info = mu.check_for_update(current_version='2.0.0',
                               opener=_release_opener('v2.0.0'))
    assert info.available is False
    assert info.current == '2.0.0'


def test_check_hits_the_configured_repo():
    seen = []
    mu.check_for_update(current_version='2.0.0', repo='robertsamples/mpact',
                        opener=_release_opener('v2.0.0', record=seen))
    assert seen == ['https://api.github.com/repos/robertsamples/mpact/releases/latest']


def test_check_is_safe_when_offline():
    def failing_opener(request, timeout=None):
        raise OSError('no network')
    info = mu.check_for_update(current_version='2.0.0', opener=failing_opener)
    assert info.available is False
    assert info.latest is None


def test_check_is_safe_when_no_releases_yet():
    # GitHub returns 404 -> urlopen raises HTTPError -> fetch returns None.
    def opener_404(request, timeout=None):
        raise OSError('HTTP 404')
    info = mu.check_for_update(current_version='2.0.0', opener=opener_404)
    assert info.available is False


def test_fetch_handles_malformed_json():
    def opener(request, timeout=None):
        return _FakeResponse(b'not json at all')
    assert mu.fetch_latest_release(opener=opener) is None


# --------------------------------------------------------------------------- #
# git update
# --------------------------------------------------------------------------- #

class _Completed:
    def __init__(self, returncode, stdout='', stderr=''):
        self.returncode = returncode
        self.stdout = stdout
        self.stderr = stderr


def test_apply_git_update_success():
    calls = []

    def runner(cmd, **kwargs):
        calls.append(cmd)
        return _Completed(0, stdout='Updating abc..def\nFast-forward\n')

    ok, output = mu.apply_git_update('/repo', runner=runner)
    assert ok is True
    assert 'Fast-forward' in output
    assert calls[0][:3] == ['git', '-C', '/repo']
    assert 'pull' in calls[0] and '--ff-only' in calls[0]


def test_apply_git_update_reports_failure():
    def runner(cmd, **kwargs):
        return _Completed(1, stderr='error: local changes would be overwritten')
    ok, output = mu.apply_git_update('/repo', runner=runner)
    assert ok is False
    assert 'local changes' in output


def test_apply_git_update_handles_missing_git():
    def runner(cmd, **kwargs):
        raise FileNotFoundError('git not found')
    ok, output = mu.apply_git_update('/repo', runner=runner)
    assert ok is False
    assert 'could not be run' in output
