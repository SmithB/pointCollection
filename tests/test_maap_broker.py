"""
The MAAP-brokered DAAC credentials in io_utils.get_s3fs(): the broker call is
retried, and when neither the broker nor earthaccess has credentials the error
says why.  No network: stand-ins replace maap-py and earthaccess.

Background (2026-10-01, 556 jobs submitted to MAAP DPS at once): 61 lost their
single broker attempt, fell back to earthaccess -- which has no login on a
worker -- and died with "'NoneType' object has no attribute
'get_s3_filesystem'", the real reason hidden by a process-wide
warnings.filterwarnings("ignore") in ps_scale_for_lat.py.
"""
import subprocess
import sys
import types
import warnings

import numpy as np
import pytest

import pointCollection as pc
from pointCollection import io_utils

CREDS = {'accessKeyId': 'id', 'secretAccessKey': 'secret', 'sessionToken': 'token',
         'expiration': '2099-01-01 00:00:00+00:00'}


@pytest.fixture
def broker(monkeypatch):
    """A stand-in maap-py whose broker fails `failures` times (with `error`),
    then answers; whose client fails to build `build_failures` times; and a
    clock that only sleep() moves, with the jitter fixed at `jitter`."""
    state = {'calls': 0, 'failures': 0, 'naps': [], 'builds': 0, 'build_failures': 0,
             'error': ConnectionError('timed out'), 'now': 0.0, 'jitter': 1.0}

    class AWS:
        def earthdata_s3_credentials(self, endpoint):
            state['calls'] += 1
            if state['calls'] <= state['failures']:
                raise state['error']
            return dict(CREDS)

    class MAAP:
        def __init__(self, maap_host=None):
            state['builds'] += 1
            if state['builds'] <= state['build_failures']:
                raise ConnectionError('config timed out')
            self.aws = AWS()

    def sleep(s):
        state['naps'].append(s)
        state['now'] += s

    package, module = types.ModuleType('maap'), types.ModuleType('maap.maap')
    module.MAAP = MAAP
    monkeypatch.setitem(sys.modules, 'maap', package)
    monkeypatch.setitem(sys.modules, 'maap.maap', module)
    monkeypatch.setenv('MAAP_PGT', 'set')
    monkeypatch.setattr('time.sleep', sleep)
    monkeypatch.setattr('time.monotonic', lambda: state['now'])
    monkeypatch.setattr('random.uniform', lambda a, b: state['jitter'] * (a + b) / 2)
    monkeypatch.setattr(io_utils, '_s3fs_with_credentials', lambda creds, **kw: ('fs', creds))
    monkeypatch.setattr(io_utils, '_S3FS_CACHE', {})
    monkeypatch.setattr(io_utils, '_BROKER_FAILURES', {})
    return state


def fake_earthaccess(monkeypatch, session):
    module = types.ModuleType('earthaccess')

    def get_s3fs_session(daac=None, **kwargs):
        if isinstance(session, Exception):
            raise session
        return session
    module.get_s3fs_session = get_s3fs_session
    monkeypatch.setitem(sys.modules, 'earthaccess', module)


def test_the_broker_call_is_retried(broker):
    broker['failures'] = 2
    fs, expires = io_utils._s3fs_from_maap('NSIDC')
    assert fs == ('fs', CREDS) and expires is not None
    assert broker['calls'] == 3
    assert broker['naps'] == list(io_utils.MAAP_BROKER_PAUSES_S[:2])
    assert broker['builds'] == 1                    # one client for every try


def test_giving_up_warns_with_the_reason(broker):
    broker['failures'] = 99
    with pytest.warns(UserWarning, match=r'in 6 attempts \(last error: ConnectionError: timed out\)'):
        assert io_utils._s3fs_from_maap('NSIDC') == (None, None)
    assert broker['calls'] == len(io_utils.MAAP_BROKER_PAUSES_S) + 1
    assert broker['naps'] == list(io_utils.MAAP_BROKER_PAUSES_S)        # none after the last
    assert sum(broker['naps']) <= io_utils.MAAP_BROKER_BUDGET_S


def test_pauses_are_jittered(broker):
    broker['failures'], broker['jitter'] = 1, 0.5
    io_utils._s3fs_from_maap('NSIDC')
    assert broker['naps'] == [io_utils.MAAP_BROKER_PAUSES_S[0] * 0.5]


def test_no_pause_past_the_budget(broker):
    # jitter at its top: 15 + 30 + 60 + 90 = 195 s, and the next 90 s would end at 285
    broker['failures'], broker['jitter'] = 99, 1.5
    with pytest.warns(UserWarning, match=r'in 5 attempts .*a 90 s pause would pass the 240 s budget'):
        assert io_utils._s3fs_from_maap('NSIDC') == (None, None)
    assert broker['naps'] == [15, 30, 60, 90]


def test_401_stops_at_once(broker):
    class HTTPError(Exception):
        response = types.SimpleNamespace(status_code=401)
    broker['failures'], broker['error'] = 99, HTTPError('401 Unauthorized')
    with pytest.warns(UserWarning, match=r'in 1 attempts .*HTTP 401: MAAP_PGT was rejected'):
        assert io_utils._s3fs_from_maap('NSIDC') == (None, None)
    assert broker['calls'] == 1 and broker['naps'] == []


def test_a_client_that_failed_to_build_is_built_again(broker):
    broker['build_failures'] = 2
    fs, _ = io_utils._s3fs_from_maap('NSIDC')
    assert fs == ('fs', CREDS)
    assert broker['builds'] == 3 and broker['calls'] == 1


def test_no_broker_and_no_earthaccess_login_names_both(broker, monkeypatch):
    # a DPS worker: this is the AttributeError earthaccess raises with no login
    broker['failures'] = 99
    fake_earthaccess(monkeypatch, AttributeError("'NoneType' object has no attribute 'get_s3_filesystem'"))
    with pytest.warns(UserWarning), pytest.raises(RuntimeError) as err:
        io_utils.get_s3fs(daac='NSIDC')
    message = str(err.value)
    assert 'no NSIDC S3 credentials' in message
    assert 'ConnectionError: timed out' in message
    assert 'earthaccess fallback has no login here (AttributeError' in message


def test_earthaccess_still_serves_when_the_broker_is_down(broker, monkeypatch):
    # the ADE with an Earthdata login: the fallback is still the answer
    broker['failures'] = 99
    fake_earthaccess(monkeypatch, 'earthaccess session')
    with pytest.warns(UserWarning, match='falling back to earthaccess'):
        assert io_utils.get_s3fs(daac='NSIDC') == 'earthaccess session'


def test_off_maap_an_earthaccess_error_is_left_alone(monkeypatch):
    monkeypatch.delenv('MAAP_PGT', raising=False)
    monkeypatch.setattr(io_utils, '_S3FS_CACHE', {})
    monkeypatch.setattr(io_utils, '_BROKER_FAILURES', {})
    fake_earthaccess(monkeypatch, ValueError('not logged in'))
    with pytest.raises(ValueError, match='not logged in'):
        io_utils.get_s3fs(daac='NSIDC')


def test_a_later_success_forgets_the_earlier_failure(broker):
    broker['failures'] = 99
    with pytest.warns(UserWarning):
        io_utils._s3fs_from_maap('NSIDC')
    assert 'NSIDC' in io_utils._BROKER_FAILURES
    broker['calls'], broker['failures'] = 0, 0
    io_utils._s3fs_from_maap('NSIDC')
    assert 'NSIDC' not in io_utils._BROKER_FAILURES


def test_importing_pointCollection_does_not_silence_warnings():
    code = ('import warnings; before = list(warnings.filters); import pointCollection;'
            'added = [f for f in warnings.filters if f not in before];'
            # a BLANKET ignore: no message pattern, every category
            'assert not [f for f in added if f[0] == "ignore" and f[1] is None and f[2] is Warning], added')
    result = subprocess.run([sys.executable, '-c', code], capture_output=True, text=True)
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize('lat', [np.array([90., 89.9, 70.]), np.array([-90., -71.]), 90.0,
                                 np.array([np.nan, np.nan])])
def test_ps_scale_for_lat_is_quiet_on_its_own(lat):
    with warnings.catch_warnings():
        warnings.simplefilter('error')
        scale = pc.ps_scale_for_lat(lat)
    assert np.shape(scale) == np.shape(lat)


def test_ps_scale_for_lat_values_are_unchanged():
    # computed with the module-level filter still in place (main, 42f66cc)
    np.testing.assert_allclose(pc.ps_scale_for_lat(np.array([90., 70., 60.])),
                               [1.03107857, 1.0, 0.96206753], rtol=1e-6)
