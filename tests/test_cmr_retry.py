"""
CMR granule searches (io_utils.search_data) are retried on transient errors:
HTTP 429/5xx and connection failures.  No network: a stand-in earthaccess
raises the way earthaccess 0.19 does -- RuntimeError(response.text) raised
FROM the requests HTTPError -- and a clock that only sleep() moves.

Background (2026-10-13..16): 59 MAAP DPS jobs lost a tile to "502 Bad Gateway"
from cmr.earthdata.nasa.gov/search, raised in find_ATL11_granules (38) and
tilingSchema.resolve_files_for_box (24); each tile succeeded when resubmitted.
"""
import sys
import types

import pytest
import requests

import pointCollection as pc
from pointCollection import io_utils


def http_error(status):
    """What earthaccess raises for a failed CMR request with this status."""
    response = requests.Response()
    response.status_code = status
    response._content = f'<title>{status}</title>'.encode()
    try:
        try:
            response.raise_for_status()
        except requests.exceptions.HTTPError as ex:
            raise RuntimeError(ex.response.text) from ex
    except RuntimeError as err:
        return err


@pytest.fixture
def cmr(monkeypatch):
    """A stand-in earthaccess whose search_data raises each error in
    state['errors'] in turn, then returns state['result']."""
    state = {'calls': [], 'errors': [], 'result': ['granule'], 'naps': [], 'now': 0.0}

    def search_data(**kwargs):
        state['calls'].append(kwargs)
        if state['errors']:
            raise state['errors'].pop(0)
        return state['result']

    def sleep(s):
        state['naps'].append(s)
        state['now'] += s

    module = types.ModuleType('earthaccess')
    module.search_data = search_data
    module.login = lambda strategy=None: None
    monkeypatch.setitem(sys.modules, 'earthaccess', module)
    monkeypatch.setattr('time.sleep', sleep)
    monkeypatch.setattr('time.monotonic', lambda: state['now'])
    monkeypatch.setattr('random.uniform', lambda a, b: (a + b) / 2)
    return state


def test_502_is_retried_then_succeeds(cmr, capsys):
    cmr['errors'] = [http_error(502), http_error(502)]
    assert io_utils.search_data(short_name='ATL11') == ['granule']
    assert len(cmr['calls']) == 3
    assert cmr['naps'] == list(io_utils.MAAP_BROKER_PAUSES_S[:2])
    err = capsys.readouterr().err
    assert 'CMR search failed (HTTP 502); retry 1/5' in err
    assert 'retry 2/5' in err


@pytest.mark.parametrize('status', [429, 500, 503, 504])
def test_other_transient_statuses_are_retried(cmr, status):
    cmr['errors'] = [http_error(status)]
    assert io_utils.search_data(short_name='ATL11') == ['granule']
    assert len(cmr['calls']) == 2


def test_connection_error_is_retried(cmr):
    cmr['errors'] = [requests.exceptions.ConnectionError('reset')]
    assert io_utils.search_data(short_name='ATL11') == ['granule']
    assert len(cmr['calls']) == 2


def test_bad_request_raises_at_once(cmr):
    cmr['errors'] = [http_error(400)]
    with pytest.raises(RuntimeError, match='400'):
        io_utils.search_data(short_name='ATL11')
    assert len(cmr['calls']) == 1 and cmr['naps'] == []


def test_other_exceptions_raise_at_once(cmr):
    cmr['errors'] = [ValueError('bad bounding box')]
    with pytest.raises(ValueError):
        io_utils.search_data(short_name='ATL11')
    assert len(cmr['calls']) == 1


def test_gives_up_with_the_last_error(cmr):
    """Unjittered, the pauses 10+20+40+60+60 fit the 240 s budget: all 6 tries
    are made, then the last 502 is raised unchanged."""
    cmr['errors'] = [http_error(502) for _ in range(10)]
    with pytest.raises(RuntimeError, match='502') as info:
        io_utils.search_data(short_name='ATL11')
    assert isinstance(info.value.__cause__, requests.exceptions.HTTPError)
    assert len(cmr['calls']) == len(io_utils.MAAP_BROKER_PAUSES_S) + 1
    assert sum(cmr['naps']) <= io_utils.MAAP_BROKER_BUDGET_S


def test_no_pause_past_the_budget(cmr, monkeypatch):
    """At the top of the jitter range the pauses would sum to 285 s, so the
    search stops before the pause that would end past 240 s."""
    monkeypatch.setattr('random.uniform', lambda a, b: b)
    cmr['errors'] = [http_error(502) for _ in range(10)]
    with pytest.raises(RuntimeError, match='502'):
        io_utils.search_data(short_name='ATL11')
    assert sum(cmr['naps']) <= io_utils.MAAP_BROKER_BUDGET_S
    assert len(cmr['calls']) < len(io_utils.MAAP_BROKER_PAUSES_S) + 1


def test_find_ATL11_granules_retries(cmr):
    from pointCollection.scripts.query_ATL11_cloud import find_ATL11_granules
    cmr['errors'] = [http_error(502)]
    assert find_ATL11_granules((-50, 60, -40, 70)) == ['granule']
    assert cmr['calls'][-1]['bounding_box'] == (-50, 60, -40, 70)
    assert len(cmr['calls']) == 2


def test_resolve_files_for_box_retries(cmr, monkeypatch):
    class Granule:
        def data_links(self, access=None):
            return ['s3://bucket/E0_N0.h5']
    cmr['result'] = [Granule()]
    cmr['errors'] = [http_error(502)]
    monkeypatch.setattr(io_utils, 'get_s3fs', lambda daac=None: 'fs')
    schema = pc.tilingSchema(tile_spacing=1.e5, mapping_function_name='round',
                             format_str='E%d_N%d', scale=1000, extension='.h5',
                             directory=None,
                             source={'type': 'EarthAccess', 'short_name': 'X', 'daac': 'NSIDC'})
    resolved, fs = schema.resolve_files_for_box([[0, 0], [0, 0]])
    assert resolved == {'E0_N0.h5': 's3://bucket/E0_N0.h5'}
    assert len(cmr['calls']) == 2
