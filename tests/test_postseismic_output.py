'''
Tests for multi-window cumulative postseismic displacement output.

`eqpostdisp` / `eqpostdsp_tspan` may be a single file/window or lists that are
paired positionally. Each output file must carry its own tstart/tend header,
and an existing file must not silently receive rows from a different window.
`eqpostts_tspan` controls the separate postseismic time-series window.
'''
import numpy as np
import pytest

from pytsfit.output import (HEADER_POSTDISP, open_postdisp_file,
                            output_postseismic_ts, pair_postseismic_outputs,
                            postdisp_header, resolve_displacement_spans,
                            resolve_postseismic_ts_span)


# ---------------------------------------------------------------------------
# pair_postseismic_outputs: old and new config forms
# ---------------------------------------------------------------------------
def test_pair_old_single_file_and_window():
    pairs = pair_postseismic_outputs('postseismic.gmtvec', [2020.557, 2026.5])
    assert pairs == [('postseismic.gmtvec', [2020.557, 2026.5])]


def test_pair_new_lists_are_paired_by_position():
    pairs = pair_postseismic_outputs(
        ['postseismic1.gmtvec', 'postseismic2.gmtvec'],
        [[2020.557, 2021.572], [2021.573, 2026.5]])
    assert pairs == [('postseismic1.gmtvec', [2020.557, 2021.572]),
                     ('postseismic2.gmtvec', [2021.573, 2026.5])]


def test_pair_accepts_numpy_arrays():
    pairs = pair_postseismic_outputs(
        np.array(['a.gmtvec', 'b.gmtvec']),
        np.array([[2020.0, 2021.0], [2021.0, 2022.0]]))
    assert pairs == [('a.gmtvec', [2020.0, 2021.0]),
                     ('b.gmtvec', [2021.0, 2022.0])]


def test_pair_empty_file_name_disables_output():
    assert pair_postseismic_outputs('', [2020.0, 2021.0]) == []
    assert pair_postseismic_outputs([], [[2020.0, 2021.0]]) == []


def test_pair_file_without_window_warns_and_returns_empty(caplog):
    with caplog.at_level('WARNING'):
        assert pair_postseismic_outputs('postseismic.gmtvec', []) == []
    assert 'eqpostdsp_tspan is empty' in caplog.text


def test_pair_length_mismatch_raises():
    with pytest.raises(ValueError, match='paired one-to-one'):
        pair_postseismic_outputs(['a.gmtvec', 'b.gmtvec'],
                                 [[2020.0, 2021.0]])


def test_pair_rejects_bad_window():
    with pytest.raises(ValueError, match='t0 < t1'):
        pair_postseismic_outputs('a.gmtvec', [2021.0, 2020.0])
    with pytest.raises(ValueError, match=r'\[t0, t1\] pair'):
        pair_postseismic_outputs('a.gmtvec', [[2020.0, 2021.0, 2022.0]])


def test_pair_rejects_duplicate_file_with_different_windows():
    with pytest.raises(ValueError, match='different file for each window'):
        pair_postseismic_outputs(['a.gmtvec', 'a.gmtvec'],
                                 [[2020.0, 2021.0], [2021.0, 2022.0]])


# ---------------------------------------------------------------------------
# Header
# ---------------------------------------------------------------------------
def test_postdisp_header_records_window_and_columns():
    header = postdisp_header([2020.557, 2021.572])
    lines = header.splitlines()
    assert lines[0].startswith('# Postseismic cumulative displacement')
    assert 'tstart(yr) = 2020.557' in lines[1]
    assert 'tend(yr) = 2021.572' in lines[1]
    assert 'tspan(yr) = 1.015' in lines[1]
    assert lines[2] == HEADER_POSTDISP.rstrip('\n')
    assert all(line.startswith('#') for line in lines)


# ---------------------------------------------------------------------------
# File handling: header on first write, window validation on append
# ---------------------------------------------------------------------------
def test_new_file_gets_window_header(tmp_path):
    path = tmp_path / 'post.gmtvec'
    fid = open_postdisp_file(str(path), [2020.557, 2021.572])
    fid.write(' 1.0 2.0\n')
    fid.close()
    text = path.read_text()
    assert 'tstart(yr) = 2020.557' in text
    assert 'tend(yr) = 2021.572' in text
    assert text.rstrip().endswith('1.0 2.0')


def test_append_same_window_succeeds_without_duplicate_header(tmp_path):
    path = tmp_path / 'post.gmtvec'
    for _ in range(2):
        fid = open_postdisp_file(str(path), [2020.557, 2021.572])
        fid.write(' 1.0 2.0\n')
        fid.close()
    text = path.read_text()
    assert text.count('tstart(yr)') == 1
    assert text.count(' 1.0 2.0') == 2


def test_append_different_window_raises(tmp_path):
    path = tmp_path / 'post.gmtvec'
    fid = open_postdisp_file(str(path), [2020.557, 2021.572])
    fid.close()
    with pytest.raises(ValueError, match='already contains the window'):
        open_postdisp_file(str(path), [2021.573, 2026.5])


def test_existing_file_without_window_header_raises(tmp_path):
    path = tmp_path / 'post.gmtvec'
    path.write_text(HEADER_POSTDISP + ' 1.0 2.0\n')
    with pytest.raises(ValueError, match='no tstart\\(yr\\)/tend\\(yr\\) header'):
        open_postdisp_file(str(path), [2020.557, 2021.572])


# ---------------------------------------------------------------------------
# Window resolution: eqpostdsp_tspan (displacement) vs eqpostts_tspan (series)
# ---------------------------------------------------------------------------
def test_displacement_spans_use_new_key():
    assert resolve_displacement_spans(
        [[2020.0, 2021.0], [2021.0, 2022.0]], None) == [
            [2020.0, 2021.0], [2021.0, 2022.0]]


def test_displacement_spans_fall_back_to_legacy_key():
    assert resolve_displacement_spans([], [2020.0, 2022.0]) == [[2020.0, 2022.0]]


def test_displacement_spans_prefer_new_key_and_warn_on_conflict(caplog):
    with caplog.at_level('WARNING'):
        spans = resolve_displacement_spans([2020.0, 2021.0], [2020.0, 2022.0])
    assert spans == [[2020.0, 2021.0]]
    assert 'legacy eqpost_tspan' in caplog.text


def test_ts_span_uses_explicit_window():
    assert resolve_postseismic_ts_span([2020.0, 2026.0], [[2020.0, 2021.0]]) \
        == [2020.0, 2026.0]


def test_ts_span_falls_back_to_single_displacement_window():
    assert resolve_postseismic_ts_span([], [[2020.0, 2021.0]]) == [2020.0, 2021.0]


def test_ts_span_falls_back_to_full_extent_of_multiple_windows(caplog):
    with caplog.at_level('WARNING'):
        span = resolve_postseismic_ts_span(
            [], [[2021.0, 2022.0], [2020.5, 2025.0]])
    assert span == [2020.5, 2025.0]
    assert 'eqpostts_tspan is not set' in caplog.text


def test_ts_span_empty_when_nothing_configured():
    assert resolve_postseismic_ts_span([], []) == []
    assert resolve_postseismic_ts_span('', '') == []
    assert resolve_displacement_spans('', '') == []


def test_ts_span_rejects_multiple_explicit_windows():
    with pytest.raises(ValueError, match='single \\[t0, t1\\] window'):
        resolve_postseismic_ts_span([[2020.0, 2021.0], [2021.0, 2022.0]])


# ---------------------------------------------------------------------------
# Time-series header distinguishes its window from the displacement window
# ---------------------------------------------------------------------------
class _FakeRun:
    '''Minimal stand-in for tsfitting used by output_postseismic_ts.'''

    flag = ['EQDECAY']

    def __init__(self, site, t):
        self.site = site
        self.t = np.asarray(t, dtype=float)
        self.obs = np.zeros_like(self.t)
        self.sigma = np.ones_like(self.t)

    def get_correct(self, mod_dict, eqcode=None):
        zeros = np.zeros_like(self.t)
        return zeros, zeros


def test_ts_file_header_records_its_own_window(tmp_path, monkeypatch):
    monkeypatch.chdir(tmp_path)
    run = _FakeRun('TEST', [2020.0, 2021.0, 2022.0])
    output_postseismic_ts(run, run, run, [2020.0, 2022.0])
    text = (tmp_path / 'TEST_pos.neu').read_text()
    assert '# time-series window: tstart(yr) = 2020.000' in text
    assert 'tend(yr) = 2022.000' in text
    assert 'tspan(yr) = 2.000' in text
    assert 'displacement window' not in text
