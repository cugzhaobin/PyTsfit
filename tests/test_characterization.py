'''
Characterization test: regenerate the full behaviour snapshot and require it to
match the committed golden baseline exactly.

This is the regression net for the refactor. Run it before and after moving
code; a byte-for-byte match of the JSON is the definition of "no behaviour
change". After an intentional behaviour change (e.g. the bug fixes), regenerate
the baseline with:

    python tests/snapshot.py tests/golden/baseline.json

and review the diff of that JSON before committing it.
'''
import json
import math
import pathlib

import pytest

import snapshot


def _leaf_equal(a, b):
    '''Compare two JSON leaf values, treating NaN == NaN.'''
    if isinstance(a, float) and isinstance(b, float):
        if math.isnan(a) and math.isnan(b):
            return True
    return a == b

HERE = pathlib.Path(__file__).resolve().parent
GOLDEN = HERE / 'golden' / 'baseline.json'


@pytest.fixture(scope='module')
def golden():
    with open(GOLDEN) as fid:
        return json.load(fid)


@pytest.fixture(scope='module')
def current(tmp_path_factory):
    td = tmp_path_factory.mktemp('snapshot')
    return snapshot.build_snapshot(td)


def _diff_paths(a, b, prefix=()):
    '''Yield (dotted.path, a_value, b_value) for every leaf where a != b.'''
    if isinstance(a, dict) and isinstance(b, dict):
        for key in a.keys() | b.keys():
            if key not in a:
                yield (prefix + (key,), '<missing>', b[key])
            elif key not in b:
                yield (prefix + (key,), a[key], '<missing>')
            else:
                yield from _diff_paths(a[key], b[key], prefix + (key,))
    elif isinstance(a, list) and isinstance(b, list):
        if len(a) != len(b):
            yield (prefix, 'len %d' % len(a), 'len %d' % len(b))
        for i, (x, y) in enumerate(zip(a, b)):
            yield from _diff_paths(x, y, prefix + (str(i),))
    else:
        if not _leaf_equal(a, b):
            yield (prefix, a, b)


def test_snapshot_matches_golden(golden, current):
    diffs = list(_diff_paths(golden, current))
    assert not diffs, (
        'Snapshot diverged from golden baseline. First differences:\n' +
        '\n'.join('  {}: {!r}  !=  {!r}'.format('.'.join(p), a, b)
                  for p, a, b in diffs[:20]) +
        ('\n  ... (%d total)' % len(diffs) if len(diffs) > 20 else '') +
        '\nRegenerate with: python tests/snapshot.py tests/golden/baseline.json'
    )
