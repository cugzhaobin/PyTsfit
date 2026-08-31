'''
Tests for the filesystem helpers in pytsfit.ui.widgets.

The Streamlit components themselves (dir_picker / file_picker /
upload_picker) need a running Streamlit runtime, so only the pure
filesystem helpers are exercised here: directory listing / pattern matching
and materialising uploaded files to disk.
'''
import os

from pytsfit.ui.widgets import list_entries, materialize_uploads


class _FakeUpload:
    '''Minimal stand-in for streamlit.runtime.uploaded_file_manager.UploadedFile.'''
    def __init__(self, name, data=b''):
        self.name = name
        self._data = data

    def getbuffer(self):
        return self._data


def _tree(tmp_path):
    (tmp_path / 'a_dir').mkdir()
    (tmp_path / 'b_dir').mkdir()
    (tmp_path / 'hidden_dir').mkdir()
    (tmp_path / 'aa.pos').write_text('x')
    (tmp_path / 'bb.neu').write_text('x')
    (tmp_path / 'cc.dat').write_text('x')
    (tmp_path / '.hidden.pos').write_text('x')
    return tmp_path


def test_list_entries_sorted_and_hidden_skipped(tmp_path):
    d = _tree(tmp_path)
    dirs, files = list_entries(str(d))
    assert dirs == ['a_dir', 'b_dir']
    assert files == ['aa.pos', 'bb.neu', 'cc.dat']


def test_list_entries_pattern_filtering(tmp_path):
    d = _tree(tmp_path)
    _, files = list_entries(str(d), patterns=('*.pos',))
    assert files == ['aa.pos']
    _, files = list_entries(str(d), patterns=('eq_rename*', 'eq*', '*'))
    assert files == ['aa.pos', 'bb.neu', 'cc.dat']


def test_list_entries_dirs_only(tmp_path):
    d = _tree(tmp_path)
    dirs, files = list_entries(str(d), dirs_only=True)
    assert dirs == ['a_dir', 'b_dir']
    assert files == []


def test_list_entries_missing_path(tmp_path):
    assert list_entries(str(tmp_path / 'nope')) == ([], [])


def test_materialize_uploads_writes_files(tmp_path):
    dest = str(tmp_path / 'out')
    uploads = [_FakeUpload('AAA.pos', b'data-a'), _FakeUpload('BBB.pos', b'data-b'),
               _FakeUpload('', b'ignored')]
    paths = materialize_uploads(uploads, dest)
    assert sorted(os.path.basename(p) for p in paths) == ['AAA.pos', 'BBB.pos']
    for p in paths:
        assert os.path.isfile(p)
    assert open(os.path.join(dest, 'AAA.pos'), 'rb').read() == b'data-a'


def test_materialize_uploads_creates_dest(tmp_path):
    dest = str(tmp_path / 'a' / 'b')
    paths = materialize_uploads([_FakeUpload('X.pos', b'x')], dest)
    assert os.path.isfile(paths[0])
