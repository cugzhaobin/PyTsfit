# ----------------------------------------------------------
# Interactive file / directory pickers for the PyTsfit Streamlit UI.
#
# Streamlit has no native file browser, so these provide:
#   * dir_picker / file_picker — a text input paired with a popover that
#     browses the local filesystem by clicking sub-directories and matching
#     files, writing the chosen path back into the input;
#   * upload_picker — materialises uploaded files into a session-scoped
#     temp directory so the rest of the app keeps treating them as local
#     files (same tsdir/eqfile logic, no path typing needed).
#
# The filesystem helpers (list_entries / materialize_uploads) are plain
# functions so they can be unit-tested without a Streamlit runtime; only the
# *_picker components talk to Streamlit.
# ----------------------------------------------------------
import os
import tempfile
import fnmatch

import streamlit as st


# --------------------------------------------------------------------------
# Pure helpers (unit-testable, no Streamlit runtime)
# --------------------------------------------------------------------------
def list_entries(path, patterns=('*',), dirs_only=False):
    '''
    List sub-directories (and matching files) under ``path``.

    Returns (dirs, files) as sorted names. Hidden entries are skipped. When
    ``dirs_only`` is True only sub-directories are returned.
    '''
    path = os.path.abspath(os.path.expanduser(path))
    if not os.path.isdir(path):
        return [], []
    names = sorted(os.listdir(path))
    dirs = [n for n in names
            if not n.startswith('.') and os.path.isdir(os.path.join(path, n))]
    if dirs_only:
        return dirs, []
    match = lambda n: any(fnmatch.fnmatch(n, p) for p in patterns)
    files = [n for n in names
             if not n.startswith('.') and os.path.isfile(os.path.join(path, n))
             and match(n)]
    return dirs, files


def materialize_uploads(uploaded_files, dest_dir):
    '''
    Write an iterable of Streamlit UploadedFile objects into ``dest_dir``.
    Returns the absolute paths written; entries with no filename are skipped.
    '''
    dest_dir = os.path.abspath(dest_dir)
    os.makedirs(dest_dir, exist_ok=True)
    paths = []
    for f in uploaded_files:
        if f is None or not f.name:
            continue
        target = os.path.join(dest_dir, os.path.basename(f.name))
        with open(target, 'wb') as fh:
            fh.write(f.getbuffer())
        paths.append(target)
    return paths


def session_dir(name):
    '''
    Return a session-scoped temp directory for a widget group, creating it
    once and reusing it for the life of the session.
    '''
    key = '_pytsfit_dir_{}'.format(name)
    d = st.session_state.get(key)
    if not d or not os.path.isdir(d):
        d = tempfile.mkdtemp(prefix='pytsfit_ui_{}_'.format(name))
        st.session_state[key] = d
    return d


# --------------------------------------------------------------------------
# Navigation callbacks (set st.session_state[key]; the popover re-renders)
# --------------------------------------------------------------------------
def _go_up(key, path):
    st.session_state[key] = os.path.dirname(os.path.abspath(path))


def _go_home(key):
    st.session_state[key] = os.path.expanduser('~')


def _pick(key, path):
    st.session_state[key] = path


def _browse(key, value, dirs_only, patterns):
    '''
    Body of the directory-browser popover. Clicking a sub-directory navigates
    deeper; clicking a file (or "select this directory") writes the chosen
    path into st.session_state[key]. The popover stays open across reruns,
    so multi-level browsing works without closing.
    '''
    path = value if os.path.isdir(value) else os.path.dirname(value)
    if not os.path.isdir(path):
        path = os.path.expanduser('~')
    path = os.path.abspath(path)

    st.caption('当前目录：`{}`'.format(path))

    c1, c2, c3 = st.columns([1, 1, 3])
    c1.button('⬆ 上级', key='{}_up'.format(key), use_container_width=True,
              on_click=_go_up, args=(key, path))
    c2.button('🏠 家', key='{}_home'.format(key), use_container_width=True,
              on_click=_go_home, args=(key,))
    if dirs_only:
        c3.button('✔ 选择此目录', key='{}_pick'.format(key), type='primary',
                  use_container_width=True, on_click=_pick, args=(key, path))

    dirs, files = list_entries(path, patterns, dirs_only=dirs_only)
    if not dirs and not files:
        st.caption('（无匹配条目）')
        return
    for d in dirs:
        st.button('📁 ' + d, key='{}_d_{}'.format(key, d),
                  use_container_width=True, on_click=_pick,
                  args=(key, os.path.join(path, d)))
    if not dirs_only:
        for f in files:
            st.button('📄 ' + f, key='{}_f_{}'.format(key, f),
                      use_container_width=True, on_click=_pick,
                      args=(key, os.path.join(path, f)))


# --------------------------------------------------------------------------
# Public widgets
# --------------------------------------------------------------------------
def dir_picker(label, key, value='', help=None):
    '''
    Directory selector: a text input paired with a popover browser. Returns
    the directory path currently held in st.session_state[key].

    The default is seeded into session_state only on first appearance; the
    widget is then created without a ``value=`` so the browser callbacks
    (which run before the script body) can rewrite the key freely.
    '''
    if key not in st.session_state:
        st.session_state[key] = value
    st.text_input(label, key=key, help=help)
    with st.popover('📂 浏览…'):
        _browse(key, st.session_state[key], dirs_only=True, patterns=())
    return st.session_state[key]


def file_picker(label, key, value='', patterns=('*',), help=None):
    '''
    File selector: a text input paired with a popover browser filtered by
    ``patterns`` (fnmatch globs). Returns the chosen path.
    '''
    if key not in st.session_state:
        st.session_state[key] = value
    st.text_input(label, key=key, help=help)
    with st.popover('📂 浏览…'):
        _browse(key, st.session_state[key], dirs_only=False, patterns=patterns)
    return st.session_state[key]


def upload_picker(label, key, dest_dir, accepted=(), multiple=True):
    '''
    Upload widget that writes accepted files into ``dest_dir`` (session
    scoped) and returns ``dest_dir`` once at least one file is present,
    otherwise None.
    '''
    up = st.file_uploader(label, key=key, type=list(accepted) or None,
                          accept_multiple_files=multiple)
    if not up:
        return None
    files = up if multiple else [up]
    materialize_uploads([f for f in files if f is not None], dest_dir)
    return dest_dir
