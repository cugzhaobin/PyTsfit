#!/usr/bin/env python
# ----------------------------------------------------------
# PyTsfit interactive UI (Streamlit).
#
# Launch with:
#     streamlit run src/pytsfit/ui/app.py
#
# Layout — the sidebar is a five-step console that drives a single-page
# pipeline on the right:
#
#   ① 数据源      选数据（本地目录/上传文件 + 目录浏览器）
#   ② 模型参数
#   ③ 拟合选项
#   ④ 站点与运行  ▶ 运行拟合（当前站点 / 批量全部）
#   ⑤ 绘图/输出
#
# All state lives in st.session_state: changing the data source invalidates
# the loaded site list, changing model/fit options invalidates the last run.
# The right side is display-only and reads that state, so the sidebar's
# actions are always reflected immediately.
#
# Thin shell only: all data loading / fitting / plotting is delegated to
# the existing pytsfit modules through pytsfit.ui.api and the pickers in
# pytsfit.ui.widgets.
# ----------------------------------------------------------
import os
import logging

import yaml
import streamlit as st

from pytsfit.qualitycontrol import DEFAULT_FIT_OPTS
from pytsfit.ui import api
from pytsfit.ui import widgets

st.set_page_config(page_title='PyTsfit UI', page_icon='📈', layout='wide')

logging.basicConfig(level=logging.WARNING)

_DEFAULT_CFG_PATH = os.path.join(os.path.dirname(__file__),
                                 '..', 'scripts', 'config.yaml')

_POS_NEU = ('pos', 'neu')

# Sidebar widgets whose keys we clear on "恢复默认".
_CONFIG_WIDGET_KEYS = ('tsdir', 'sitefile', 'eqfile', 'eq_upload',
                       'prior_velfile', 'prior_offsetfile', 'prior_periodfile',
                       'src_mode', 'tspan0', 'tspan1')


# --------------------------------------------------------------------------
# Defaults (mirror the CLI's config.yaml so the UI starts from the same
# settings the command line uses).
# --------------------------------------------------------------------------
def load_default_cfg():
    cfg = {
        'dict_input': {'eqfile': './eq_rename', 'sitefile': '', 'tsdir': '../pos/',
                       'tsformat': 'pos', 'timespan': [1998, 2024],
                       'prior_velfile': '', 'prior_offsetfile': '', 'prior_periodfile': ''},
        'dict_param': {'constant': True, 'linear': True, 'annual': True, 'semiannual': True,
                       'break': True, 'eqoffset_ne': True, 'eqoffset_up': True,
                       'eqpost_ne': False, 'eqpost_up': False},
        'dict_plot': {'detrend': True, 'debreak': True, 'deeqoffset': False,
                      'depost': False, 'deseason': False},
        'dict_fit': dict(DEFAULT_FIT_OPTS),
        'dict_output': {'eqoffset': '', 'break': '', 'eqpostdisp': '',
                        'eqpostts': False, 'eqpost_tspan': [2015, 2023]},
    }
    try:
        with open(_DEFAULT_CFG_PATH) as fid:
            loaded = yaml.load(fid, Loader=yaml.FullLoader)
        for section in cfg:
            if section in loaded:
                cfg[section].update({k: v for k, v in loaded[section].items()
                                     if k in cfg[section]})
    except FileNotFoundError:
        pass
    return cfg


def _input_fingerprint(di):
    '''Fingerprint of the inputs that determine the site list.'''
    return (di.get('tsdir'), di.get('tsformat'), di.get('sitefile'))


def _fit_fingerprint(dp, fo):
    '''Fingerprint of the model + fit options, for stale-run detection.'''
    return (repr(sorted(dp.items())), repr(sorted(fo.items())))


# Cached data loaders. The cache key includes the input fingerprint so a
# changed data source can never serve a stale site's series.
@st.cache_data(show_spinner=False)
def _cached_load_site(tsdir, tsformat, site, token):
    return api.load_site(tsdir, tsformat, site)


@st.cache_data(show_spinner=False)
def _cached_station_metadata(tsdir, tsformat, sites, token):
    return api.station_metadata(tsdir, tsformat, tuple(sites))


def _inputs_valid(di, dp):
    if not os.path.isdir(di['tsdir']):
        st.error('数据目录不存在：{}'.format(di['tsdir']))
        return False
    need_eq = (dp['eqoffset_ne'] or dp['eqoffset_up']
               or dp['eqpost_ne'] or dp['eqpost_up'] or dp['break'])
    if need_eq and not os.path.isfile(di['eqfile']):
        st.error('已启用同震/间断/震后模型，但 eqfile 不存在：{}'.format(di['eqfile']))
        return False
    return True


# --------------------------------------------------------------------------
# Sidebar — the five-step console
# --------------------------------------------------------------------------
def sidebar(cfg):
    with st.sidebar:
        st.title('⚙️ PyTsfit 控制台')
        di = _section_data_source(cfg)
        dp = _section_model(cfg)
        fo = _section_fit_options(cfg)
        _section_run(di, dp, fo)
        dplt, dout = _section_plot_output(cfg)
        st.divider()
        c1, c2 = st.columns(2)
        c1.button('恢复默认', key='btn_reset_cfg', use_container_width=True,
                  on_click=_reset_config)
        c2.button('清空结果', key='btn_reset_data', use_container_width=True,
                  on_click=_reset_data)
    return di, dp, fo, dplt, dout


def _section_data_source(cfg):
    di = dict(cfg['dict_input'])
    with st.expander('① 数据源', expanded=True):
        source = st.radio('数据来源', ('本地目录', '上传文件'),
                          key='src_mode', horizontal=True)

        if source == '本地目录':
            di['tsdir'] = widgets.dir_picker(
                '数据目录 tsdir', key='tsdir', value=di['tsdir'],
                help='存放 *.pos / *.neu 时间序列的目录')
            fmt_index = _POS_NEU.index(di['tsformat']) \
                if di['tsformat'] in _POS_NEU else 0
            di['tsformat'] = st.selectbox('格式 tsformat', _POS_NEU,
                                          index=fmt_index)
            di['sitefile'] = widgets.file_picker(
                '站点列表 sitefile（可选）', key='sitefile', value=di['sitefile'])
        else:
            up_dir = widgets.session_dir('ts')
            di['tsdir'] = up_dir
            st.caption('上传的序列写入临时目录：`{}`'.format(up_dir))
            widgets.upload_picker('上传 *.pos / *.neu 序列（可多选）',
                                  key='ts_upload', dest_dir=up_dir,
                                  accepted=_POS_NEU)
            fmt_index = _POS_NEU.index(di['tsformat']) \
                if di['tsformat'] in _POS_NEU else 0
            di['tsformat'] = st.selectbox('格式 tsformat', _POS_NEU,
                                          index=fmt_index)
            di['sitefile'] = ''

        di['eqfile'] = widgets.file_picker(
            '地震目录 eqfile', key='eqfile', value=di['eqfile'],
            patterns=('eq_rename*', 'eq*', '*'))
        with st.popover('⬆ 上传 eq_rename'):
            up_eq = widgets.session_dir('eq')
            paths = widgets.upload_picker(
                '选择要上传的 eq_rename', key='eq_upload', dest_dir=up_eq,
                accepted=('eq_rename', 'eq*'), multiple=False)
            if paths and st.session_state.get('eqfile') != paths[0]:
                st.button('使用此文件', key='eq_use',
                          on_click=_set_cfg_key, args=('eqfile', paths[0]))

        with st.expander('先验文件（可选）'):
            di['prior_velfile'] = widgets.file_picker(
                '先验速度文件', key='prior_velfile', value=di['prior_velfile'])
            di['prior_offsetfile'] = widgets.file_picker(
                '先验偏移文件', key='prior_offsetfile', value=di['prior_offsetfile'])
            di['prior_periodfile'] = widgets.file_picker(
                '先验周期文件', key='prior_periodfile', value=di['prior_periodfile'])

        t0, t1 = di['timespan']
        c1, c2 = st.columns(2)
        di['timespan'] = [c1.number_input('起始年', value=float(t0),
                                          format='%.1f', key='tspan0'),
                          c2.number_input('结束年', value=float(t1),
                                          format='%.1f', key='tspan1')]

        # Status + load action
        loaded = st.session_state.get('loaded_cfg')
        sites = st.session_state.get('sites', [])
        current = _input_fingerprint(di)
        if loaded is not None and loaded == current:
            st.caption('● 已加载 **{}** 个站点'.format(len(sites)))
        else:
            st.caption('○ 配置已变更，需重新加载站点')
        if st.button('🔄 加载站点',
                     type='primary' if loaded != current else 'secondary',
                     key='btn_load_sites', use_container_width=True):
            _load_sites(di)
            st.rerun()
    return di


def _section_model(cfg):
    dp = {}
    p = cfg['dict_param']
    with st.expander('② 模型参数', expanded=True):
        dp['constant'] = st.checkbox('常数项', value=p['constant'])
        dp['linear'] = st.checkbox('线性速度', value=p['linear'])
        dp['annual'] = st.checkbox('周年周期', value=p['annual'])
        dp['semiannual'] = st.checkbox('半周年周期', value=p['semiannual'])
        dp['break'] = st.checkbox('非地震间断 (break)', value=p['break'])
        st.caption('同震偏移')
        dp['eqoffset_ne'] = st.checkbox('水平 NE', value=p['eqoffset_ne'])
        dp['eqoffset_up'] = st.checkbox('垂直 Up', value=p['eqoffset_up'])
        st.caption('震后形变')
        dp['eqpost_ne'] = st.checkbox('水平 NE', value=p['eqpost_ne'])
        dp['eqpost_up'] = st.checkbox('垂直 Up', value=p['eqpost_up'])
    return dp


def _section_fit_options(cfg):
    fo = dict(cfg['dict_fit'])
    with st.expander('③ 拟合选项', expanded=False):
        fo['sigma_scale'] = st.selectbox(
            'sigma_scale 误差缩放', ('none', 'nrms', 'realistic'),
            index=('none', 'nrms', 'realistic').index(fo.get('sigma_scale', 'nrms')))
        fo['outlier'] = st.checkbox('迭代离群点剔除 (outlier)',
                                    value=fo.get('outlier', False))
        fo['nsigma'] = st.number_input('nsigma（剔除阈值）',
                                       value=float(fo.get('nsigma', 4.0)),
                                       min_value=0.0, step=0.5)
        fo['outlier_scale'] = st.selectbox(
            'outlier_scale', ('mad', 'sigma'),
            index=('mad', 'sigma').index(fo.get('outlier_scale', 'mad')))
        fo['max_iter'] = st.number_input('max_iter',
                                         value=int(fo.get('max_iter', 10)),
                                         min_value=1, step=1)
        fo['restore_factor'] = st.number_input(
            'restore_factor', value=float(fo.get('restore_factor', 0.9)),
            min_value=0.0, max_value=1.0, step=0.05)
        fo['max_sigma'] = st.number_input('max_sigma（mm，0=关）',
                                          value=float(fo.get('max_sigma') or 0.0),
                                          min_value=0.0, step=1.0)
        c1, c2 = st.columns(2)
        fo['min_rsig'] = c1.number_input('min_rsig',
                                         value=int(fo.get('min_rsig', 30)),
                                         min_value=1, step=1)
        fo['min_sigscale'] = c2.number_input('min_sigscale',
                                             value=int(fo.get('min_sigscale', 10)),
                                             min_value=1, step=1)
    return fo


def _section_run(di, dp, fo):
    with st.expander('④ 站点与运行', expanded=True):
        sites = st.session_state.get('sites', [])
        if not sites:
            st.caption('请先在 ① 加载站点。')
            return
        mode = st.radio('运行模式', ('当前站点', '批量全部'),
                        key='run_mode', horizontal=True)
        site = st.session_state.get('selected_site')
        if mode == '当前站点':
            site = _site_selector(sites)
        else:
            st.caption('将批量拟合全部 **{}** 个站点。'.format(len(sites)))
        if st.button('▶ 运行拟合', type='primary', key='btn_run',
                     use_container_width=True):
            _run_fit(di, dp, fo, mode, site, sites)
        _run_status(dp, fo)


def _run_status(dp, fo):
    fr = st.session_state.get('fit_result')
    if fr is None:
        st.caption('尚未运行拟合。')
        return
    if _fit_fingerprint(dp, fo) != st.session_state.get('fit_cfg'):
        st.caption('○ 模型/拟合选项已变更，需重新运行')
    if fr['mode'] == 'single':
        try:
            summary = api.fit_summary(fr['runs'])
            wrms = {r['Component']: r['WRMS (mm)'] for _, r in summary.iterrows()}
            txt = 'N {:.2f} | E {:.2f} | U {:.2f}'.format(
                wrms.get('N', float('nan')), wrms.get('E', float('nan')),
                wrms.get('U', float('nan')))
        except Exception:
            txt = ''
        st.caption('✓ 单站 **{}**：WRMS {}'.format(fr['site'], txt))
    else:
        n_ok, n_fail = len(fr['results']), len(fr['failures'])
        st.caption('✓ 批量：成功 **{}/{}**'.format(n_ok, n_ok + n_fail))
        if n_fail:
            st.warning('失败站点：{}'.format(', '.join(s for s, _ in fr['failures'])))


def _section_plot_output(cfg):
    dplt = dict(cfg['dict_plot'])
    dout = dict(cfg['dict_output'])
    with st.expander('⑤ 绘图 / 输出', expanded=False):
        st.caption('obs/mod 图去趋势项')
        dplt['detrend'] = st.checkbox('去趋势 detrend', value=dplt.get('detrend', True))
        dplt['debreak'] = st.checkbox('去间断 debreak', value=dplt.get('debreak', True))
        dplt['deeqoffset'] = st.checkbox('去同震偏移', value=dplt.get('deeqoffset', False))
        dplt['depost'] = st.checkbox('去震后', value=dplt.get('depost', False))
        dplt['deseason'] = st.checkbox('去季节', value=dplt.get('deseason', False))
        st.divider()
        dout['velfile'] = st.checkbox('速度表 velfile', value=True)
        dout['eqoffset'] = st.checkbox('同震偏移表 eqoffset',
                                       value=bool(cfg['dict_output'].get('eqoffset')))
        dout['break'] = st.checkbox('间断表 break',
                                    value=bool(cfg['dict_output'].get('break')))
        dout['eqpostdisp'] = st.checkbox('震后位移表 eqpostdisp',
                                         value=bool(cfg['dict_output'].get('eqpostdisp')))
        dout['eqpostts'] = st.checkbox('震后时间序列 eqpostts',
                                       value=bool(cfg['dict_output'].get('eqpostts', False)))
        t0, t1 = cfg['dict_output'].get('eqpost_tspan', [2015, 2023])
        c1, c2 = st.columns(2)
        dout['eqpost_tspan'] = [c1.number_input('震后起年', value=float(t0), format='%.1f'),
                                c2.number_input('震后止年', value=float(t1), format='%.1f')]
        st.caption('单站输出下载见右侧 ④ 参数表。')
    return dplt, dout


# --------------------------------------------------------------------------
# Sidebar actions
# --------------------------------------------------------------------------
def _set_cfg_key(key, value):
    st.session_state[key] = value


def _reset_config():
    _reset_data()
    for k in _CONFIG_WIDGET_KEYS:
        st.session_state.pop(k, None)


def _reset_data():
    for k in ('sites', 'loaded_cfg', 'fit_result', 'fit_cfg', 'selected_site'):
        st.session_state.pop(k, None)


def _load_sites(di):
    with st.spinner('正在扫描站点…'):
        sites = api.list_sites(di['tsdir'], di['tsformat'], di['sitefile'])
    st.session_state['sites'] = sites
    st.session_state['loaded_cfg'] = _input_fingerprint(di)
    current = st.session_state.get('selected_site')
    if current not in sites:
        st.session_state['selected_site'] = sites[0] if sites else None
    st.session_state.pop('fit_result', None)
    st.session_state.pop('fit_cfg', None)


def _site_selector(sites):
    current = st.session_state.get('selected_site')
    if current not in sites:
        current = sites[0]
    index = sites.index(current) if current in sites else 0
    chosen = st.selectbox('当前站点', sites, index=index)
    st.session_state['selected_site'] = chosen
    return chosen


def _run_fit(di, dp, fo, mode, site, sites):
    st.session_state.pop('fit_result', None)
    try:
        if mode == '当前站点':
            if not site:
                st.warning('请先选择站点。')
                return
            with st.spinner('正在拟合 {} 的 N/E/U 分量…'.format(site)):
                res = api.run_fit(di, dp, fo, site)
            st.session_state['fit_result'] = {
                'mode': 'single', 'site': site,
                'data': res['data'], 'runs': res['runs'], 'failures': [],
            }
        else:
            prog = st.progress(0.0, text='批量拟合 0/{}'.format(len(sites)))

            def _cb(i, n):
                prog.progress(i / n, text='批量拟合 {}/{}'.format(i, n))

            batch = api.run_fit_batch(di, dp, fo, sites, progress_cb=_cb)
            prog.progress(1.0, text='批量拟合完成')
            st.session_state['fit_result'] = {
                'mode': 'batch', 'site': site,
                'results': batch['results'], 'failures': batch['failures'],
            }
    except Exception as exc:
        st.error('拟合失败：{}'.format(exc))
        return
    st.session_state['fit_cfg'] = _fit_fingerprint(dp, fo)


# --------------------------------------------------------------------------
# Main page — single pipeline driven by session state
# --------------------------------------------------------------------------
def main():
    cfg = load_default_cfg()
    di, dp, fo, dplt, dout = sidebar(cfg)

    st.title('PyTsfit — GNSS 时间序列拟合')
    st.caption('左侧控制台完成 ①选数据 → ②配模型 → ④运行；本页按步骤展示结果。')

    if not _inputs_valid(di, dp):
        st.stop()

    _section_site_map(di)
    _section_raw_ts(di)
    _section_fit_results(dplt)
    _section_param_table(di, dp, dout)


def _section_site_map(di):
    st.subheader('① 站点分布')
    sites = st.session_state.get('sites', [])
    if not sites:
        st.info('未发现站点。请检查 ① 数据源 后点击"加载站点"。')
        return
    current = st.session_state.get('selected_site')
    if current in sites:
        st.caption('当前站点：**{}**（也可在地图上点击切换）'.format(current))
    token = st.session_state.get('loaded_cfg', ())
    meta = _cached_station_metadata(di['tsdir'], di['tsformat'],
                                    tuple(sites), token)
    fig = api.make_sitemap(meta)
    try:
        event = st.plotly_chart(fig, on_select='rerun', selection_mode='points')
        if event and event.selection and event.selection.points:
            pts = event.selection.points
            cd = pts[0].get('customdata') if pts else None
            if cd:
                cd = cd[0] if isinstance(cd, (list, tuple)) else cd
                if str(cd) in sites:
                    st.session_state['selected_site'] = str(cd)
                    st.rerun()
    except TypeError:
        st.plotly_chart(fig)


def _section_raw_ts(di):
    st.subheader('② 原始时间序列')
    sites = st.session_state.get('sites', [])
    site = st.session_state.get('selected_site')
    if not site or site not in sites:
        st.info('先在左侧 ① 加载站点。')
        return
    st.caption('站点：**{}**'.format(site))
    token = st.session_state.get('loaded_cfg', ())
    try:
        data = _cached_load_site(di['tsdir'], di['tsformat'], site, token)
    except Exception as exc:
        st.error(str(exc))
        return
    st.plotly_chart(api.make_raw_ts_figure(data), width='stretch')


def _section_fit_results(dplt):
    st.subheader('③ 拟合结果')
    fr = st.session_state.get('fit_result')
    if fr is None:
        st.info('尚未运行拟合。在左侧 ④ 点击"运行拟合"。')
        return
    if fr['mode'] == 'batch':
        _render_batch_results(fr, dplt)
    else:
        _render_single_results(fr, dplt)


def _render_single_results(fr, dplt):
    runs = fr['runs']
    if all(not hasattr(runs[c], 'param') or len(runs[c].param) == 0
           for c in api.COMPONENTS):
        st.warning('所选模型在该站点未拟合出任何参数（请检查模型开关）。')
        return
    st.plotly_chart(api.make_obsmod_figure(runs, dplt), width='stretch')
    st.plotly_chart(api.make_residual_figure(runs), width='stretch')


def _render_batch_results(fr, dplt):
    results = fr['results']
    if fr['failures']:
        st.warning('{} 个站点拟合失败：{}'.format(
            len(fr['failures']), ', '.join(s for s, _ in fr['failures'])))
    if not results:
        st.warning('批量拟合无成功站点。')
        return
    st.markdown('**查看单个站点拟合图**')
    detail_site = st.selectbox('选择站点', sorted(results), index=0)
    _render_single_results({'runs': results[detail_site]['runs']}, dplt)


def _section_param_table(di, dp, dout):
    st.subheader('④ 参数表与导出')
    fr = st.session_state.get('fit_result')
    if fr is None:
        st.info('运行拟合后显示参数表。')
        return
    if fr['mode'] == 'batch':
        _render_batch_params(fr)
    else:
        _render_single_params(fr, dp, dout)


def _render_single_params(fr, dp, dout):
    runs = fr['runs']
    table = api.param_table(runs)
    if table.empty:
        st.warning('没有可显示的参数。')
    else:
        st.dataframe(table, width='stretch', hide_index=True)
    st.markdown('**拟合质量**')
    st.dataframe(api.fit_summary(runs), width='stretch', hide_index=True)

    st.markdown('**输出下载（与 do_pytsfit 相同的输出格式）**')
    try:
        csvs = api.outputs_to_csv(runs, dp, dout)
    except Exception as exc:
        st.caption('导出失败：{}'.format(exc))
        return
    labels = {
        'velfile': ('速度表 velfile', '{}.vel'.format(fr['site'])),
        'eqoffset': ('同震偏移表', '{}.eqoffset'.format(fr['site'])),
        'break': ('间断表 break', '{}.break'.format(fr['site'])),
        'eqpostdisp': ('震后位移表', '{}.postdisp'.format(fr['site'])),
        'eqpostts': ('震后时间序列', '{}.postts'.format(fr['site'])),
    }
    for key, (label, fname) in labels.items():
        text = csvs.get(key, '')
        if text.strip():
            st.download_button(label, data=text, file_name=fname,
                               mime='text/plain')
        else:
            st.caption('{}：该模型项未估计，无输出。'.format(label))


def _render_batch_params(fr):
    results = fr['results']
    if not results:
        st.warning('批量拟合无成功站点。')
        return
    st.markdown('**全站速度汇总（mm/yr）**')
    table = api.batch_velocity_table(results)
    st.dataframe(table, width='stretch')
    text = api.batch_output_velo(results, fmt='DETAIL')
    if text.strip():
        st.download_button('⬇ 下载全部速度表（velfile 格式）',
                           data=text, file_name='batch.vel',
                           mime='text/plain')
    st.caption('各站点参数表：在 ③ 选择站点查看单站拟合图；单站下载按钮在单站模式下使用。')


if __name__ == '__main__':
    main()
