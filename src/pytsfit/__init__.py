'''
pytsfit: fitting position time series of GNSS data.

The module layout follows responsibilities:

    data       -- file readers (neuData, posData)
    models     -- events and priors (earthquake, eqcatalog, eqPost, ...)
    tsfitting  -- the fitting engine
    output     -- output_* writers and plotting
    PyTsfit    -- compatibility shim re-exporting the above

``geotools`` and ``GPSTime`` are bundled here and imported by the fitting code
via relative imports, so the package no longer depends on same-named modules on
sys.path.
'''

from . import geotools as gt
from . import GPSTime as gpstime

from .data import neuData, posData
from .models import (earthquake, eqcatalog, eqPost, eqPostList,
                     offset, breakcatalog, correction, build_param_dict)
from .tsfitting import tsfitting
from . import qualitycontrol
from .qualitycontrol import (DEFAULT_FIT_OPTS, merge_fit_opts, robust_scale,
                             realistic_sigma, flag_outliers)
from .output import (output_velo, output_postseismic_velo, output_postseismic_disp,
                     output_eqoffset, output_break, output_postseismic_ts,
                     output_period, output_summary, plot_obs_mod,
                     output_obs_mod, output_param)

__all__ = [
    # geotools / GPSTime (kept importable for backward compatibility)
    'gt', 'gpstime',
    # data
    'neuData', 'posData',
    # models
    'earthquake', 'eqcatalog', 'eqPost', 'eqPostList',
    'offset', 'breakcatalog', 'correction', 'build_param_dict',
    # fitting engine
    'tsfitting',
    # quality control
    'qualitycontrol', 'DEFAULT_FIT_OPTS', 'merge_fit_opts', 'robust_scale',
    'realistic_sigma', 'flag_outliers',
    # output / plotting
    'output_velo', 'output_postseismic_velo', 'output_postseismic_disp',
    'output_eqoffset', 'output_break', 'output_postseismic_ts',
    'output_period', 'output_summary', 'plot_obs_mod',
    'output_obs_mod', 'output_param',
]
