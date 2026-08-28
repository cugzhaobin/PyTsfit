'''
Deterministic input builders for the characterization tests.

Everything here writes plain text files in exactly the formats the readers in
PyTsfit expect, so the tests exercise the real parsers rather than hand-built
objects. No randomness that is not seeded, no dependence on the wall clock.
'''
import numpy as np

# Column widths of a PBO .pos data record, as consumed by posData via
# np.genfromtxt(delimiter=...). Sums to 254 characters per line.
POS_WIDTHS = (9, 7, 11, 15, 15, 15, 9, 9, 9, 7, 7, 7, 19, 16,
              11, 12, 10, 10, 11, 9, 9, 7, 7, 7, 6)

# posData hard-codes skip_header=37, so the header must be exactly this long.
POS_HEADER_LINES = 37

SITE = 'TEST'
LAT, LON, HEI = 35.5000, 103.2000, 1520.35


def _pos_header(site=SITE, lat=LAT, lon=LON, hei=HEI):
    '''
    Build a 37-line PBO header.

    posData picks the site out of ``line[16:20]`` of a line containing 'ID',
    and lat/lon/height out of ``split()[4:7]`` of the line starting with 'NEU'.
    Only one line may contain the substring 'ID', because the reader keeps
    overwriting self.site for every match.
    '''
    lines = [
        'PBO Station Position Time Series. Reference Frame : IGS14',
        'Format Version: 1.1.0',
        '4-character ID: {:4s}'.format(site),
        'Station name  : SYNTHETIC TEST STATION',
        'First Epoch   : 20100101 120000',
        'Last Epoch    : 20200101 120000',
        'Release Date  : 20200102 000000',
        'XYZ Reference position :  -1281261.00000  5640610.00000  3682825.00000 (IGS14)',
        'NEU Reference position :  {:14.9f} {:14.9f} {:10.5f} (IGS14/2010)'.format(lat, lon, hei),
        'Start Field Description',
    ]
    while len(lines) < POS_HEADER_LINES - 1:
        lines.append('Field {:d} : synthetic filler line'.format(len(lines)))
    lines.append('End Field Description')
    assert len(lines) == POS_HEADER_LINES, len(lines)
    # Guard the fragile bits of the reader contract.
    assert sum(1 for ln in lines if 'ID' in ln) == 1
    assert sum(1 for ln in lines if ln.startswith('NEU')) == 1
    return lines


def _pos_record(mjd, dn, de, du, sn, se, su):
    '''
    Format one fixed-width .pos data record.

    Displacements and sigmas go in as metres; posData converts to millimetres.
    '''
    fields = ['' for _ in POS_WIDTHS]
    fields[0]  = '20100101'
    fields[1]  = '120000'
    fields[2]  = '{:.4f}'.format(mjd)
    fields[3]  = '{:.5f}'.format(-1281261.0)
    fields[4]  = '{:.5f}'.format(5640610.0)
    fields[5]  = '{:.5f}'.format(3682825.0)
    fields[6]  = '{:.5f}'.format(0.001)
    fields[7]  = '{:.5f}'.format(0.001)
    fields[8]  = '{:.5f}'.format(0.002)
    fields[9]  = '{:.3f}'.format(0.0)
    fields[10] = '{:.3f}'.format(0.0)
    fields[11] = '{:.3f}'.format(0.0)
    fields[12] = '{:.9f}'.format(LAT)
    fields[13] = '{:.9f}'.format(LON)
    fields[14] = '{:.5f}'.format(HEI)
    fields[15] = '{:.5f}'.format(dn)
    fields[16] = '{:.5f}'.format(de)
    fields[17] = '{:.5f}'.format(du)
    fields[18] = '{:.5f}'.format(sn)
    fields[19] = '{:.5f}'.format(se)
    fields[20] = '{:.5f}'.format(su)
    fields[21] = '{:.3f}'.format(0.0)
    fields[22] = '{:.3f}'.format(0.0)
    fields[23] = '{:.3f}'.format(0.0)
    fields[24] = 'final'

    out = ''
    for value, width in zip(fields, POS_WIDTHS):
        assert len(value) <= width, (value, width)
        out += value.rjust(width)
    assert len(out) == sum(POS_WIDTHS)
    return out


def synthetic_series(n=1200, start_mjd=55197.0, step=3.0):
    '''
    Build a synthetic three-component series with a known signal:
    constant + linear trend + annual + semi-annual + a coseismic step at the
    2010-04-13 event + a logarithmic postseismic decay after it.

    Returns (mjd, N, E, U, SN, SE, SU) with displacements in metres.
    '''
    rng = np.random.default_rng(20240101)
    mjd = start_mjd + step * np.arange(n)
    # Decimal year, computed independently of GPSTime so the test signal does
    # not depend on the code under test.
    t = 2010.0 + (mjd - 55197.0) / 365.25

    t_eq = 2010.0 + (55299.0 - 55197.0) / 365.25   # matches the 'TA' event below
    step_fn = np.heaviside(t - t_eq, 0.0)
    decay = np.zeros_like(t)
    post = t > t_eq
    decay[post] = np.log(1.0 + (t[post] - t_eq) * 365.25 / 40.0)

    def component(const, vel, a_sin, a_cos, s_sin, s_cos, coseis, postamp, noise):
        return (const
                + vel * (t - t.mean())
                + a_sin * np.sin(2 * np.pi * t) + a_cos * np.cos(2 * np.pi * t)
                + s_sin * np.sin(4 * np.pi * t) + s_cos * np.cos(4 * np.pi * t)
                + coseis * step_fn
                + postamp * decay
                + rng.normal(0.0, noise, size=t.size))

    n_mm = component(12.0, -3.5, 1.8, -0.9, 0.6, 0.4, 18.0, 4.0, 1.2)
    e_mm = component(-7.0, 24.0, -1.1, 1.4, -0.5, 0.3, -26.0, -6.0, 1.4)
    u_mm = component(3.0, 0.8, 5.5, -3.2, 1.9, -1.1, 9.0, 2.5, 4.0)

    sn = np.full(t.size, 1.3)
    se = np.full(t.size, 1.5)
    su = np.full(t.size, 4.6)
    return mjd, n_mm, e_mm, u_mm, sn, se, su


def write_pos(path, site=SITE):
    '''
    Write a synthetic PBO .pos file and return the arrays that went into it.
    '''
    mjd, n_mm, e_mm, u_mm, sn, se, su = synthetic_series()
    lines = _pos_header(site=site)
    for i in range(mjd.size):
        lines.append(_pos_record(mjd[i], n_mm[i] / 1e3, e_mm[i] / 1e3, u_mm[i] / 1e3,
                                 sn[i] / 1e3, se[i] / 1e3, su[i] / 1e3))
    with open(path, 'w') as fid:
        fid.write('\n'.join(lines) + '\n')
    return mjd, n_mm, e_mm, u_mm, sn, se, su


def write_neu(path, site=SITE):
    '''
    Write a synthetic .neu file in the format neuData expects. The header line
    starting with 'NEU' carries site/lat/lon in fields 5, 6 and 7.
    '''
    mjd, n_mm, e_mm, u_mm, sn, se, su = synthetic_series()
    decyr = 2010.0 + (mjd - 55197.0) / 365.25
    with open(path, 'w') as fid:
        # neuData reads the first line containing 'NEU' and takes
        # split()[5]=site, [6]=lat, [7]=lon, so it must have >= 8 fields.
        fid.write('# NEU Reference position : {:s} {:.5f} {:.5f} {:.5f}\n'.format(site, LAT, LON, HEI))
        for i in range(decyr.size):
            fid.write('{:12.5f} {:12.4f} {:12.4f} {:12.4f} {:10.4f} {:10.4f} {:10.4f}\n'.format(
                decyr[i], n_mm[i], e_mm[i], u_mm[i], sn[i], se[i], su[i]))
    return decyr, n_mm, e_mm, u_mm, sn, se, su


# eq_rename in GAMIT/GLOBK format. 'TA' sits inside the synthetic time span and
# close enough to the station to be selected; 'FA' is far away and must be
# rejected by the distance test; 'LT' is after the series ends.
EQ_RENAME = '''*
* Synthetic eq_rename for the characterization tests
*
# near event, inside the time span
  eq_def TA 35.60 103.30 800 12 2010 4 13 12 0
  eq_rename TA
  eq_log TA 30 60
# far event, outside the 200 km radius
  eq_def FA 20.00 130.00 200 10 2012 6 15 3 30
  eq_rename FA
# late event, after the series ends
  eq_def LT 35.55 103.25 900 8 2024 1 1 0 0
  eq_rename LT
#
# non-earthquake breaks
  break TEST 2012  3 17  0  0
  break TEST 2016  9  2  0  0
# commented out, must be ignored
# break TEST 2011  1  1  0  0
  break OTHR 2013  5  5  0  0
'''


def write_eq_rename(path):
    with open(path, 'w') as fid:
        fid.write(EQ_RENAME)
    return path


def write_correction(veldir, site=SITE):
    '''
    Write velocity / offset / seasonal constraint files for class correction.
    Returns (velfile, offsetfile, periodfile).
    '''
    velfile = str(veldir / 'model.vel')
    offsetfile = str(veldir / 'offset.dat')
    periodfile = str(veldir / 'season.dat')

    # Lon Lat E N Se Sn Cne Site U Su
    with open(velfile, 'w') as fid:
        fid.write('# Lon Lat E N Se Sn Cne Site U Su\n')
        fid.write('{:10.5f} {:10.5f} {:8.2f} {:8.2f} {:6.2f} {:6.2f} {:5.2f} {:4s} {:8.2f} {:6.2f}\n'.format(
            LON, LAT, 24.0, -3.5, 0.25, 0.22, 0.0, site, 0.8, 0.55))
        # A second, distinct site: correction() iterates over velsite, which
        # crashes on a single-row 0-d array. Two rows is the realistic case.
        fid.write('{:10.5f} {:10.5f} {:8.2f} {:8.2f} {:6.2f} {:6.2f} {:5.2f} {:4s} {:8.2f} {:6.2f}\n'.format(
            LON + 0.5, LAT + 0.5, -11.0, 8.0, 0.30, 0.28, 0.0, 'OTHR', -1.2, 0.60))

    # E N U Site decyr
    with open(offsetfile, 'w') as fid:
        fid.write('# E N U Site decyr\n')
        fid.write('{:8.2f} {:8.2f} {:8.2f} {:4s} {:12.5f}\n'.format(
            -26.0, 18.0, 9.0, site, 2010.27995))
        # A second row keeps offsetdata 2-D; a single row would make
        # genfromtxt return a 1-D array and break offsetdata[indx][0, ...].
        fid.write('{:8.2f} {:8.2f} {:8.2f} {:4s} {:12.5f}\n'.format(
            0.0, 0.0, 0.0, 'OTHR', 2013.00000))

    # EAsin EAcos ESsin EScos NAsin NAcos NSsin NScos UAsin UAcos USsin UScos Site
    with open(periodfile, 'w') as fid:
        fid.write('{:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f} '
                  '{:8.2f} {:8.2f} {:8.2f} {:8.2f} {:4s}\n'.format(
                      -1.1, 1.4, -0.5, 0.3, 1.8, -0.9, 0.6, 0.4,
                      5.5, -3.2, 1.9, -1.1, site))
        fid.write('{:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f} {:8.2f} '
                  '{:8.2f} {:8.2f} {:8.2f} {:8.2f} {:4s}\n'.format(
                      0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8,
                      0.9, 1.0, 1.1, 1.2, 'OTHR'))

    return velfile, offsetfile, periodfile


def write_correction_no_offset(veldir, site=SITE):
    '''
    Same constraints as write_correction but with no offset file. This avoids
    the BREAK branch of setBoundAndInit, which today indexes offsetdata[:,4]
    (offsetdata only has three columns) and raises.
    '''
    velfile, _, periodfile = write_correction(veldir, site=site)
    return velfile, '', periodfile
