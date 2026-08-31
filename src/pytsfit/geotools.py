# Written by Zhao Bin, Insitutute of Seismology, CEA. Jan 23th, 2016
# Mod by Zhao Bin, June 1, 2016 for xy2llh function
import numpy as np
from numpy import sin, cos, sqrt, deg2rad, rad2deg
from numpy import linspace, array, zeros, savetxt
from pyproj import Proj

__all__ = ['llh2localxy','polyconic', 'xy2llh', 'vpvs2nu', 'vpvsrho2K', 'vsrho2mu', 'makegrids']


def distance3d(pt1, pt2):
    '''
    compute 3D distance
    Input: 
        pt1  -- [lon, lat, dep], unit: degree, degree, km
        pt2  -- [lon, lat, dep], unit: degree, degree, km
    Output:
        distance in kilometer
    '''
    [e,n]  = llh2localxy([pt1[1], pt1[0]], [pt2[1], pt2[0]]) 
    u      = pt1[2]-pt2[2]
    dist3d = np.sqrt(e**2 + n**2 + u**2) 
    return dist3d


def llh2utm(llh, llh_org):
    '''
    convert coordinate from [lat, lon] to [east, north] w.r.t original coordinate
    Written by Zhao Bin, Institute of Seismology, CEA. Aug 8, 2017
    This function need pyproj 
    
    Input:
      llh    - inout latitude and longitude in degree
      ll_org - refernence point (lat, lon)
    Output:
      xy     - east(km), north(km)
     
    '''
    P      = Proj(ellps='WGS84', proj='tmerc', lon_0=llh_org[1])
    en     = P(llh[1], llh[0])
    en_org = P(llh_org[1], llh_org[0])
    e_loc  = (en[0]-en_org[0])/1000.0
    n_loc  = (en[1]-en_org[1])/1000.0

    return np.array([e_loc, n_loc])


def utm2llh(xy, origin):
    '''
    convert coordinate from [east, north] to [lat, lon] w.r.t original coordinate
    Written by Zhao Bin, Institute of Seismology, CEA. Aug 8, 2017
    This function need pyproj 
    
    Input:
      xy     - xy[:,0] = north, xy[:,1] = east, unit is kilometer
      origin - [lat, lon]
    Output:
      llh    - llh[:,0] = lat, llh[:,1] = lon
    '''
    llh_org = origin
    P      = Proj(ellps='WGS84', proj='tmerc', lon_0=llh_org[1])
    en_org = P(llh_org[1], llh_org[0])

    xy = np.mat(xy)
    [m, n] = xy.shape
    # Now we convert the matrix to array
    xy = np.array(xy)
    llh = np.zeros((m, 2))

    for i in range(len(xy)):
        n_abs  = xy[i,0]*1000.0 + en_org[1]
        e_abs  = xy[i,1]*1000.0 + en_org[0]
        new    = P(e_abs, n_abs, inverse=True)
        llh[:,0] = new[1]
        llh[:,1] = new[0]

    return llh


def llh2localxy(llh, ll_org):
    '''
    convert coordinate from lat lon to XY
    
    Input:
      llh    - inout latitude and longitude in degree
      ll_org - refernence point (lat, lon)
    Output:
      xy     - east(km), north(km)
     
    '''

 
    lat = 3600*llh[0]
    lon = 3600*llh[1]
    
    Lat_Orig  = 3600.0 * ll_org[0]
    Diff_long = 3600.0 * ll_org[1] - lon
    xy = polyconic(lat, Diff_long, Lat_Orig)

    xy[0] = -xy[0]/1000.0
    xy[1] =  xy[1]/1000.0

    return xy


def llh2local(llh, origin):
    '''
    Converts from longitude and latitude to local coorindates given an origin.  llh and origin should be in decimal degrees.
    Written by Yi Qujie, Institute of Seismology. 2017
    Modified by Zhao Bin, Institute of Seismology, CEA. April 20, 2018
    
    INPUT:
        llh    = np.array([lat, lon, up])
        origin = np.array([lat, lon])
    OUTPUT:
        xyz    = np.array([north, east, up])
    '''

    # Set ellipsoid constants (WGS84)
    a = 6378137.0
    e = 0.08209443794970

    # Convert to radians
    llh    = np.dot(llh, np.diag([np.deg2rad(1), np.deg2rad(1), 1]))
    origin = np.deg2rad(origin)

    # Do the projection
    xyz    = np.zeros(llh.shape)
    z      = (llh[:,0] !=0 )
    for i in range(0, llh.shape[0]):
        dlambda = llh[i,1] - origin[1]
        M       = a*((1-e**2/4-3*e**4/64-5*e**6/256)*llh[i,0] - (3*e**2/8+3*e**4/32+45*e**6/1024)*np.sin(2*llh[i,0]) + (15*e**4/256 +45*e**6/1024)*np.sin(4*llh[i,0]) - (35*e**6/3072)*np.sin(6*llh[i,0]))
        M0      = a*((1-e**2/4-3*e**4/64-5*e**6/256)*origin[0] - (3*e**2/8+3*e**4/32+45*e**6/1024)*np.sin(2*origin[0]) + (15*e**4/256 +45*e**6/1024)*np.sin(4*origin[0]) - (35*e**6/3072)*np.sin(6*origin[0]))
        if z[i]:
           N    = a/np.sqrt(1-e**2*np.sin(llh[i,0])**2)
           E    = dlambda*np.sin(llh[i,0])
           xyz[i,0] = N*np.sin(E)/np.tan(llh[i,0])
           xyz[i,1] = M - M0 + N*(1-np.cos(E))/np.tan(llh[i,0])
       # Handle special case of latitude = 0
        else:
           xyz[i,0] = a*dlambda
           xyz[i,1] = -M0
    xyz[:,2] = llh[:,2]
    # Convert to km
    xyz = np.dot(xyz, np.diag([0.001,0.001,1]))
    return xyz

def local2llh(xyh, origin):
    '''
    Converts from local coorindates to longitude and latitude given the [ lat,lon] of an origin. 'origin' should be in 
    decimal degrees. Note that heights are ignored and that xy is in km.  Output is [lon, lat, height] in decimal     
    degrees. This is an iterative solution for the inverse of a polyconic projection.
    
    Written by Yi Qujie, Institute of Seismology. 2017
    Modified by Zhao Bin, Institute of Seismology, CEA. April 20, 2018
    
    INPUT: 
        xyh    = np.array([north, east, up])
        origin = np.array([lat, lon])
    OUTPUT:
        llh    = np.array([lat, lon, up])
    
    '''
    
    a = 6378137.0
    e = 0.08209443794970

    # Convert to radians / meters
    xyh    = xyh*1e3
    origin = np.deg2rad(origin)
    llh    = np.zeros(xyh.shape)

    # Iterate to perform inverse projection
    M0     = a*((1-e**2/4-3*e**4/64-5*e**6/256)*origin[0] - (3*e**2/8+3*e**4/32+45*e**6/1024)*np.sin(2*origin[0]) + (15*e**4/256 +45*e**6/1024)*np.sin(4*origin[0]) - (35*e**6/3072)*np.sin(6*origin[0]))
    z      = (xyh[:,1] != -M0)
    for i in range(len(z)):
        if z[i]:
            A    = (M0 + xyh[i,1])/a
            B    = xyh[i,0]**2/a**2 + A**2
            llh[i,1] = A
            delta    = 1e10
            c        = 0
            while abs(delta) > 1e-8:
                C    = np.sqrt((1-e**2*np.sin(llh[i,1])**2))*np.tan(llh[i,1])
                M    = a*((1-e**2/4-3*e**4/64-5*e**6/256)*llh[i,1] - (3*e**2/8+3*e**4/32+45*e**6/1024)*np.sin(2*llh[i,1]) + (15*e**4/256 +45*e**6/1024)*np.sin(4*llh[i,1]) - (35*e**6/3072)*np.sin(6*llh[i,1]))
                Mn   = 1-e**2/4-3*e**4/64-5*e**6/256 -2*(3*e**2/8+3*e**4/32+45*e**6/1024)*np.cos(2*llh[i,1]) +  4*(15*e**4/256 +45*e**6/1024)*np.cos(4*llh[i,1]) -6*(35*e**6/3072)*np.cos(6*llh[i,1])
                Ma   = M/a
                delta= -(A*(C*Ma+1)-Ma-0.5*(Ma**2+B)*C)/(e**2*np.sin(2*llh[i,1])*(Ma**2+B-2*A*Ma)/(4*C)+(A-Ma)*(C*Mn-2/np.sin(2*llh[i,1]))-Mn)
                llh[i,1] = llh[i,1] + delta
                c        = c+1
                if c > 100:
                    print('Convergence failure.')
            llh[i,0] = (np.arcsin(xyh[i,0]*C/a))/np.sin(llh[i,1])+origin[1]
        else:
            # Handle special case of latitude = 0
            llh[i,0] = xyh[i,0]/a + origin[1]
            llh[i,1] = 0

    # Convert back to decimal degrees
    llh[:,0] = np.rad2deg(llh[:,0])
    llh[:,1] = np.rad2deg(llh[:,1])
    llh[:,2] = xyh[:,2]/1e3
    return llh


def geod2xy(o_lon, o_lat, lon, lat):
    rlat  = deg2rad(o_lat)
    scale = cos(rlat)*(6378.139*(1.0-sin(rlat)**2/298.247))
    scale = deg2rad(scale)

    x     = (lon - o_lon)*scale
    y     = (lat - o_lat)*111.32

    return x, y

    
# Written by Zhao Bin, Insitutute of Seismology, CEA. Jan 23th, 2016
def polyconic(Lat, Diff_long, Lat_Orig):
    p1 = Lat_Orig
    p2 = Lat
    il = Diff_long

    arcone = 4.8481368e-6
    esq = 6.7686580e-3
    la = 6378206.4
    a0 = 6367399.7
    a2 = 32433.888
    a4 = 34.4187
    a6 = .0454
    a8 = 6.0e-5

    ip = p2 - p1
    sinp2 = sin(p2 * arcone)
    cosp2 = cos(p2 * arcone)
    theta = il * sinp2
    a = sqrt(1.0 - (esq * (2. * sinp2))) / (la * arcone)
    cot = cosp2 / sinp2
    x = (cot * sin(theta * arcone)) / (a * arcone)
    ipr = ip *arcone
    pr = ((p2 + p1) /2.) * arcone
    y = ((((a0*ipr) - ((a2*cos(2.*pr))*sin(ipr))) + ((a4*cos(4.*pr))*sin(2.*ipr))) - ((a6*cos(6.*pr))*sin(3.*ipr))) + ((a8*cos(8.*pr))*sin(4.*ipr))
    xy = array([x, y])
 
    return xy
    
    
    
def xy2llh(xy, origin):
    '''
    Given the model calculated relative to the lat/lon origin, return lat/lon/change
    Mod by Zhao Bin, June 1, 2016 to fix debug

    Input parameters:
    xy     -  numpy array with two column, north and east, unit is kilometer
    origin -  numpy array with two column, lat and lon in degree

    Output:
    llh    -  numpy array, lat and lon in degree
    '''


    #N = radius of earth, e = eccentricity
    N = 6378206.4
    e = 1./298.257

    W = sqrt(1 - e**2*sin(deg2rad(origin[0])))
    Rp = (N/W)*cos(deg2rad(origin[0]))
    Rm = N*(1-e**2)/W**3

    # We first convert array to matrix, so we can get the shape of xy
    xy = np.mat(xy)
    [m, n] = xy.shape

    # Now we convert the matrix to array
    xy = np.array(xy)
    llh = np.zeros((m, 2))

    llh[:,0] = origin[0] + rad2deg(1.0e3*xy[:,0]/Rm)
    llh[:,1] = origin[1] + rad2deg(1.0e3*xy[:,1]/Rp)
    
    return llh 


def vpvs2nu(vp, vs):
    '''
    compute Poisson's ratio using vp and vs
    
    Input:
       vp   - velocity of P wave
       vs   - velocity of S wave
    Output:
       nu   - Poisson's ratio
    '''

    nu = (vp**2-2*vs**2)/(vp**2-vs**2)*0.5
    return nu

def vpvsrho2K(vp, vs, rho):
    '''
    compute bulk modulus using vp and vs and rho

    Input:
       vp   - velocity of P wave
       vs   - velocity of S wave
       rho  - density of rock
    Output:
       K    - bulk modulus
    '''

    K = rho*(vp**2 - 4*vs**2/3)
    return K


def vsrho2mu(vs, rho):
    '''
    compute shear modulus mu using Vs and rho

    Input:
        vs   - velocity of S save
        rho  - density of crust
    Output:
        mu   - shear modulus

    '''
    mu = rho*vs**2
    return mu

def lameshear2nu(lame, shear):
    '''
    compute Poisson's ratio using Lame's coefficient and shear modulus
    Written by Zhao Bin, Instutite of Seismology, CEA. Mar 3, 2017

    Input:
        lame  - Lame coefficient
        shear - shear modulus
    Output:
        nu    - Possion's ratio
    '''
    vp2_vs2 = (lame+2*shear)/shear
    nu = (vp2_vs2-2)/(vp2_vs2-1)/2.0
    return nu

def rotate2d(xy, angle):
    angle   = deg2rad(angle)
    localxy = np.zeros(2)
    localxy[0] =-xy[0]*sin(angle) + xy[1]*cos(angle)
    localxy[1] = xy[0]*cos(angle) + xy[1]*sin(angle)
    return localxy


def transect(geodata, slon, slat, azimuth, width):
    '''
    make a profile section with a width from start point and end point
    Written by Zhao Bin, Institute of Seismology @UC Berkeley May 3, 2016
    MOD     by Zhao Bin, Institute of Seismology, Aug 2, 2017
    MOD     by Zhao Bin, Institute of Seismology, Oct 12, 2017

    input:
    geodata    -- numpy array, containing at least two column lon and lat
               -- also contain other information.
               -- it is convient for InSAR, earthquake catlog, GPS data 
    slon       -- original longitude in degree
    slat       -- original latitude in degree
    azimuth    -- profile azimuth from north, clockwise is positive
    width      -- in kilometer 

    output:
    '''
    # origial point
    origin = array([slat, slon])
    angle  = deg2rad(90-azimuth)
    localxy = zeros((len(geodata),2))
    for i in range(0, len(geodata)):
        llh = array([geodata[i,1], geodata[i,0]]) 
        xy  = llh2localxy(llh, origin)

        localxy[i,0] =-xy[0]*sin(angle) + xy[1]*cos(angle)
        localxy[i,1] = xy[0]*cos(angle) + xy[1]*sin(angle)
    combdata = np.hstack((localxy, geodata))
    indx     = np.where(np.abs(combdata[:,0])<width)[0]
    combdata = combdata[indx]
    m,n      = combdata.shape
    np.savetxt('transect.dat', combdata, fmt=n*"%10.4f")

def distazimuth2point(lon1, lat1, distance, azimuth, R):
    '''
    compute end point (longitude and latitude) given a starting point, distance, and azimuth
    By Zhao Bin, Institute of Seismology, CEA. @ UC Berkeley, June 9, 2016

    Input:
    lon1      -- starting point longitude in degree
    lat1      -- starting point latitude in degree
    distance  -- distance in kilometer
    azimuth   -- clockwise from north in degree
    R         -- radius of the Earth/shpere in kilometer
  
    Output:
    lon2      -- ending point longitude in degree
    lat2      -- ending point latitude in degree

    Note: http://www.codeguru.com/cpp/cpp/algorithms/article.php/c5115/Geographic-Distance-and-Azimuth-Calculations.htm
    '''
    b   = float(distance)/float(R)
    a   = cos(b)*cos(deg2rad(90-lat1)) + sin(deg2rad(90-lat1))*sin(b)*cos(deg2rad(azimuth))
    a   = np.arccos(a)
    B   = sin(b)*sin(deg2rad(azimuth))/sin(a)
    B   = np.arcsin(B)

    lat2 = 90-np.rad2deg(a)
    lon2 = np.rad2deg(B)+lon1

    return lon2, lat2


def makegrids(slon, elon, slat, elat, nx, ny):
   '''
   make one dimension profile or two dimension grid
   
   Input:
   slon   -- float, start longitude in degree
   elon   -- float, end longitude in degree
   slat   -- float, start latitude in degree
   elat   -- float, end latitude in degree
   nx     -- int, number of points along longitude
   ny     -- int, number of points along latitude

   Output:
   grid.txt  -- lon and lat
   '''

   lon = linspace(slon, elon, nx)
   lat = linspace(slat, elat, ny)

   grid = zeros((nx*ny,2))
   k    = 0
   for i in range(0, nx):
       for j in range(0, ny):
          pt   = array([lon[j], lat[i]])
          grid[k,:] = pt
          k=k+1

   filename = 'grid.txt'
   savetxt(filename, grid, fmt="%10.2f\t%10.2f")


def make_grid_3D(slon, elon, slat, elat, nx, ny):
   '''
   make one dimension profile or two dimension grid
   
   Input:
   slon   -- float, start longitude in degree
   elon   -- float, end longitude in degree
   slat   -- float, start latitude in degree
   elat   -- float, end latitude in degree
   nx     -- int, number of points along longitude
   ny     -- int, number of points along latitude

   Output:
   grid.txt  -- lon and lat
   '''

   lon = linspace(slon, elon, nx)
   lat = linspace(slat, elat, ny)

   grid = zeros((nx*ny,3))
   k    = 0
   for i in range(0, nx):
       for j in range(0, ny):
          pt   = array([lat[i], lon[j], 0.0])
          grid[k,:] = pt
          k=k+1

   filename = 'grid.txt'
   savetxt(filename, grid, fmt="%10.2f\t%10.2f\t%10.2f")




def Moment(Len, Wid, slip, shearmodulus):
    '''
    Estimate moment and moment magnitude
    Mod by Zhao Bin, Apr. 24, 2019. Now 6.067 is used.

    Input:
        Len        -- Length of fault in km
        Wid        -- Width of fault in km
        shearmodulus:    
        
    Output:
        Mo:
        Mw:
    
    '''

    import numpy as np

    Mo_total = 1000*Len * 1000*Wid * slip * shearmodulus
    Mw_total = 2.0/3.0*np.log10(Mo_total) - 6.067

#   print Mo_total, Mw_total
    return Mo_total, Mw_total

def vrot(east, north, azimuth):
    '''
    rotate east-north coordinate with azimuth
    By Zhao Bin, Institute of Seismology, CEA. @ UC Berkeley. August 16, 2016

    Input:
    east    : x
    north   : y
    azimuth : counterclockwise from east in degree

    Output:
    '''
    angle = deg2rad(azimuth)

    a     =  east*cos(angle) + north*sin(angle)
    b     = -east*sin(angle) + north*cos(angle)
    return a, b


def mdiag(a, b):
    '''
    create a new diag matrix using a and b.
    '''
    n = len(a)
    m = len(b)
    c = np.zeros((n,m))
    d = c.T

    up = np.hstack((a,c))
    dn = np.hstack((d,b))
    diag = np.vstack((up, dn))
    return diag

def ml2mw(ml):
    '''
    Written by Zhao Bin, Instutute of Seismology, CEA. Oct. 12 2017
    Convert local magnitude to moment magnitude
    The relationship is from Meng et al., 2015, EPSL
    '''
    if ml <= 4.41:
        mw = ml*2.0/3.0 + 1.56
    elif ml > 4.41 and ml < 5.71:
        mw = ml+0.09

    return mw

def mb2ml(mb):
    '''
    Written by Zhao Bin, Institute of Seismology, CEA. Oct. 20 2017
    Convert Mb to Ml
    The relationship is from https://wenku.baidu.com/view/b4e1dd37f12d2af90242e65c.html
    '''
    ml = (1.17*mb+0.67)/1.13
    return ml


def mw2mo(mw):
    '''
    Written by Zhao Bin, Instutute of Seismology, CEA. Oct. 12 2017
    Convert local magnitude to moment magnitude
    unit for moment magnitude is Nm
    '''
    Mo = 10**((mw+6.0)*1.5)
    return Mo

def mw2ra(mw):
    '''
    Written by Zhao Bin, Instutute of Seismology, CEA. Oct. 12 2017
    Convert Mw to rupture area following the relationship by Wells et al., 1994
    unit for area is km^2
    '''
    ra = 10**(0.91*mw-3.49)
    return ra

def mw2as(mw):
    '''
    Written by Zhao Bin, Instutute of Seismology, CEA. Oct. 12 2017
    Convert Mw to average slip
    unit for average slip is mm
    '''
    ra = mw2ra(mw)
    mo = mw2mo(mw)
    meanslip = mo/ra/1.0E3/3E10
    return meanslip

if __name__ == '__main__':
    pass
