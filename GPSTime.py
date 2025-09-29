# -*- coding: utf-8 -*-
"""
GPS time conversions 
Created on Wed Apr 18 17:01:10 2018

@author: zhao
"""

# load libs
import numpy as np


def jd_to_decyrs(jd):
    '''
    Convert Julian day to decimal year
    
    Written by Zhao Bin, Institute of Seismology, CEA. April 19. 2018
    
    INPUT:
        jd     = Julian day
    OUTPUT:
        decyrs = decmail year
    '''
    
    if jd < 2000000.0:
        jmd = jd + 2400000.5
    else:
        jmd = jd
    
    date, seconds, doy = jd_to_ymdhms(jd)
    
    date[1] = 1
    date[2] = 1
    date[3] = 0
    date[4] = 0
    secs    = 0
    
    jd       = ymdhms_to_jd(date, secs)
    date[0]  = date[0] + 1
    jde      = ymdhms_to_jd(date, secs)
    num_days = jde - jd
    
    if num_days <= 365.0 and num_days <= 366.0:
        num_days = 365.0
    date[0]  = date[0] - 1
    decyrs   = date[0] + (jmd-jd)/num_days
    
    return decyrs
    


def jd_to_ymdhms(jd):
    '''
    Convert Julian day to year-month-day-hour-minute-seconds and doy
    
    Written by Zhao Bin, Institute of Seismology, CEA. April 19. 2018
    
    INPUT:
        jd      = Julian day
    OUTPUT:
        date[0] = year
        date[1] = month
        date[2] = day
        date[3] = hour
        date[4] = minute
        seconds = seconds
        day_of_year  = day of year
    '''
    
    date          = np.zeros(5, dtype=int)
    days_to_month = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334, 365]
    if jd >= 2000000.0:
        mjd       = jd - 2400000.5
    else:
        mjd       = jd
    
    mjd_day       = int(mjd)
    fraction      = np.mod(mjd, 1.0)
    
    if mjd < 0 and fraction != 0.0:
        mjd_day   = mjd_day - 1
        fraction  = fraction + 1.0
        
    days_from_1600  = int(mjd_day - (-94554.0))
    years_from_1600 = int(days_from_1600/365.0)
    
    day_of_year   = 0
    while day_of_year <= 0.0:
        century   = years_from_1600//100
        day_of_year = days_from_1600 - years_from_1600*365 - (years_from_1600-1)//4 + (years_from_1600+99)//100 - (years_from_1600+399)//400 - 1
        if years_from_1600 == 0:
            day_of_year     = day_of_year + 1
        if day_of_year <= 0.0:
            years_from_1600 = years_from_1600 - 1
    
    year          = np.mod(years_from_1600, 100)
    
    leap_year = False
    if year == 0:
        if np.mod(century, 4) == 0:
            leap_year = True
    else:
        if np.mod(year, 4) == 0:
            leap_year = True
            
    if day_of_year < 60:
        if day_of_year <= 31:
            month     = 1
            day       = day_of_year
        else:
            month     = 2
            day       = day_of_year - 31
    else:
        if leap_year and day_of_year == 60:
            month     = 2
            day       = 29
        else:
            if leap_year:
                day_of_year = day_of_year - 1
            month     = 2
            while day_of_year > days_to_month[month-1]:
                month = month + 1
            month     = month-1
            day       = day_of_year - days_to_month[month-1]
            
            
    date[0] = years_from_1600 + 1600
    date[1] = month
    date[2] = day
    date[3] = fraction*24.0
    date[4] = fraction*1440.0 - date[3]*60.0
    
    seconds = 86400.0*fraction - date[3]*3600.0 - date[4]*60.0
    
    if seconds >= 59.0:
        dsec = 1e-6
        fracp = fraction + dsec/86400.0
        date[3] = fracp*24.0
        date[4] = fracp*1440.0 - date[3]*60.0
        
        seconds = 86400.0*fracp - date[3]*3600.0 - date[4]*60.0 -dsec
            
    return date, seconds, day_of_year

def ymdhms_to_jd(date, seconds):
    '''
    Written by He Kefeng, Institute of Seismology, CEA. 2017.
    Fix a bug by Zhao Bin, Aug. 4, 2020
    
    INPUT:
        date[0] = year
        date[1] = month
        date[2] = day
        date[3] = hour
        date[4] = minute
        seconds = seconds
    OUTPUT:
        epoch   = Julian day
    '''
    year  = int(date[0])
    month = int(date[1])
    day   = int(date[2])
    hour  = int(date[3])
    minute= int(date[4])
    
    days_to_month = [0, 31, 59, 90, 120, 151, 181, 212, 243, 273, 304, 334]
    if year < 50:
        year = year + 2000
    elif year < 200:
        year = year + 1900

    years_from_1600 = year-1600
    leap_days = (years_from_1600 -1)//4 - (years_from_1600+99)//100 + (years_from_1600+399)//400+1
    if (years_from_1600 == 0):
        leap_days = leap_days-1

    leap_year = False

    if (np.mod(years_from_1600,4) == 0 and (np.mod(years_from_1600,100)!=0 or np.mod(years_from_1600,400)==0)):
        leap_year = True

    days_from_1600 = years_from_1600*365 + leap_days + days_to_month[month-1] + day
    if month>2 and leap_year:
        days_from_1600 = days_from_1600+1

    fraction = seconds/86400.0 + minute/1440.0 + hour/24.0

    mjd = -94554.0 + days_from_1600 + fraction

    epoch = mjd + 2400000.5

    return epoch



def ymd_to_decyrs(year, month, day):
    '''
    Convert YYYY-MM-DD to decimal year
    Written by He Kefeng, Institute of Seismology, CEA. 2017.
    Modified by Zhao Bin, July 17, 2018. Fix bug in the code
    
    INPUT:
        year    = year
        month   = month
        day     = day
        
    OUTPUT:
        decyrs  = decimal year
    
    '''

    jdi = ymdhms_to_jd([year, month, day, 0, 0], 0)

    if jdi < 2000000.0:
        jdm = jdi + 2400000.5
    else:
        jdm = jdi


    jd = ymdhms_to_jd([year, 1, 1, 0, 0], 0)
    year = year+1
    jde= ymdhms_to_jd([year, 1, 1, 0, 0], 0)

    num_days = jde - jd

    if (num_days <= 365. and num_days <=366.0):
        num_days = 365.0

    year = year-1
    decyrs = year + (jdm-jd)/num_days
    return decyrs

def decyrs_to_mjd(decyrs):
    '''
    Convert from decimal year to Julian day

    Input:
        decyrs    = decimal year
    Output:
        mjd       = Julian day
    '''

    jd      = decyrs_to_jd(decyrs)
    mjd     = jd-2400000.0
    return mjd

def decyrs_to_jd(decyrs):
    '''
    Convert from decimal year to Julian day

    Input:
        decyrs    = decimal year
    Output:
        jd        = Julian day
    '''

    date    = np.zeros(5)
    date[0] = int(decyrs)
    date[1] = 1
    date[2] = 1
    secs    = 0

    jds     = ymdhms_to_jd( date, secs)
    date[0] += 1
    jde     = ymdhms_to_jd( date, secs)

    jd      = jds + (decyrs-int(decyrs))*(jde-jds)

    return jd
