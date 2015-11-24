#!/usr/bin/env python

from datetime import datetime, timedelta
import os
import urllib.request

indate = datetime(2015,11,20,21)

left = -125.5
right = -65.5
top = 50
bottom = 24.5


cores = ('arw','nmb')
perts = ('ctl','n1','n2','n3','n4','n5','n6','p1','p2','p3','p4','p5','p6')

base = 'http://nomads.ncep.noaa.gov/cgi-bin/filter_srefbc.pl?file=sref_{core:s}.t{init:%H}z.pgrb212.{pert:s}.grib2&lev_2_m_above_ground=on&lev_10_m_above_ground=on&lev_surface=on&lev_850_mb=on&lev_mean_sea_level=on&var_HGT=on&var_PRMSL=on&var_TMP=on&var_UGRD=on&var_VGRD=on&var_APCP=on&subregion=&leftlon={left:f}&rightlon={right:f}&toplat={top:f}&bottomlat={bottom:f}&dir=%2Fsref.{init:%Y%m%d}%2F{init:%H}%2Fpgrb_biasc'


os.system('rm -f *.nc')
os.system('rm -f *.grib2')

for core in cores:
    for pert in perts:
        curdat = {'init': indate,
                  'top': top,
                  'bottom': bottom,
                  'left' : left,
                  'right' : right,
                  'core' : core,
                  'pert' : pert}
        fullpath = base.format(**curdat)
        filename = 'sref_{core:s}_{pert:s}.grb2'.format(**curdat)
        ncfile = filename.replace('grb2','nc')
        if os.path.exists(ncfile):
            continue
        # Retrieve the file
        urllib.request.urlretrieve(fullpath, filename)
        # Convert to netcdf
        os.system('wgrib2 {:s} -netcdf {:s}'.format(filename, ncfile))

