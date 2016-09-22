import field2d_earth
import GeoPolygon
import raypath
import time

# 
# basins=GeoPolygon.GeoPolygonLst()
# basins.ReadGeoPolygonLst('basin1')
# 
minlat=20.
maxlat=52.
minlon=80.
maxlon=134.
# 
field=field2d_earth.Field2d(minlon=minlon, maxlon=maxlon, dlon=0.1, minlat=minlat, maxlat=maxlat, dlat=0.1, period=10.)
# field.read_dbase(datadir='./fmst_dbase_0.2')
# dset=raypath.rayASDF('../rays_0.2.h5')
# dset.read_raydat('/projects/life9360/code/Pyfmst/fmm_working_0.2/gmtplot/rays.dat')
# dset.get_pre_amp(field2d=field, outfname='pre_amp_0.2.lst')

field.read_dbase(datadir='./fmst_dbase_0.2')
dset=raypath.rayASDF('../rays_0.2.h5')
# dset.read_raydat('/projects/life9360/code/Pyfmst/fmm_working/gmtplot/rays.dat')
dset.get_pre_amp(field2d=field, outfname='pre_amp_0.2_0.2.lst')