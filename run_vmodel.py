import pyfmst
import GeoPolygon
import field2d_earth
basins=GeoPolygon.GeoPolygonLst()
basins.ReadGeoPolygonLst('basin1')
minlat=20.
maxlat=53.
minlon=80.
# maxlon=135.
maxlon=135. 
vmodel=pyfmst.vmodel(minlat=minlat, maxlat=maxlat, minlon=minlon, maxlon=maxlon, dlat=0.5, dlon=0.5, v0=3.0)
vmodel.read_cv('./China_2015_disp_v1.0/10.phase.map')
vmodel.smooth(100)
vmodel.read_cv('./China_2015_disp_v1.0/10.phase.map')
vmodel.write4field2d('10.phase_extended_map')

minlat=20.2
maxlat=51.8
minlon=80.2
maxlon=133.8


field=field2d_earth.Field2d(minlon=minlon, maxlon=maxlon, dlon=0.1, minlat=minlat, maxlat=maxlat, dlat=0.1, period=10., fieldtype='Tph')

field.read(fname='10.phase_extended_map')
workingdir='./field_working'
field.interp_surface(workingdir=workingdir, outfname='Tph_10sec')
field.write('inPhV.lst', 'txt')
# field.plot_field(contour=False, geopolygons=basins)