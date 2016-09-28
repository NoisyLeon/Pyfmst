import field2d_earth
import GeoPolygon
import pyasdf
# basins=GeoPolygon.GeoPolygonLst()
# basins.ReadGeoPolygonLst('basin1')

# minlat=23.
# maxlat=52.
# minlon=85.
# maxlon=133.

minlat=20.
maxlat=52.
minlon=80.
maxlon=134.


field=field2d_earth.Field2d(minlon=minlon, maxlon=maxlon, dlon=0.5, minlat=minlat, maxlat=maxlat, dlat=0.5, period=10.)

field.read(fname='./Tph_10sec_0.5.lst')
workingdir='./field_working_0.5'
field.interp_surface(workingdir=workingdir, outfname='Tph_10sec')
field.check_curvature(workingdir=workingdir)
field.gradient_qc(workingdir=workingdir, evlo=129.0, evla=41.306, nearneighbor=False)
# field.reset_reason()
# field.plot_field(contour=True, geopolygons=basins)
# field.np2ma()
# inraydset=pyasdf.ASDFDataSet('../ray_0.5.h5')
field.plot_diffa_special(factor=2, rayASDF='../rays_0.5.h5', rayfactor=1)
# evlo=129.0
# evla=41.306
# field.replace_appv_lplc(evlo=evlo, evla=evla, invfname='inPhV.lst', cdist=300.)
# field.plot_appV(geopolygons=basins)
# field.plot_lplc(vmin=-0.005, vmax=0.005, geopolygons=basins)
# field.write_dbase(outdir='./fmst_dbase_0.1')

