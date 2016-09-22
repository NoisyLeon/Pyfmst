import field2d_earth
import matplotlib.pyplot as plt
import GeoPolygon
basins=GeoPolygon.GeoPolygonLst()
basins.ReadGeoPolygonLst('basin1')

minlat=23.
maxlat=51.
minlon=86.
maxlon=132.

field=field2d_earth.Field2d(minlon=minlon, maxlon=maxlon, dlon=0.2, minlat=minlat, maxlat=maxlat, dlat=0.2, period=10., fieldtype='Amp')
# Tfield=field2d_earth.Field2d(minlon=minlon, maxlon=maxlon, dlon=0.2, minlat=minlat, maxlat=maxlat, dlat=0.2, period=10., fieldtype='Tph')

# Tfield.read_dbase(datadir='./output_ses3d_all6')
# field.read(fname='./stf_10_20sec/Amp_10.0.txt')
field.read(fname='../SES3DPy/stf_10sec_all/Amp_10.0.txt')
field.ZarrIn=field.ZarrIn/10.
# # field.read(fname='../Pyfmst/Tph_10sec_0.5.lst')
# field.add_noise(sigma=5.)
workingdir='./field_working'
field.interp_surface(workingdir=workingdir, outfname='Amp_10sec')
# field.check_curvature(workingdir=workingdir, threshold=20.)
# field.gradient_qc(workingdir=workingdir, evlo=129.0, evla=41.306, nearneighbor=False)
# field.reset_reason()
# field.plot_field_sta(contour=False, geopolygons=basins)
# # field.np2ma()
# # field.plot_diffa()
# # field.plot_appV()
# field.Laplacian()
# # field.np2ma()
# field.plot_field(contour=False, geopolygons=basins)
# field.plot_lplcC(infield=Tfield)

evlo=129.0
evla=41.306
field.get_dist_data(evlo=evlo, evla=evla, outfname='./compare_amp/obs_amp_1075_1125.lst', mindist=1075, maxdist=1125, minazi=-1., maxazi=361., geopolygons=basins)