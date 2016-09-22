import field2d_earth
import GeoPolygon
basins=GeoPolygon.GeoPolygonLst()
basins.ReadGeoPolygonLst('basin1')



minlat=20.
maxlat=52.
minlon=80.
maxlon=134.


field=field2d_earth.Field2d(minlon=minlon, maxlon=maxlon, dlon=0.2, minlat=minlat, maxlat=maxlat, dlat=0.2, period=10., fieldtype='Amp')


# field.read_dbase(datadir='./output')
# # field.read(fname='./stf_10_20sec/Tph_10.0.txt')
# field.read(fname='./stf_10sec_all/Tph_10.0.txt')
field.read(fname='./pre_amp_0.5_0.5.lst')
# field.add_noise(sigma=5.)
workingdir='./field_working'
field.interp_surface(workingdir=workingdir, outfname='Amp_10sec')
# field.check_curvature(workingdir=workingdir)
# field.gradient_qc(workingdir=workingdir, evlo=129.0, evla=41.306, nearneighbor=False)
# field.reset_reason()
# field.plot_field(contour=False, geopolygons=basins)
# field.np2ma()
# field.plot_diffa()
# field.plot_appV(geopolygons=basins)
# field.plot_lplc(vmin=0.005, vmax=-0.005)
# field.write_dbase(outdir='./fmst_dbase')
# field.get_distArr(evlo=129.0,evla=41.306)
# field.write_dbase(outdir='./output_ses3d_all6')
evlo=129.0
evla=41.306
field.get_dist_data(evlo=evlo, evla=evla, outfname='./compare_amp/pre_amp_1675_1725.lst', mindist=1675, maxdist=1725, minazi=-1., maxazi=361., geopolygons=basins)
field.get_dist_data(evlo=evlo, evla=evla, outfname='./compare_amp/pre_amp_1075_1125.lst', mindist=1075, maxdist=1125, minazi=-1., maxazi=361., geopolygons=basins)
import matplotlib.pyplot as plt
plt.show()