import field2d_earth
import GeoPolygon
basins=GeoPolygon.GeoPolygonLst()
basins.ReadGeoPolygonLst('basin1')

minlat=23.
maxlat=52.
minlon=85.
maxlon=133.

# minlat=20.
# maxlat=52.
# minlon=80.
# maxlon=134.


field=field2d_earth.Field2d(minlon=minlon, maxlon=maxlon, dlon=0.2, minlat=minlat, maxlat=maxlat, dlat=0.2, period=10.)


# field.read_dbase(datadir='./output')
# # field.read(fname='./stf_10_20sec/Tph_10.0.txt')
# field.read(fname='./stf_10sec_all/Tph_10.0.txt')
# field.read(fname='./Tph_10sec_0.5.lst')
field.read(fname='./Tph_10sec_0.5.lst')
# field.add_noise(sigma=5.)
workingdir='./field_working'
field.interp_surface(workingdir=workingdir, outfname='Tph_10sec')
field.check_curvature(workingdir=workingdir)
field.gradient_qc(workingdir=workingdir, evlo=129.0, evla=41.306, nearneighbor=False)
# field.reset_reason()
field.Laplacian_Green()
# field.plot_field(contour=True, geopolygons=basins)
field.np2ma()
# field.plot_diffa()
# field.plot_appV(geopolygons=basins, vmin=2.9, vmax=3.4)
field.plot_lplc(vmin=0.005, vmax=-0.005)
# field.write_dbase(outdir='./fmst_dbase_0.2')
# field.get_distArr(evlo=129.0,evla=41.306)
# field.write_dbase(outdir='./output_ses3d_all6')



# field.read_dbase(datadir='./output_ses3d')
# field.np2ma()
# fieldFMM=field.copy()
# fieldFMM.read_dbase(datadir='./output_FMM')
# fieldFMM.np2ma()
# # fieldFMM.plot_diffa()
# field.compare(fieldFMM)
# # field.histogram()
# field.mean_std()
# field.plot_compare()

# field.Laplacian('convolve')
# field.plot_lplc()
