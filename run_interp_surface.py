import field2d_earth

minlat=19
maxlat=52.5
minlon=74.5
maxlon=143.5




field=field2d_earth.Field2d(minlon=minlon, maxlon=maxlon, dlon=0.5, minlat=minlat, maxlat=maxlat, dlat=0.5, period=10.)
field.read(fname=('./China_2015_disp_v1.0/10.phase.map'))
field.interp_surface(workingdir='.', outfname='10.phase_interp.map')

# field.plot_appV()
# field.plot_lplc()
# field.write_dbase(outdir='./output', gmt=True)

# field.Laplacian('convolve')
# field.plot_lplc()
# field.plot_field(contour=True)