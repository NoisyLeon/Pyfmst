import field2d_earth

minlat=25.
maxlat=50.
minlon=235.
maxlon=295.


# field.read(fname='../Pyfmst/US_phV/8sec_lf')
# field.interp_surface(workingdir='./field_working_US', outfname='8sec_v')


field=field2d_earth.Field2d(minlon=minlon, maxlon=maxlon, dlon=0.5, minlat=minlat, maxlat=maxlat, dlat=0.5, period=8.)
# field.LoadFile(fname='./stf_100_10sec_all/Tph_10.0.txt')
# field.LoadFile(fname='./stf_10sec_all/Tph_10.0.txt')
# field.interp_surface(workingdir='./field_working', outfname='Tph_10sec')

# field.read_dbase(datadir='./output')
# field.read(fname='./stf_10sec_all/Tph_10.0.txt')
field.read(fname='./Tph_8sec_US_0.5.lst')
# field.add_noise(sigma=5.)
field.interp_surface(workingdir='./field_working_US', outfname='Tph_8sec')
field.check_curvature(workingdir='./field_working')
# 40.195000, 247.1867
field.gradient_qc(workingdir='./field_working', evlo=247.1867, evla=40.195000, nearneighbor=False)
field.np2ma()
field.plot_diffa(projection='merc')
# field.plot_appV()
# field.plot_lplc()
# field.write_dbase(outdir='./output', gmt=True)

# field.Laplacian('convolve')
# field.plot_lplc()
# field.plot_field(contour=True)