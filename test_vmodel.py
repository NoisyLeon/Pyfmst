import pyfmst

minlat=19
maxlat=52.5
minlon=74.5
maxlon=143.5

vmodel=pyfmst.vmodel(minlat=minlat, maxlat=maxlat, minlon=minlon, maxlon=maxlon, dlat=0.5, dlon=0.5, v0=3.0)
vmodel.read_cv('10.phase_interp.map.HD')
vmodel.smooth(1)
# vmodel.read_cv('/projects/life9360/China_2015_disp_v1.0/10.phase.map')
vmodel.plot()
vmodel.write(outfname='10.phase_smooth.map')
