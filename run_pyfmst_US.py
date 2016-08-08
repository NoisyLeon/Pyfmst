import pyfmst


minlat=20.
maxlat=53.
minlon=80.
# maxlon=135.
maxlon=135.
outdir='./fmm_working'
SLst=pyfmst.StaLst()
# SLst.HomoStaLst(minlat=minlat+0.5, Nlat=32, minlon=minlon+0.5, Nlon=54, dlat=1., dlon=1.)
SLst.HomoStaLst(minlat=minlat+0.5, Nlat=64, minlon=minlon+0.5, Nlon=108, dlat=.5, dlon=.5)
# # SLst.HomoStaLst(minlat = 22, Nlat = 150, minlon=85, Nlon = 240, dlat=0.2, dlon=0.2)
SLst.write(outdir=outdir)
# 
ELst=pyfmst.StaLst()
# # ELst.append(pyfmst.StaInfo('E001',41.306, 129.))
ELst.append(pyfmst.StaInfo('E001',41.306, 129.))
ELst.write(outdir=outdir, issource=True)
# 
SLst.write_otime(outdir=outdir, ELst=ELst)
# 
# # minlat=22.
# # maxlat=52.
# # minlon=85.
# # maxlon=133.
# 
# 
vmodel=pyfmst.vmodel(minlat=minlat, maxlat=maxlat, minlon=minlon, maxlon=maxlon, dlat=0.5, dlon=0.5, v0=3.0)
# # # vmodel=pyfmst.vmodel(minlat=-45, maxlat=-5, minlon=110, maxlon=160, dlat=1., dlon=1.)
vmodel.read_cv('/projects/life9360/China_2015_disp_v1.0/10.phase.map')
vmodel.smooth(100)
vmodel.read_cv('/projects/life9360/China_2015_disp_v1.0/10.phase.map')
vmodel.plot()
vmodel.write(outdir=outdir)

# Convert rtravel.out to txt field2d file
# SLst.travel2field2d(datadir=outdir, outfname='./Tph_10sec_0.5.lst')