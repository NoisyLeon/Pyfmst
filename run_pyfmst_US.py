import pyfmst


minlat=25.
maxlat=50.
minlon=235.
maxlon=295.
outdir='./fmm_working_US'
SLst=pyfmst.StaLst()
# SLst.HomoStaLst(minlat=minlat+0.5, Nlat=32, minlon=minlon+0.5, Nlon=54, dlat=1., dlon=1.)
SLst.HomoStaLst(minlat=minlat+0.5, Nlat=50, minlon=minlon+0.5, Nlon=120, dlat=.5, dlon=.5)
# # SLst.HomoStaLst(minlat = 22, Nlat = 150, minlon=85, Nlon = 240, dlat=0.2, dlon=0.2)
SLst.write(outdir=outdir)
# 
ELst=pyfmst.StaLst()
# # ELst.append(pyfmst.StaInfo('E001',41.306, 129.))
ELst.append(pyfmst.StaInfo('E001',47.1523, 249.313))
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
vmodel.read_cv('./US_phV/8sec_v.HD')
# vmodel.smooth(100)
# vmodel.read_cv('./US_phV/8sec_v2_iso')
vmodel.plot()
vmodel.write(outdir=outdir)

# Convert rtravel.out to txt field2d file
# SLst.travel2field2d(datadir=outdir, outfname='./Tph_8sec_US_0.5.lst')