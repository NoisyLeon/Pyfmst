import pyfmst

outdir='./test'
SLst=pyfmst.StaLst()
SLst.HomoStaLst(minlat=-44, Nlat=38, minlon=112, Nlon=48, dlat=1., dlon=1.)
# SLst.HomoStaLst(minlat = 22, Nlat = 150, minlon=85, Nlon = 240, dlat=0.2, dlon=0.2)
SLst.write(outdir=outdir)

ELst=pyfmst.StaLst()
# ELst.append(pyfmst.StaInfo('E001',41.306, 129.))
ELst.append(pyfmst.StaInfo('E001',-33, 129.))
ELst.write(outdir=outdir, issource=True)

SLst.write_otime(outdir=outdir, ELst=ELst)

minlat=22.
maxlat=52.
minlon=85.
maxlon=133.


# vmodel=pyfmst.vmodel(minlat=minlat, maxlat=maxlat, minlon=minlon, maxlon=maxlon, dlat=0.5, dlon=0.5)
vmodel=pyfmst.vmodel(minlat=-45, maxlat=-5, minlon=110, maxlon=160, dlat=1., dlon=1.)
vmodel.write('./test')
vmodel.read_cv('./testmodel.txt')

# Convert rtravel.out to txt field2d file
SLst.travel2field2d(datadir=outdir, outfname='./T.lst')