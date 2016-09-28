#!/usr/bin/env python
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.mlab import griddata
import numpy.ma as ma
import scipy.ndimage.filters 
# from skimage.filters import roberts, sobel, scharr, prewitt
from scipy.ndimage import convolve
import matplotlib
import multiprocessing
from functools import partial
from lasif import colors
import os
from subprocess import call
import obspy.geodetics
from mpl_toolkits.basemap import Basemap, shiftgrid, cm
from pyproj import Geod
import random
import copy
import pycpt
import pyasdf

lon_diff_weight_2 = np.array([[1., 0., -1.]])/2.
lat_diff_weight_2 = lon_diff_weight_2.T
lon_diff_weight_4 = np.array([[-1., 8., 0., -8., 1.]])/12.
lat_diff_weight_4 = lon_diff_weight_4.T
lon_diff_weight_6 = np.array([[1./60., 	-3./20.,  3./4.,  0., -3./4., 3./20.,  -1./60.]])
lat_diff_weight_6 = lon_diff_weight_6.T

lon_diff2_weight_2 = np.array([[1., -2., 1.]])
lat_diff2_weight_2 = lon_diff2_weight_2.T
lon_diff2_weight_4 = np.array([[-1., 16., -30., 16., -1.]])/12.
lat_diff2_weight_4 = lon_diff2_weight_4.T
lon_diff2_weight_6 = np.array([[1./90., 	-3./20.,  3./2.,  -49./18., 3./2., -3./20.,  1./90.]])
lat_diff2_weight_6 = lon_diff2_weight_6.T

geodist = Geod(ellps='WGS84')

def discrete_cmap(N, base_cmap=None):
    """Create an N-bin discrete colormap from the specified input map"""
    # Note that if base_cmap is a string or None, you can simply do
    #    return plt.cm.get_cmap(base_cmap, N)
    # The following works for string, None, or a colormap instance:
    base = plt.cm.get_cmap(base_cmap)
    color_list = base(np.linspace(0, 1, N))
    cmap_name = base.name + str(N)
    return base.from_list(cmap_name, color_list, N)

class Field2d(object):
    """
    An object to analyze 2D spherical field data on Earth
    ===========================================================================
    Parameters:
    dlon, dlat      - grid interval
    Nlon, Nlat      - grid number in longitude, latitude 
    lonArr, latArr  - arrays for grid location
    fieldtype       - field type (Tph, Tgr, Amp)
    ---------------------------------------------------------------------------
    Note: meshgrid's default indexing is 'xy', which means:
    lons, lats = np.meshgrid[lon, lat]
    in lons[i, j] or lats[i, j],  i->lat, j->lon
    ===========================================================================
    """
    def __init__(self, minlon, maxlon, dlon, minlat, maxlat, dlat, period, fieldtype='Tph'):
        self.Nlon=round((maxlon-minlon)/dlon)+1
        self.Nlat=round((maxlat-minlat)/dlat)+1
        self.dlon=dlon
        self.dlat=dlat
        self.lon=np.arange(self.Nlon)*self.dlon+minlon
        self.lat=np.arange(self.Nlat)*self.dlat+minlat
        self.lonArr, self.latArr = np.meshgrid(self.lon, self.lat)
        self.minlon=minlon
        self.maxlon=self.lon.max()
        self.minlat=minlat
        self.maxlat=self.lat.max()
        self._get_dlon_dlat_km()
        self.period=period
        self.fieldtype=fieldtype
        self.Zarr=np.zeros((self.Nlat, self.Nlon))
        return
    
    def copy(self):
        return copy.deepcopy(self)
    
    def _get_dlon_dlat_km(self):
        """Get longitude and latitude interval in km
        """
        self.dlon_km=np.array([])
        self.dlat_km=np.array([])
        for lat in self.lat:
            dist_lon, az, baz = obspy.geodetics.gps2dist_azimuth(lat, 0., lat, self.dlon)
            dist_lat, az, baz = obspy.geodetics.gps2dist_azimuth(lat, 0., lat+self.dlat, 0.)
            self.dlon_km=np.append(self.dlon_km, dist_lon/1000.)
            self.dlat_km=np.append(self.dlat_km, dist_lat/1000.)
        self.dlon_kmArr=(np.tile(self.dlon_km, self.Nlon).reshape(self.Nlon, self.Nlat)).T
        self.dlat_kmArr=(np.tile(self.dlat_km, self.Nlon).reshape(self.Nlon, self.Nlat)).T
        return
    
    def read(self, fname):
        """read field file
        """
        try:
            Inarray=np.loadtxt(fname)
            with open(fname) as f:
                inline = f.readline()
                if inline.split()[0] =='#':
                    evlostr = inline.split()[1]
                    evlastr = inline.split()[2]
                    if evlostr.split('=')[0] =='enx':
                        self.evlo = float(evlostr.split('=')[1])
                    if evlastr.split('=')[0] =='eny':
                        self.evla = float(evlastr.split('=')[1])
        except:
            Inarray=np.load(fname)
        self.lonArrIn=Inarray[:,0]
        self.latArrIn=Inarray[:,1]
        self.ZarrIn=Inarray[:,2]
        return
    
    def add_noise(self, sigma=0.5):
        """Add Gaussian noise with standard deviation = sigma to the input data
        """
        for i in xrange(self.ZarrIn.size):
            self.ZarrIn[i]=self.ZarrIn[i] + random.gauss(0, sigma)
        return
    
    def load_field(self, inField):
        """Load field data from an input object
        """
        self.lonArrIn=inField.lonArr
        self.latArrIn=inField.latArr
        self.ZarrIn=inField.Zarr
        return
    
    def write(self, fname, fmt='npy'):
        """Save field file
        """
        OutArr=np.append(self.lonArr, self.latArr)
        OutArr=np.append(OutArr, self.Zarr)
        OutArr=OutArr.reshape(3, self.Nlon*self.Nlat)
        OutArr=OutArr.T
        if fmt=='npy':
            np.save(fname, OutArr)
        elif fmt=='txt':
            np.savetxt(fname, OutArr)
        else:
            raise TypeError('Wrong output format!')
        return
    
    def write_dbase(self, outdir, gmt=False):
        """Write field data to an output directory
        """
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        outfdiffa=outdir+'/'+self.fieldtype+'_'+str(self.period)+'_diffa.lst'
        outffield=outdir+'/'+self.fieldtype+'_'+str(self.period)+'_field.lst'
        outfappV=outdir+'/'+self.fieldtype+'_'+str(self.period)+'_appV.lst'
        outflplc=outdir+'/'+self.fieldtype+'_'+str(self.period)+'_lplc.lst'
        outfreason=outdir+'/'+self.fieldtype+'_'+str(self.period)+'_reason.lst'
        outfgrad=outdir+'/'+self.fieldtype+'_'+str(self.period)+'_grad.lst'
        outfdist=outdir+'/'+self.fieldtype+'_'+str(self.period)+'_dist.lst'
        if gmt:
            lonArr=self.lonArr[self.reason_n==0]
            latArr=self.latArr[self.reason_n==0]
            self._write_txt(fname=outffield, outlon=lonArr, outlat=latArr, outZ=self.Zarr[self.reason_n==0])
            try:
                self._write_txt(fname=outfdiffa, outlon=lonArr, outlat=latArr, outZ=self.diffaArr[self.reason_n==0])
            except:
                pass
            try:
                self._write_txt(fname=outfappV, outlon=lonArr, outlat=latArr, outZ=self.appV[self.reason_n==0])
            except:
                pass
            try:
                self._write_txt(fname=outflplc, outlon=lonArr, outlat=latArr, outZ=self.lplc[self.reason_n==0])
            except:
                pass
            try:
                outArr=np.append(lonArr, latArr)
                outArr=np.append(outArr, self.grad[1][self.reason_n==0])
                outArr=np.append(outArr, self.grad[0][self.reason_n==0])
                outArr=outArr.reshape((4,latArr.size))
                outArr=outArr.T
                np.savetxt(outfgrad, outArr, fmt='%g')
            except:
                pass
            try:
                self._write_txt(fname=outfdist, outlon=lonArr, outlat=latArr, outZ=self.distArr[self.reason_n==0])
            except:
                pass
        else:
            lonArr=self.lonArr
            latArr=self.latArr
            self._write_txt(fname=outffield, outlon=lonArr, outlat=latArr, outZ=self.Zarr)
            try:
                self._write_txt(fname=outfdiffa, outlon=lonArr, outlat=latArr, outZ=self.diffaArr)
            except:
                pass
            try:
                self._write_txt(fname=outfappV, outlon=lonArr, outlat=latArr, outZ=self.appV)
            except:
                pass
            try:
                self._write_txt(fname=outflplc, outlon=lonArr, outlat=latArr, outZ=self.lplc)
            except:
                pass
            try:
                self._write_txt(fname=outfreason, outlon=lonArr, outlat=latArr, outZ=self.reason_n)
            except:
                pass
            try:
                outArr=np.append(lonArr, latArr)
                outArr=np.append(outArr, self.grad[1])
                outArr=np.append(outArr, self.grad[0])
                outArr=outArr.reshape((4,latArr.size))
                outArr=outArr.T
                np.savetxt(outfgrad, outArr, fmt='%g')
            except:
                pass
            try:
                self._write_txt(fname=outfdist, outlon=lonArr, outlat=latArr, outZ=self.distArr)
            except:
                pass
        return
    
    def read_dbase(self, datadir):
        """Read field data from an input directory
        """
        infdiffa=datadir+'/'+self.fieldtype+'_'+str(self.period)+'_diffa.lst'
        inffield=datadir+'/'+self.fieldtype+'_'+str(self.period)+'_field.lst'
        infappV=datadir+'/'+self.fieldtype+'_'+str(self.period)+'_appV.lst'
        inflplc=datadir+'/'+self.fieldtype+'_'+str(self.period)+'_lplc.lst'
        infreason=datadir+'/'+self.fieldtype+'_'+str(self.period)+'_reason.lst'
        infgrad=datadir+'/'+self.fieldtype+'_'+str(self.period)+'_grad.lst'
        infdist=datadir+'/'+self.fieldtype+'_'+str(self.period)+'_dist.lst'
        inArr=np.loadtxt(inffield)
        Zarr=inArr[:,2]
        if Zarr.size==(self.Nlat-1)*(self.Nlon-1):
            self.cut_edge(1,1)
        elif Zarr.size==(self.Nlat-4)*(self.Nlon-4):
            self.cut_edge(2,2)
        elif Zarr.size!=self.Nlat*self.Nlon:
            raise ValueError('Incompatible Nlon/Nlat with input data base!')
        self.Zarr=Zarr.reshape(self.Nlat, self.Nlon)
        if os.path.isfile(infdiffa):
            inArr=np.loadtxt(infdiffa)
            self.diffaArr=inArr[:,2]
            self.diffaArr=self.diffaArr.reshape(self.Nlat, self.Nlon)
        if os.path.isfile(infappV):
            inArr=np.loadtxt(infappV)
            self.appV=inArr[:,2]
            self.appV=self.appV.reshape(self.Nlat, self.Nlon)
        if os.path.isfile(inflplc):
            inArr=np.loadtxt(inflplc)
            self.lplc=inArr[:,2]
            self.lplc=self.lplc.reshape(self.Nlat, self.Nlon)
        if os.path.isfile(infreason):
            inArr=np.loadtxt(infreason)
            self.reason_n=inArr[:,2]
            self.reason_n=self.reason_n.reshape(self.Nlat, self.Nlon)
        if os.path.isfile(infgrad):
            inArr=np.loadtxt(infgrad)
            self.grad=[]
            self.grad.append(inArr[:,3])
            self.grad.append(inArr[:,2])
            self.grad[0]=self.grad[0].reshape(self.Nlat, self.Nlon)
            self.grad[1]=self.grad[1].reshape(self.Nlat, self.Nlon)
        if os.path.isfile(infdist):
            inArr=np.loadtxt(infdist)
            self.distArr=inArr[:,2]
            self.distArr=self.distArr.reshape(self.Nlat, self.Nlon)
        
        return
        
    
    def _write_txt(self, fname, outlon, outlat, outZ):
        outArr=np.append(outlon, outlat)
        outArr=np.append(outArr, outZ)
        outArr=outArr.reshape((3,outZ.size))
        outArr=outArr.T
        np.savetxt(fname, outArr, fmt='%g')
        return
    
    def np2ma(self):
        """Convert all the data array to masked array according to reason_n array.
        """
        try:
            reason_n=self.reason_n
        except:
            raise AttrictError('No reason_n array!')
        self.Zarr=ma.masked_array(self.Zarr, mask=np.zeros(reason_n.shape) )
        self.Zarr.mask[reason_n!=0]=1
        try:
            self.diffaArr=ma.masked_array(self.diffaArr, mask=np.zeros(reason_n.shape) )
            self.diffaArr.mask[reason_n!=0]=1
        except:
            pass
        try:
            self.appV=ma.masked_array(self.appV, mask=np.zeros(reason_n.shape) )
            self.appV.mask[reason_n!=0]=1
        except:
            pass
        try:
            self.grad[0]=ma.masked_array(self.grad[0], mask=np.zeros(reason_n.shape) )
            self.grad[0].mask[reason_n!=0]=1
            self.grad[1]=ma.masked_array(self.grad[1], mask=np.zeros(reason_n.shape) )
            self.grad[1].mask[reason_n!=0]=1
        except:
            pass
        try:
            self.lplc=ma.masked_array(self.lplc, mask=np.zeros(reason_n[1:-1, 1:-1].shape) )
            self.lplc.mask[reason_n[1:-1, 1:-1]!=0]=1
        except:
            print 'No Laplacian array!'
            pass
        return
    
    def ma2np(self):
        """Convert all the maksed data array to numpy array
        """
        self.Zarr=ma.getdata(self.Zarr)
        try:
            self.diffaArr=ma.getdata(self.diffaArr)
        except:
            pass
        try:
            self.appV=ma.getdata(self.appV)
        except:
            pass
        try:
            self.lplc=ma.getdata(self.lplc)
        except:
            pass
        return
    
    def cut_edge(self, nlon, nlat):
        """Cut edge
        =======================================================================================
        Input Parameters:
        nlon, nlon  - number of edge point in longitude/latitude to be cutted
        =======================================================================================
        """
        self.Nlon=self.Nlon-2*nlon
        self.Nlat=self.Nlat-2*nlat
        self.minlon=self.minlon + nlon*self.dlon
        self.maxlon=self.maxlon - nlon*self.dlon
        self.minlat=self.minlat + nlat*self.dlat
        self.maxlat=self.maxlat - nlat*self.dlat
        self.lon=np.arange(self.Nlon)*self.dlon+self.minlon
        self.lat=np.arange(self.Nlat)*self.dlat+self.minlat
        self.lonArr, self.latArr = np.meshgrid(self.lon, self.lat)
        self.Zarr=self.Zarr[nlat:-nlat, nlon:-nlon]
        try:
            self.reason_n=self.reason_n[nlat:-nlat, nlon:-nlon]
        except:
            pass
        self._get_dlon_dlat_km()
        return
    
    def gradient(self, method='default', edge_order=1, order=2):
        """Compute gradient of the field
        =============================================================================================================
        Input Parameters:
        edge_order  - edge_order : {1, 2}, optional, only has effect when method='default'
                        Gradient is calculated using Nth order accurate differences at the boundaries
        method      - method: 'default' : use numpy.gradient 'convolve': use convolution
        order       - order of finite difference scheme, only has effect when method='convolve'
        =============================================================================================================
        """
        Zarr=self.Zarr
        if method=='default':
            # self.dlat_kmArr : dx here in numpy gradient since Zarr is Z[ilat, ilon]
            self.grad=np.gradient( self.Zarr, self.dlat_kmArr, self.dlon_kmArr, edge_order=edge_order)
            self.grad[0]=self.grad[0][1:-1, 1:-1]
            self.grad[1]=self.grad[1][1:-1, 1:-1]
        elif method=='fft':
            diff_lon=self.fftDiff2(m=1, n=0)
            diff_lat=self.fftDiff2(m=0, n=1)
            self.grad=[]
            self.grad.append(diff_lat)
            self.grad.append(diff_lon)
        elif method == 'convolve':
            dlat_km=self.dlat_kmArr
            dlon_km=self.dlon_kmArr
            if order==2:
                diff_lon=convolve(Zarr, lon_diff_weight_2)/dlon_km
                diff_lat=convolve(Zarr, lat_diff_weight_2)/dlat_km
            elif order==4:
                diff_lon=convolve(Zarr, lon_diff_weight_4)/dlon_km
                diff_lat=convolve(Zarr, lat_diff_weight_4)/dlat_km
            elif order==6:
                diff_lon=convolve(Zarr, lon_diff_weight_6)/dlon_km
                diff_lat=convolve(Zarr, lat_diff_weight_6)/dlat_km
            self.grad=[]
            self.grad.append(diff_lat[1:-1, 1:-1])
            self.grad.append(diff_lon[1:-1, 1:-1])
        self.proAngle=np.arctan2(self.grad[0], self.grad[1])/np.pi*180.
        return
    
    def get_appV(self):
        """Get the apparent velocity from gradient
        """
        slowness=np.sqrt ( self.grad[0] ** 2 + self.grad[1] ** 2)
        slowness[slowness==0]=0.3
        self.appV = 1./slowness
        return
    
    def Laplacian(self, method='convolve', order=4, verbose=False):
        """Compute Laplacian of the field
        =============================================================================================================
        Input Parameters:
        edge_order  - edge_order : {1, 2}, optional, only has effect when method='default'
                        Gradient is calculated using Nth order accurate differences at the boundaries
        method      - method: 'default' : use numpy.gradient 'convolve': use convolution
        order       - order of finite difference scheme, only has effect when method='convolve'
        =============================================================================================================
        """
        Zarr=self.Zarr
        if method == 'default':
            dlat_km=self.dlat_kmArr[1:-1, 1:-1]
            dlon_km=self.dlon_kmArr[1:-1, 1:-1]
            Zarr_latp=Zarr[2:, 1:-1]
            Zarr_latn=Zarr[:-2, 1:-1]
            Zarr_lonp=Zarr[1:-1, 2:]
            Zarr_lonn=Zarr[1:-1, :-2]
            Zarr=Zarr[1:-1, 1:-1]
            self.lplc=(Zarr_latp+Zarr_latn-2*Zarr) / (dlat_km**2) + (Zarr_lonp+Zarr_lonn-2*Zarr) / (dlon_km**2)
        elif method == 'convolve':
            dlat_km=self.dlat_kmArr
            dlon_km=self.dlon_kmArr
            if order==2:
                diff2_lon=convolve(Zarr, lon_diff2_weight_2)/dlon_km/dlon_km
                diff2_lat=convolve(Zarr, lat_diff2_weight_2)/dlat_km/dlat_km
            elif order==4:
                diff2_lon=convolve(Zarr, lon_diff2_weight_4)/dlon_km/dlon_km
                diff2_lat=convolve(Zarr, lat_diff2_weight_4)/dlat_km/dlat_km
            elif order==6:
                diff2_lon=convolve(Zarr, lon_diff2_weight_6)/dlon_km/dlon_km
                diff2_lat=convolve(Zarr, lat_diff2_weight_6)/dlat_km/dlat_km
            self.lplc=diff2_lon+diff2_lat
            self.lplc=self.lplc[1:-1, 1:-1]
        if verbose:
            print 'max lplc:',self.lplc.max(), 'min lplc:',self.lplc.min()
        return
    
    def Laplacian_Green(self):
        """Compute Laplacian of the field using Green's theorem( 2D Gauss's theorem )
        """
        try:
            grad_y=self.grad[0]; grad_x=self.grad[1]
        except:
            self.gradient('default')
            self.cut_edge(1,1)
            grad_y=self.grad[0]; grad_x=self.grad[1]
        grad_xp=grad_x[1:-1, 2:];  grad_xn=grad_x[1:-1, :-2]
        grad_yp=grad_y[2:, 1:-1];  grad_yn=grad_y[:-2, 1:-1]
        dlat_km=self.dlat_kmArr[1:-1, 1:-1]; dlon_km=self.dlon_kmArr[1:-1, 1:-1]
        loopsum=(grad_xp-grad_xn)*dlat_km + (grad_yp-grad_yn)*dlon_km
        area=dlat_km*dlon_km
        lplc = loopsum/area
        # # lplc=np.zeros((self.Nlat-2, self.Nlon-2))
        # # for ilon in np.arange(int(self.Nlon-2))+1:
        # #     for ilat in np.arange(int(self.Nlat-2))+1:
        # #         dlat_km=self.dlat_kmArr[ilat, ilon]; dlon_km=self.dlon_kmArr[ilat, ilon]
        # #         area=dlat_km*dlon_km
        # #         loopsum=grad_x[ilat, ilon+1]*dlat_km-grad_x[ilat, ilon-1]*dlat_km \
        # #                 + grad_y[ilat+1, ilon]*dlon_km - grad_y[ilat-1, ilon]*dlon_km
        # #         lplc[ilat-1, ilon-1] = loopsum/area
        self.lplc=lplc
        return
    
    def interp_surface(self, workingdir, outfname, tension=0.0):
        """Interpolate input data to grid point with gmt surface command
        =======================================================================================
        Input Parameters:
        workingdir  - working directory
        outfname    - output file name for interpolation
        tension     - input tension for gmt surface(0.0-1.0)
        =======================================================================================
        """
        if not os.path.isdir(workingdir):
            os.makedirs(workingdir)
        OutArr=np.append(self.lonArrIn, self.latArrIn)
        OutArr=np.append(OutArr, self.ZarrIn)
        OutArr=OutArr.reshape(3, self.lonArrIn.size)
        OutArr=OutArr.T
        np.savetxt(workingdir+'/'+outfname, OutArr, fmt='%g')
        fnameHD=workingdir+'/'+outfname+'.HD'
        tempGMT=workingdir+'/'+outfname+'_GMT.sh'
        grdfile=workingdir+'/'+outfname+'.grd'
        with open(tempGMT,'wb') as f:
            REG='-R'+str(self.minlon)+'/'+str(self.maxlon)+'/'+str(self.minlat)+'/'+str(self.maxlat)
            f.writelines('gmtset MAP_FRAME_TYPE fancy \n')
            f.writelines('surface %s -T%g -G%s -I%g %s \n' %( workingdir+'/'+outfname, tension, grdfile, self.dlon, REG ))
            f.writelines('grd2xyz %s %s > %s \n' %( grdfile, REG, fnameHD ))
        call(['bash', tempGMT])
        os.remove(grdfile)
        os.remove(tempGMT)
        Inarray=np.loadtxt(fnameHD)
        ZarrIn=Inarray[:,2]
        self.Zarr=(ZarrIn.reshape(self.Nlat, self.Nlon))[::-1, :]
        return
    
    def check_curvature(self, workingdir, threshold=0.005):
        """
        Check and discard those points with large curvatures.
        Points at boundaries will be discarded.
        Two interpolation schemes with different tension (0, 0.2) will be applied to the quality controlled field data file. 
        =====================================================================================================================
        Input parameters:
        workingdir  - working directory
        threshold   - threshold value for Laplacian
        ---------------------------------------------------------------------------------------------------------------------
        Output format:
        workingdir/fieldtype_per_v1.lst         - output field file with data point passing curvature checking
        workingdir/fieldtype_per_v1.lst.HD      - interpolated travel time file 
        workingdir/fieldtype_per_v1.lst.HD_0.2  - interpolated travel time file with tension=0.2
        =====================================================================================================================
        """
        # Compute Laplacian
        self.Laplacian(method='convolve', order=4)
        self.cut_edge(1,1)
        # quality control
        LonLst=self.lonArr.reshape(self.lonArr.size)
        LatLst=self.latArr.reshape(self.latArr.size)
        TLst=self.Zarr.reshape(self.Zarr.size)
        lplc = self.lplc.reshape(self.lplc.size)
        index = np.where((lplc>-threshold)*(lplc<threshold))[0]
        LonLst=LonLst[index]
        LatLst=LatLst[index]
        TLst=TLst[index]
        # output to txt file
        outfname=workingdir+'/'+self.fieldtype+'_'+str(self.period)+'_v1.lst'
        TfnameHD=outfname+'.HD'
        self._write_txt(fname=outfname, outlon=LonLst, outlat=LatLst, outZ=TLst)
        # interpolate with gmt surface
        tempGMT=workingdir+'/'+self.fieldtype+'_'+str(self.period)+'_v1_GMT.sh'
        grdfile=workingdir+'/'+self.fieldtype+'_'+str(self.period)+'_v1.grd'
        with open(tempGMT,'wb') as f:
            REG='-R'+str(self.minlon)+'/'+str(self.maxlon)+'/'+str(self.minlat)+'/'+str(self.maxlat)
            f.writelines('gmtset MAP_FRAME_TYPE fancy \n')
            f.writelines('surface %s -T0.0 -G%s -I%g %s \n' %( outfname, grdfile, self.dlon, REG ))
            f.writelines('grd2xyz %s %s > %s \n' %( grdfile, REG, TfnameHD ))
            f.writelines('surface %s -T0.2 -G%s -I%g %s \n' %( outfname, grdfile+'.T0.2', self.dlon, REG ))
            f.writelines('grd2xyz %s %s > %s \n' %( grdfile+'.T0.2', REG, TfnameHD+'_0.2' ))
        call(['bash', tempGMT])
        os.remove(grdfile+'.T0.2')
        os.remove(grdfile)
        os.remove(tempGMT)
        return 
        
    def gradient_qc(self, workingdir, evlo, evla, nearneighbor=True, cdist=None):
        """
        Generate Slowness Maps from Travel Time Maps.
        Two interpolated travel time file with different tension will be used for quality control.
        =====================================================================================================================
        Input parameters:
        workingdir      - working directory
        evlo, evla      - event location
        nearneighbor    - do near neighbor quality control or not
        cdist           - distance for quality control, default is 12*period
        
        Output format:
        outdir/slow_azi_stacode.pflag.txt.HD.2.v2 - Slowness map
        =====================================================================================================================
        """
        if cdist==None:
            cdist=12.*self.period
        # Read data,
        # v1: data that pass check_curvature criterion
        # v1HD and v1HD02: interpolated v1 data with tension = 0. and 0.2
        fnamev1=workingdir+'/'+self.fieldtype+'_'+str(self.period)+'_v1.lst'
        fnamev1HD=fnamev1+'.HD'
        fnamev1HD02=fnamev1HD+'_0.2'
        InarrayV1=np.loadtxt(fnamev1)
        loninV1=InarrayV1[:,0]
        latinV1=InarrayV1[:,1]
        fieldin=InarrayV1[:,2]
        Inv1HD=np.loadtxt(fnamev1HD)
        lonv1HD=Inv1HD[:,0]
        latv1HD=Inv1HD[:,1]
        fieldv1HD=Inv1HD[:,2]
        Inv1HD02=np.loadtxt(fnamev1HD02)
        lonv1HD02=Inv1HD02[:,0]
        latv1HD02=Inv1HD02[:,1]
        fieldv1HD02=Inv1HD02[:,2]
        # Set field value to be zero if there is large difference between v1HD and v1HD02
        diffArr = fieldv1HD-fieldv1HD02
        fieldArr=fieldv1HD*((diffArr<2.)*(diffArr>-2.)) 
        fieldArr = (fieldArr.reshape(self.Nlat, self.Nlon))[::-1,:]
        # reason_n -> 0: accepted point 1: data point the has large difference between v1HD and v1HD02
        # 2: data point that does not have near neighbor points at all E/W/N/S directions
        reason_n=np.ones(fieldArr.size)
        reason_n1=reason_n*(diffArr>2.)
        reason_n2=reason_n*(diffArr<-2.)
        reason_n=reason_n1+reason_n2
        reason_n = (reason_n.reshape(self.Nlat, self.Nlon))[::-1,:]
        # Nested loop, may need modification to speed the code up
        if nearneighbor:
            print 'Start near neighbor quality control checking'
            for ilat in xrange(self.Nlat):
                for ilon in xrange(self.Nlon):
                    if reason_n[ilat, ilon]==1:
                        continue
                    lon=self.lon[ilon]
                    lat=self.lat[ilat]
                    dlon_km=self.dlon_km[ilat]
                    dlat_km=self.dlat_km[ilat]
                    difflon=abs(loninV1-lon)/self.dlon*dlon_km
                    difflat=abs(latinV1-lat)/self.dlat*dlat_km
                    index=np.where((difflon<cdist)*(difflat<cdist))[0]
                    marker_EN=np.zeros((2,2))
                    marker_nn=4
                    tflag = False
                    for iv1 in index:
                        lon2=loninV1[iv1]
                        lat2=latinV1[iv1]
                        if lon2-lon<0:
                            marker_E=0
                        else:
                            marker_E=1
                        if lat2-lat<0:
                            marker_N=0
                        else:
                            marker_N=1
                        if marker_EN[marker_E , marker_N]==1:
                            continue
                        dist, az, baz=obspy.geodetics.gps2dist_azimuth(lat,lon,lat2,lon2)
                        dist=dist/1000.
                        if dist< cdist*2 and dist >= 1:
                            marker_nn=marker_nn-1
                            if marker_nn==0:
                                tflag = True
                                break
                            marker_EN[marker_E, marker_N]=1
                    if tflag==False:
                        fieldArr[ilat, ilon]=0
                        reason_n[ilat, ilon] = 2
            print 'End near neighbor quality control checking'
        # Start to Compute Gradient
        self.Zarr=fieldArr
        self.gradient('default')
        self.cut_edge(1,1)
        # if one field point has zero value, reason_n for four near neighbor points will all be set to 4
        index0=np.where(self.Zarr==0)
        ilatArr=index0[0]+1
        ilonArr=index0[1]+1
        reason_n[ilatArr+1, ilonArr]=4
        reason_n[ilatArr-1, ilonArr]=4
        reason_n[ilatArr, ilonArr+1]=4
        reason_n[ilatArr, ilonArr-1]=4
        reason_n=reason_n[1:-1,1:-1]
        # if slowness is too large/small, reason_n will be set to 3
        slowness=np.sqrt(self.grad[0]**2+self.grad[1]**2)
        # reason_n[slowness>100]=3
        if self.fieldtype=='Tph' or self.fieldtype=='Tgr':
            reason_n[slowness>0.6]=3
            reason_n[slowness<0.2]=3
        print 'Computing deflections'
        indexvalid=np.where(reason_n==0)
        diffaArr=np.zeros(reason_n.shape)
        latinArr=self.lat[indexvalid[0]]
        loninArr=self.lon[indexvalid[1]]
        evloArr=np.ones(loninArr.size)*evlo
        evlaArr=np.ones(loninArr.size)*evla
        g = Geod(ellps='WGS84')
        az, baz, distevent = g.inv(loninArr, latinArr, evloArr, evlaArr)
        distevent=distevent/1000.
        az = az + 180.
        az = 90.-az
        baz = 90.-baz
        az[az>180.]=az[az>180.] - 360.
        az[az<-180.]=az[az<-180.] + 360.
        baz[baz>180.]=baz[baz>180.] - 360.
        baz[baz<-180.]=baz[baz<-180.] + 360.
        diffaArr[indexvalid[0], indexvalid[1]] = \
            self.proAngle[indexvalid[0], indexvalid[1]] - az
        # if epicentral distance is too small, reason_n will be set to 5, and diffaArr will be 0.
        tempArr = diffaArr[indexvalid[0], indexvalid[1]]
        tempArr[distevent<cdist+50.] = 0.
        diffaArr[indexvalid[0], indexvalid[1]] = tempArr
        diffaArr[diffaArr>180.]=diffaArr[diffaArr>180.]-360.
        diffaArr[diffaArr<-180.]=diffaArr[diffaArr<-180.]+360.
        tempArr = reason_n[indexvalid[0], indexvalid[1]]
        tempArr[distevent<cdist+50.] = 5
        reason_n[indexvalid[0], indexvalid[1]] = tempArr
        self.diffaArr=diffaArr
        self.get_appV()
        self.reason_n=reason_n
        try:
            self.lplc=self.lplc[1:-1, 1:-1]
        except:
            pass
        return 
    
    def _get_basemap(self, projection='lambert', geopolygons=None, resolution='i'):
        """Plot data with contour
        """
        # fig=plt.figure(num=None, figsize=(12, 12), dpi=80, facecolor='w', edgecolor='k')
        lat_centre = (self.maxlat+self.minlat)/2.0
        lon_centre = (self.maxlon+self.minlon)/2.0
        if projection=='merc':
            m=Basemap(projection='merc', llcrnrlat=self.minlat-5., urcrnrlat=self.maxlat+5., llcrnrlon=self.minlon-5.,
                      urcrnrlon=self.maxlon+5., lat_ts=20, resolution=resolution)
            # m.drawparallels(np.arange(self.minlat,self.maxlat,self.dlat), labels=[1,0,0,1])
            # m.drawmeridians(np.arange(self.minlon,self.maxlon,self.dlon), labels=[1,0,0,1])
            m.drawparallels(np.arange(-80.0,80.0,5.0), labels=[1,0,0,1])
            m.drawmeridians(np.arange(-170.0,170.0,5.0), labels=[1,0,0,1])
            m.drawstates(color='g', linewidth=2.)
        elif projection=='global':
            m=Basemap(projection='ortho',lon_0=lon_centre, lat_0=lat_centre, resolution=resolution)
            m.drawparallels(np.arange(-80.0,80.0,10.0), labels=[1,0,0,1])
            m.drawmeridians(np.arange(-170.0,170.0,10.0), labels=[1,0,0,1])
        
        elif projection=='regional_ortho':
            m1 = Basemap(projection='ortho', lon_0=self.minlon, lat_0=self.minlat, resolution='l')
            m = Basemap(projection='ortho', lon_0=self.minlon, lat_0=self.minlat, resolution=resolution,\
                llcrnrx=0., llcrnry=0., urcrnrx=m1.urcrnrx/mapfactor, urcrnry=m1.urcrnry/3.5)
            m.drawparallels(np.arange(-80.0,80.0,10.0), labels=[1,0,0,0],  linewidth=2,  fontsize=20)
            # m.drawparallels(np.arange(-90.0,90.0,30.0),labels=[1,0,0,0], dashes=[10, 5], linewidth=2,  fontsize=20)
            # m.drawmeridians(np.arange(10,180.0,30.0), dashes=[10, 5], linewidth=2)
            m.drawmeridians(np.arange(-170.0,170.0,10.0),  linewidth=2)
        elif projection=='lambert':
            distEW, az, baz=obspy.geodetics.gps2dist_azimuth(self.minlat, self.minlon,
                                self.minlat, self.maxlon) # distance is in m
            distNS, az, baz=obspy.geodetics.gps2dist_azimuth(self.minlat, self.minlon,
                                self.maxlat+2., self.minlon) # distance is in m
            m = Basemap(width=distEW, height=distNS, rsphere=(6378137.00,6356752.3142), resolution='l', projection='lcc',\
                lat_1=self.minlat, lat_2=self.maxlat, lon_0=lon_centre, lat_0=lat_centre+1)
            m.drawparallels(np.arange(-80.0,80.0,10.0), linewidth=1, dashes=[2,2], labels=[1,1,0,0], fontsize=15)
            m.drawmeridians(np.arange(-170.0,170.0,10.0), linewidth=1, dashes=[2,2], labels=[0,0,1,0], fontsize=15)
            # m.drawparallels(np.arange(-80.0,80.0,10.0), linewidth=0.5, dashes=[2,2], labels=[1,0,0,0], fontsize=5)
            # m.drawmeridians(np.arange(-170.0,170.0,10.0), linewidth=0.5, dashes=[2,2], labels=[0,0,0,1], fontsize=5)
        m.drawcoastlines(linewidth=1.0)
        m.drawcountries(linewidth=1.)
        
        # m.drawmapboundary(fill_color=[1.0,1.0,1.0])
        m.fillcontinents(lake_color='#99ffff',zorder=0.2)
        # m.drawlsmask(land_color='0.8', ocean_color='#99ffff')
        m.drawmapboundary(fill_color="white")
        try:
            geopolygons.PlotPolygon(inbasemap=m)
        except:
            pass
        return m
    
    def _get_basemap_1(self, projection='lambert', geopolygons=None, resolution='i'):
        """Plot data with contour
        """
        # fig=plt.figure(num=None, figsize=(12, 12), dpi=80, facecolor='w', edgecolor='k')
        lat_centre = (self.maxlat+self.minlat)/2.0
        lon_centre = (self.maxlon+self.minlon)/2.0
        if projection=='merc':
            m=Basemap(projection='merc', llcrnrlat=self.minlat-5., urcrnrlat=self.maxlat+5., llcrnrlon=self.minlon-5.,
                      urcrnrlon=self.maxlon+5., lat_ts=20, resolution=resolution)
            # m.drawparallels(np.arange(self.minlat,self.maxlat,self.dlat), labels=[1,0,0,1])
            # m.drawmeridians(np.arange(self.minlon,self.maxlon,self.dlon), labels=[1,0,0,1])
            m.drawparallels(np.arange(-80.0,80.0,5.0), labels=[1,0,0,1])
            m.drawmeridians(np.arange(-170.0,170.0,5.0), labels=[1,0,0,1])
            m.drawstates(color='g', linewidth=2.)
        elif projection=='global':
            m=Basemap(projection='ortho',lon_0=lon_centre, lat_0=lat_centre, resolution=resolution)
            m.drawparallels(np.arange(-80.0,80.0,10.0), labels=[1,0,0,1])
            m.drawmeridians(np.arange(-170.0,170.0,10.0), labels=[1,0,0,1])
        
        elif projection=='regional_ortho':
            m1 = Basemap(projection='ortho', lon_0=self.minlon, lat_0=self.minlat, resolution='l')
            m = Basemap(projection='ortho', lon_0=self.minlon, lat_0=self.minlat, resolution=resolution,\
                llcrnrx=0., llcrnry=0., urcrnrx=m1.urcrnrx/mapfactor, urcrnry=m1.urcrnry/3.5)
            m.drawparallels(np.arange(-80.0,80.0,10.0), labels=[1,0,0,0],  linewidth=2,  fontsize=20)
            # m.drawparallels(np.arange(-90.0,90.0,30.0),labels=[1,0,0,0], dashes=[10, 5], linewidth=2,  fontsize=20)
            # m.drawmeridians(np.arange(10,180.0,30.0), dashes=[10, 5], linewidth=2)
            m.drawmeridians(np.arange(-170.0,170.0,10.0),  linewidth=2)
        elif projection=='lambert':
            distEW, az, baz=obspy.geodetics.gps2dist_azimuth(self.minlat, self.minlon,
                                self.minlat, self.maxlon) # distance is in m
            distNS, az, baz=obspy.geodetics.gps2dist_azimuth(self.minlat, self.minlon,
                                self.maxlat+2, self.minlon) # distance is in m
            m = Basemap(width=distEW, height=distNS, rsphere=(6378137.00,6356752.3142), resolution='l', projection='lcc',\
                lat_1=self.minlat, lat_2=self.maxlat, lon_0=lon_centre, lat_0=lat_centre+1.5)
            # m.drawparallels(np.arange(-80.0,80.0,10.0), linewidth=2, dashes=[2,2], labels=[1,0,0,0], fontsize=15)
            # m.drawmeridians(np.arange(-170.0,170.0,10.0), linewidth=2, dashes=[2,2], labels=[0,0,1,1], fontsize=15)
            m.drawparallels(np.arange(-80.0,80.0,10.0), linewidth=0.5, dashes=[2,2], labels=[0,0,0,0], fontsize=5)
            m.drawmeridians(np.arange(-170.0,170.0,10.0), linewidth=0.5, dashes=[2,2], labels=[0,0,0,1], fontsize=5)
        m.drawcoastlines(linewidth=1.0)
        m.drawcountries(linewidth=1.)
        
        # m.drawmapboundary(fill_color=[1.0,1.0,1.0])
        m.fillcontinents(lake_color='#99ffff',zorder=0.2)
        # m.drawlsmask(land_color='0.8', ocean_color='#99ffff')
        m.drawmapboundary(fill_color="white")
        try:
            geopolygons.PlotPolygon(inbasemap=m)
        except:
            pass
        return m
    

    def plot_field(self, projection='lambert', contour=True, geopolygons=None, showfig=True, vmin=None, vmax=None):
        """Plot data with contour
        """
        m=self._get_basemap(projection=projection, geopolygons=geopolygons)
        x, y=m(self.lonArr, self.latArr)
        if self.fieldtype=='Ms':
            cmap=pycpt.load.gmtColormap('./GMT_panoply.cpt')
            cmap =discrete_cmap(int((vmax-vmin)/0.1)+1, cmap)
            im=m.pcolormesh(x, y, self.Zarr, cmap=cmap, shading='gouraud', vmin=vmin, vmax=vmax)
            cb = m.colorbar(im, "bottom", size="3%", pad='2%', ticks=np.arange(20)*0.1+2.5)
        else:
            im=m.pcolormesh(x, y, self.Zarr, cmap='gist_ncar_r', shading='gouraud', vmin=vmin, vmax=vmax)
            cb = m.colorbar(im, "bottom", size="3%", pad='2%', ticks=np.arange(17)*100.)
        cb.ax.tick_params(labelsize=10)
        if self.fieldtype=='Tph' or self.fieldtype=='Tgr':
            cb.set_label('Travel Time (sec)', fontsize=15, rotation=0)
        if self.fieldtype=='Amp':
            cb.set_label('nm', fontsize=12, rotation=0)
        if contour:
            # levels=np.linspace(ma.getdata(self.Zarr).min(), ma.getdata(self.Zarr).max(), 20)
            levels=np.linspace(ma.getdata(self.Zarr).min(), ma.getdata(self.Zarr).max(), 60)
            m.contour(x, y, self.Zarr, colors='k', levels=levels, linewidths=0.5)
        if showfig:
            plt.show()
        return
    
    def plot_lplc(self, projection='lambert', contour=False, geopolygons=None, vmin=None, vmax=None, showfig=True):
        """Plot data with contour
        """
        m=self._get_basemap(projection=projection, geopolygons=geopolygons)
        if self.lonArr.shape[0]-2==self.lplc.shape[0] and self.lonArr.shape[1]-2==self.lplc.shape[1]:
            self.cut_edge(1,1)
        elif self.lonArr.shape[0]!=self.lplc.shape[0] or self.lonArr.shape[1]!=self.lplc.shape[1]:
            raise ValueError('Incompatible shape for lplc and lon/lat array!')
        x, y=m(self.lonArr, self.latArr)
        # cmap =discrete_cmap(int(vmax-vmin)/2+1, 'seismic')
        m.pcolormesh(x, y, self.lplc, cmap='seismic', shading='gouraud', vmin=vmin, vmax=vmax)
        cb=m.colorbar()
        cb.ax.tick_params(labelsize=15) 
        levels=np.linspace(self.lplc.min(), self.lplc.max(), 100)
        if contour:
            plt.contour(x, y, self.lplc, colors='k', levels=levels)
        if showfig:
            plt.show()
        return
    
    def plot_lplcC(self, infield=None, projection='lambert', contour=False, geopolygons=None, vmin=-0.012, vmax=0.012, period=10., showfig=True):
        """Plot data with contour
        """
        m=self._get_basemap(projection=projection, geopolygons=geopolygons)
        if self.lonArr.shape[0]-2==self.lplc.shape[0] and self.lonArr.shape[1]-2==self.lplc.shape[1]:
            self.cut_edge(1,1)
        elif self.lonArr.shape[0]!=self.lplc.shape[0] or self.lonArr.shape[1]!=self.lplc.shape[1]:
            raise ValueError('Incompatible shape for lplc and lon/lat array!')
        w=2*np.pi/period
        Zarr=self.Zarr.copy()
        Zarr[self.Zarr==0]=-1
        lplcC=self.lplc/Zarr/w**2
        if infield!=None:
            lplcC=lplcC*(infield.appV[1:-1,1:-1])**3/2.
            vmin=-0.2
            vmax=0.2
            # appV=infield.appV[1:-1,1:-1]
            # lplcC=1./np.sqrt(1./appV**2 - lplcC)
            # vmin=2.9
            # vmax=3.4
        lplcC=ma.masked_array(lplcC, mask=np.zeros(self.Zarr.shape) )
        lplcC.mask[self.reason_n!=0]=1
        x, y=m(self.lonArr, self.latArr)
        cmap =discrete_cmap(int((vmax-vmin)*80)/2+1, 'seismic')
        # cmap =discrete_cmap(int((vmax-vmin)*1000)/2+1, 'seismic')
        # cmap=pycpt.load.gmtColormap('./GMT_panoply.cpt')
        # cmap =discrete_cmap(int((vmax-vmin)*1000), 'seismic'2
        # cmap = colors.get_colormap('tomo_80_perc_linear_lightness')
        im=m.pcolormesh(x, y, lplcC, cmap=cmap, shading='gouraud', vmin=vmin, vmax=vmax)
        cb = m.colorbar(im, "right", size="3%", pad='2%')
        cb.ax.tick_params(labelsize=5)
        cb.set_label(r"$\frac{\mathrm{km}}{\mathrm{s}}$", fontsize=8, rotation=0)
        if showfig:
            plt.show()
        return
    
    def plot_diffa(self, projection='lambert', factor=3, prop=True, geopolygons=None, cmap='seismic', vmin=-20, vmax=20, showfig=True, rayASDF='', rayfactor=30):
        """Plot data with contour
        """
        m=self._get_basemap(projection=projection, geopolygons=geopolygons)
        if self.lonArr.shape[0]-2==self.diffaArr.shape[0] and self.lonArr.shape[1]-2==self.diffaArr.shape[1]:
            self.cut_edge(1,1)
        elif self.lonArr.shape[0]!=self.diffaArr.shape[0] or self.lonArr.shape[1]!=self.diffaArr.shape[1]:
            raise ValueError('Incompatible shape for deflection and lon/lat array!')
        x, y=m(self.lonArr, self.latArr)
        cmap=pycpt.load.gmtColormap('./GMT_panoply.cpt')
        cmap='bwr'
        cmap =discrete_cmap(int(vmax-vmin)/4, cmap)
        im=m.pcolormesh(x, y, self.diffaArr, cmap=cmap, shading='gouraud', vmin=vmin, vmax=vmax)
        cb = m.colorbar(im, "bottom", size="3%", pad='2%')
        cb.ax.tick_params(labelsize=10)
        cb.set_label('degree', fontsize=12, rotation=0)
        if os.path.isfile(rayASDF):
            try:
                dset=pyasdf.ASDFDataSet(rayASDF)
                dlon=self.dlon*rayfactor
                dlat=self.dlat*rayfactor
                Nlon=round((self.maxlon-self.minlon)/dlon)+1
                Nlat=round((self.maxlat-self.minlat)/dlat)+1
                lons=np.arange(Nlon)*dlon+self.minlon
                lats=np.arange(Nlat)*dlat+self.minlat
                for lon in lons:
                    for lat in lats:
                        # if lat > 40.:
                        #     continue
                        rayname='L'+str(int(lon*100))+'L'+str(int(lat*100))
                        # print rayname
                        if rayname in dset.auxiliary_data.Raypath:
                            raylons, raylats=dset.auxiliary_data.Raypath[rayname].data.value
                            rayx, rayy=m(raylons, raylats)
                            m.plot(rayx, rayy, 'g', lw=1)
            except:
                pass
        if prop:
            self.plot_propagation(inbasemap=m, factor=factor)
        if showfig:
            plt.show()
        return
    
    
    def plot_diffa_special(self, projection='lambert', factor=3, prop=True, geopolygons=None, cmap='seismic', vmin=-16, vmax=16, showfig=True, rayASDF='', rayfactor=30):
        """Plot data with contour
        """
        m=self._get_basemap(projection=projection, geopolygons=geopolygons)
        if self.lonArr.shape[0]-2==self.diffaArr.shape[0] and self.lonArr.shape[1]-2==self.diffaArr.shape[1]:
            self.cut_edge(1,1)
        elif self.lonArr.shape[0]!=self.diffaArr.shape[0] or self.lonArr.shape[1]!=self.diffaArr.shape[1]:
            raise ValueError('Incompatible shape for deflection and lon/lat array!')
        x, y=m(self.lonArr, self.latArr)
        cmap=pycpt.load.gmtColormap('./GMT_panoply.cpt')
        cmap='bwr'
        cmap =discrete_cmap(int(vmax-vmin)/4, cmap)
        im=m.pcolormesh(x, y, self.diffaArr, cmap=cmap, shading='gouraud', vmin=vmin, vmax=vmax)
        cb = m.colorbar(im, "bottom", size="3%", pad='2%')
        cb.ax.tick_params(labelsize=10)
        cb.set_label('Great Circle Deflection (deg)', fontsize=15, rotation=0)
        cb.ax.tick_params(labelsize=15)
        #
        evx, evy=m(129., 41.306)
        plt.plot(evx, evy, 'y*', markersize=20)
        #
        if os.path.isfile(rayASDF):
            try:
                dset=pyasdf.ASDFDataSet(rayASDF)
                dlon=self.dlon*rayfactor
                dlat=self.dlat*rayfactor
                Nlon=round((self.maxlon-self.minlon)/dlon)+1
                Nlat=round((self.maxlat-self.minlat)/dlat)+1
                lons=np.arange(Nlon)*dlon+self.minlon
                lats=np.arange(Nlat)*dlat+self.minlat
                for lon in lons:
                    for lat in lats:
                        evlo=129.0; evla=41.306
                        az, baz, dist = geodist.inv(evlo, evla, lon, lat) 
                        if az>-105.5 or az < -106.5:
                            continue
                        rayname='L'+str(int(lon*100))+'L'+str(int(lat*100))
                        if rayname in dset.auxiliary_data.Raypath:
                            raylons, raylats=dset.auxiliary_data.Raypath[rayname].data.value
                            rayx, rayy=m(raylons, raylats)
                            m.plot(rayx, rayy, 'g', lw=1)
            except:
                pass
        if prop:
            self.plot_propagation(inbasemap=m, factor=factor)
        if showfig:
            plt.show()
        return
    
    def plot_propagation(self, projection='lambert', inbasemap=None, factor=3, showfig=False):
        """Plot propagation direction
        """
        if inbasemap==None:
            m=self._get_basemap(projection=projection)
        else:
            m=inbasemap
        if self.lonArr.shape[0]-2==self.grad[0].shape[0] and self.lonArr.shape[1]-2==self.grad[0].shape[1]:
            self.cut_edge(1,1)
        elif self.lonArr.shape[0]!=self.grad[0].shape[0] or self.lonArr.shape[1]!=self.grad[0].shape[1]:
            raise ValueError('Incompatible shape for gradient and lon/lat array!')
        normArr = np.sqrt ( ma.getdata(self.grad[0] )** 2 + ma.getdata(self.grad[1]) ** 2)
        x, y=m(self.lonArr, self.latArr)
        U=self.grad[1]/normArr
        V=self.grad[0]/normArr
        if factor!=None:
            x=x[0:self.Nlat:factor, 0:self.Nlon:factor]
            y=y[0:self.Nlat:factor, 0:self.Nlon:factor]
            U=U[0:self.Nlat:factor, 0:self.Nlon:factor]
            V=V[0:self.Nlat:factor, 0:self.Nlon:factor]
        Q = m.quiver(x, y, U, V, scale=50, width=0.001)
        if showfig:
            plt.show()
        return
    
    def plot_appV(self, projection='lambert', geopolygons=None, showfig=True,  vmin=None, vmax=None):
        """Plot data with contour
        """
        m=self._get_basemap(projection=projection, geopolygons=geopolygons)
        x, y=m(self.lonArr, self.latArr)
        cmap = colors.get_colormap('tomo_80_perc_linear_lightness')
        im=m.pcolormesh(x, y, self.appV, cmap=cmap, shading='gouraud', vmin=vmin, vmax=vmax)
        # ###
        # stlaLst=np.arange(5)*0.25+33.25
        # stlo=110.
        # for stla in stlaLst:
        #     xs, ys=m(stlo, stla)
        #     plt.plot(xs,ys,'^g', markersize=10)
        # stlaLst=np.arange(5)*0.25+29
        # stlo=115.
        # for stla in stlaLst:
        #     xs, ys=m(stlo, stla)
        #     plt.plot(xs,ys,'^b', markersize=15)
        ###
        ###
        
        # xs, ys=m(106.5, 33.5)
        # plt.plot(xs,ys,'^g', markersize=10)
        try:
            vrange=vmin+np.arange((vmax-vmin)/0.1+1)*0.1
            cb = m.colorbar(im, "bottom", size="3%", pad='2%', ticks=vrange)
        except:
            cb = m.colorbar(im, "bottom", size="3%", pad='2%')
        cb.ax.tick_params(labelsize=15)
        # cb.set_label(r"$\frac{\mathrm{km}}{\mathrm{s}}$", fontsize=15, rotation=0)
        cb.set_label("Apparent Phase Velocity (km/sec)", fontsize=15, rotation=0)
        if showfig:
            plt.show()
        return
    
    
    def get_distArr(self, evlo, evla):
        """Get epicentral distance array
        """
        evloArr=np.ones(self.lonArr.shape)*evlo
        evlaArr=np.ones(self.lonArr.shape)*evla
        g = Geod(ellps='WGS84')
        az, baz, distevent = g.inv(self.lonArr, self.latArr, evloArr, evlaArr)
        distevent=distevent/1000.
        self.distArr=distevent
        return
    
    
    def reset_reason(self, dist):
        self.reason_n=np.zeros(self.Zarr.shape)
        # self.reason_n[self.distArr>dist]=7
        self.reason_n[self.distArr<dist]=7
        self.reason_n[self.reason_n==5]=0
        # evlo=129.0
        # evla=41.306
        # evloArr=np.ones(self.lonArr.shape)*evlo
        # evlaArr=np.ones(self.lonArr.shape)*evla
        # g = Geod(ellps='WGS84')
        # az, baz, distevent = g.inv(self.lonArr, self.latArr, evloArr, evlaArr)
        # distevent=distevent/1000.
        # self.dd=distevent
        # self.reason_n[self.lonArr<90.]=6
        # self.reason_n[(self.lonArr<100.)*(self.latArr<30.)]=6
        # self.reason_n[(self.lonArr<100.)*(self.latArr<30.)]=6
        # self.reason_n[(self.latArr<23.)]=6
        # self.reason_n[(self.latArr<25.)*(self.lonArr>110.)]=6
        # self.reason_n[(self.latArr<30.)*(self.lonArr>120.)]=6
        # self.reason_n[(self.latArr<35.)*(self.lonArr>130.)]=6
        # self.reason_n[(self.lonArr>132.)]=6
        # self.reason_n[(self.lonArr>132.)]=6
        # self.reason_n[(self.latArr>43.)*(self.lonArr<=120.)]=6
        # self.reason_n[(self.latArr>48.)*(self.lonArr>120.)]=6
        return
    
    
    
    def compare(self, inField):
        self.compareArr=self.diffaArr-inField.diffaArr
        self.compareArr=ma.masked_array(self.compareArr, mask=np.zeros(self.reason_n.shape) )
        self.compareArr.mask[self.reason_n!=0]=1
        return
    
    def plot_compare(self, projection='lambert', geopolygons=None, cmap='seismic', vmin=-20, vmax=20):
        """Plot data with contour
        """
        m=self._get_basemap(projection=projection, geopolygons=geopolygons)
        if self.lonArr.shape[0]-2==self.compareArr.shape[0] and self.lonArr.shape[1]-2==self.compareArr.shape[1]:
            self.cut_egde(1,1)
        elif self.lonArr.shape[0]!=self.compareArr.shape[0] or self.lonArr.shape[1]!=self.compareArr.shape[1]:
            raise ValueError('Incompatible shape for deflection and lon/lat array!')
        x, y=m(self.lonArr, self.latArr)
        cmap =discrete_cmap(int(vmax-vmin)/4+1, cmap)
        m.pcolormesh(x, y, self.compareArr, cmap=cmap, shading='gouraud', vmin=vmin, vmax=vmax)
        m.colorbar()
        plt.show()
        return
    
    def histogram(self):
        import matplotlib.mlab as mlab
        indiffa=self.compareArr[self.reason_n==0]
        # the histogram of the data
        n, bins, patches = plt.hist(indiffa, 50, normed=1, facecolor='green', alpha=0.75)
        
        # # add a 'best fit' line
        # y = mlab.normpdf( bins, mu, sigma)
        # l = plt.plot(bins, y, 'r--', linewidth=1)
        
        plt.xlabel('Deflection')
        plt.ylabel('Probability')
        plt.title('ses3d')
        plt.axis([-15, 15, 0, 0.25])
        plt.grid(True)
        
        plt.show()
        
    def mean_std(self):
        inarr=self.compareArr[self.reason_n==0]
        # print np.where(inarr<-15.)[0].size + np.where(inarr>15.)[0].size
        # inarr=inarr[(inarr<15.)*(inarr>-15)]
        print 'mean',np.mean(inarr),'std:',np.std(inarr)
        # print 'mean',np.mean(self.compareArr[self.reason_n==0]),'std:',np.std(self.compareArr[self.reason_n==0])
    
    
    def plot_field_sta(self, projection='lambert', contour=True, geopolygons=None, showfig=True, vmin=None, vmax=None):
        """Plot data with contour
        """
        self.reason_n=np.zeros(self.Zarr.shape)
        self.get_distArr(evlo=129.029, evla=41.306)
        # self.reason_n[self.distArr<150.]=1
        self.Zarr[self.distArr<180.]=9999
        # self.np2ma()
        m=self._get_basemap(projection=projection, geopolygons=geopolygons)
        x, y=m(self.lonArr, self.latArr)
        im=m.pcolormesh(x, y, self.Zarr, cmap='gist_ncar_r', shading='gouraud', vmin=vmin, vmax=vmax)
        cb = m.colorbar(im, "right", size="3%", pad='2%')
        cb.ax.tick_params(labelsize=10)
        ####################################
        import stations
        sfile1='./ms_experiment/az_235.0_236.0/recfile_1'
        sfile2='./ms_experiment/az_253.0_254.0/recfile_1'
        sfile3='./ms_experiment/az_300.0_301.0/recfile_1'
        # sfile1='./ms_experiment/dist_475.0_525.0/recfile_1'
        # sfile2='./ms_experiment/dist_1075.0_1125.0/recfile_1'
        # sfile3='./ms_experiment/dist_1675.0_1725.0/recfile_1'
        SLst1=stations.StaLst()
        SLst2=stations.StaLst()
        SLst3=stations.StaLst()
        SLst1.read(sfile1)
        SLst2.read(sfile2)
        SLst3.read(sfile3)
        lons=np.array([])
        lats=np.array([])
        for sta in SLst1:
            lons=np.append(lons, sta.lon)
            lats=np.append(lats, sta.lat)
        xsta, ysta = m(lons, lats)
        m.plot(xsta, ysta, 'b^', markersize=10)
        lons=np.array([])
        lats=np.array([])
        for sta in SLst2:
            lons=np.append(lons, sta.lon)
            lats=np.append(lats, sta.lat)
        xsta, ysta = m(lons, lats)
        m.plot(xsta, ysta, 'g^', markersize=10)
        lons=np.array([])
        lats=np.array([])
        for sta in SLst3:
            lons=np.append(lons, sta.lon)
            lats=np.append(lats, sta.lat)
        xsta, ysta = m(lons, lats)
        m.plot(xsta, ysta, 'r^', markersize=10)
        ####################################
        if self.fieldtype=='Tph' or self.fieldtype=='Tgr':
            cb.set_label('sec', fontsize=3, rotation=0)
        if self.fieldtype=='Amp':
            cb.set_label('nm', fontsize=10, rotation=0)
        if contour:
            levels=np.linspace(ma.getdata(self.Zarr).min(), ma.getdata(self.Zarr).max(), 20)
            m.contour(x, y, self.Zarr, colors='k', levels=levels, linewidths=0.5)
        if showfig:
            plt.show()
        return

        
    def fftDiff(self, m, n):
        """Compute derivative with FFT
        """
        try:
            h = np.fft.fft2(self.Zarr)
        except:
            h = np.fft.fft2(self.Zarr)
        hshift = np.fft.fftshift(h)
        Nlon=self.Nlon
        Nlat=self.Nlat
        if Nlon % 2 ==0:
            u=np.arange(Nlon) - Nlon/2.
        else:
            u=np.arange(Nlon) - (Nlon-1)/2.
        if Nlat % 2 ==0:
            v=np.arange(Nlat) - Nlat/2.
        else:
            v=np.arange(Nlat) - (Nlat-1)/2.
        U,V=np.meshgrid(u,v)
        hdiff =  ((1j*2*np.pi*U/Nlon)**m)*((1j*2*np.pi*V/Nlat)**n) * hshift
        out_diff = np.real( np.fft.ifft2( np.fft.ifftshift(hdiff) ) )
        out_diff = out_diff[:self.Nlat, :self.Nlon]/(self.dlon_kmArr**m)/(self.dlat_kmArr**n)
        return out_diff
    
    def fftDiff2(self, m, n):
        """Compute derivative with FFT, extend grid number to be power of 2
        """
        Nlon=1<<(self.Nlon-1).bit_length()
        Nlat=1<<(self.Nlat-1).bit_length()
        # h = np.fft.fft2(self.Zarr, s=[Nlon, Nlat] )
        h = np.fft.fft2(ma.getdata(self.Zarr), s=[Nlon, Nlat] )
        hshift = np.fft.fftshift(h)
        u = np.arange(Nlon) - Nlon/2.
        v = np.arange(Nlat) - Nlat/2.
        U,V = np.meshgrid(u,v)
        hdiff = ( (1j*2*np.pi*U/Nlon)**m )*( (1j*2*np.pi*V/Nlat)**n ) * hshift
        out_diff = np.real( np.fft.ifft2( np.fft.ifftshift(hdiff) ) )
        out_diff = out_diff[:self.Nlat, :self.Nlon]/(self.dlon_kmArr**m)/(self.dlat_kmArr**n)
        return out_diff
    
    
            
                
                    
    

    

