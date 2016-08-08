#!/usr/bin/env python

import numpy as np
import obspy
import warnings
from lasif import colors
import os

class StaInfo(object):
    """
    An object contains a station information
    ========================================================================
    General Parameters:
    stacode    - station name
    network    - network
    lon,lat    - position for station
    ========================================================================
    """
    def __init__(self, stacode=None, lat=None, lon=None, network='LF'):

        self.stacode=stacode
        self.network=network
        self.lon=lon
        self.lat=lat

    def get_contents(self):
        if self.stacode==None:
            print 'StaInfo NOT Initialized yet!'
            return
        print 'Network:%16s' %(self.network)
        print 'Station:%20s' %(self.stacode)
        print 'Longtitude:%17.3f' %(self.lon)
        print 'Latitude:  %17.3f' %(self.lat)
        return
    
    def GetPoslon(self):
        if self.lon<0:
            self.lon=self.lon+360.
        return
    
    def GetNeglon(self):
        if self.lon>180.:
            self.lon=self.lon-360.
        return
        
class StaLst(object):
    """
    An object contains a station list(a list of StaInfo object) information several methods for station list related analysis.
        stations: list of StaInfo
    """
    def __init__(self,stations=None):
        self.stations=[]
        if isinstance(stations, StaInfo):
            stations = [stations]
        if stations:
            self.stations.extend(stations)

    def __add__(self, other):
        """
        Add two StaLst with self += other.
        """
        if isinstance(other, StaInfo):
            other = StaLst([other])
        if not isinstance(other, StaLst):
            raise TypeError
        stations = self.stations + other.stations
        return self.__class__(stations=stations)

    def __len__(self):
        """
        Return the number of StaInfos in the StaLst object.
        """
        return len(self.stations)

    def __getitem__(self, index):
        """
        __getitem__ method of StaLst objects.
        :return: StaInfo objects
        """
        if isinstance(index, slice):
            return self.__class__(stations=self.stations.__getitem__(index))
        else:
            return self.stations.__getitem__(index)

    def append(self, station):
        """
        Append a single StaInfo object to the current StaLst object.
        """
        if isinstance(station, StaInfo):
            self.stations.append(station)
        else:
            msg = 'Append only supports a single StaInfo object as an argument.'
            raise TypeError(msg)
        return self

    def read(self, datadir, network='LF', issource=False):
        """
        Read Sation List from a txt file
        stacode longitude latidute network
        """
        if issource:
            infname=datadir+'/sources.dat'
        else:
            infname=datadir+'/receivers.dat'
        with open(infname, 'r') as f:
            Sta=[]
            f.readline()
            L=0
            for lines in f.readlines():
                L+=1
                lines=lines.split()
                lat=float(lines[0])
                lon=float(lines[1])
                stacode='S%06d' %L
                netsta=network+'.'+stacode
                if Sta.__contains__(netsta):
                    index=Sta.index(netsta)
                    if abs(self[index].lon-lon) >0.01 and abs(self[index].lat-lat) >0.01:
                        raise ValueError('Incompatible Station Location:' + netsta+' in Station List!')
                    else:
                        print 'Warning: Repeated Station:' +netsta+' in Station List!'
                        continue
                Sta.append(netsta)
                self.append(StaInfo (stacode=stacode, network=network, lon=lon, lat=lat ))
        return

    def HomoStaLst(self, minlat, Nlat, minlon, Nlon, dlat, dlon, net='LF'):
        """
        Generate equal grid interval station list
        ========================================================
        Input Parameters:
        minlat, minlon  - minimum latitude/longitude
        Nlat, Nlon      - number of latitude/longitude grid
        dlat, dlon      - latitude/longitude interval
        net             - network
        prx             - prefix for station name
        ========================================================
        """
        L=1
        for ilon in xrange(Nlon):
            for ilat in xrange(Nlat):
                lon=minlon+ilon*dlon
                lat=minlat+ilat*dlat
                stacode='S%06d' %L
                self.stations.append(StaInfo (stacode=stacode, network=net, lon=lon, lat=lat))
                if ilon == Nlon -1 and ilat == Nlat -1:
                    print 'maxlat=', lat, 'maxlon=',lon
                L+=1
        return
    
    def write(self, outdir, issource=False):
        """Write station list to output directory
        """
        L=len(self.stations)
        if issource:
            outfname=outdir+'/sources.dat'
        else:
            outfname=outdir+'/receivers.dat'
        with open(outfname,'wb') as f:
            f.writelines('%g\n' %L )
            for station in self.stations:
                f.writelines('%5.1f %5.1f\n' %(station.lat, station.lon) )
        return
    
    def write_otime(self, outdir, ELst, err=0.1, v0=3.0):
        
        with open(outdir+'/otimes.dat','wb') as f:
            for event in ELst:
                for sta in self.stations:
                    dist, az, baz=obspy.geodetics.gps2dist_azimuth(event.lat, event.lon, sta.lat, sta.lon)
                    T=dist/v0/1000.
                    f.writelines('1 %9.5f%7.4f\n' %(T, err) )
                    
        return
    
    def travel2field2d(self, datadir, outfname):
        lons=np.array([])
        lats=np.array([])
        Tarr=np.array([])
        with open(datadir+'/rtravel.out','r') as f:
            for sta in self.stations:
                lines=f.readline()
                lines=lines.split()
                T=float(lines[1])
                lons=np.append(lons, sta.lon)
                lats=np.append(lats, sta.lat)
                Tarr=np.append(Tarr, T)
        OutArr=np.append(lons, lats)
        OutArr=np.append(OutArr, Tarr)
        OutArr=OutArr.reshape(3, Tarr.size)
        OutArr=OutArr.T
        np.savetxt(outfname, OutArr, fmt='%g')
        return
                
                
                
    

class vmodel(object):
    """
    An object to handle input model for SPECFEM2D.
    ===========================================================================
    Parameters:
    minlat, maxlat, minlon, maxlon  - bound of study region
    dlat, dlon                      - latitude/longitude interval                                    
    v0                              - background velocity
    ===========================================================================
    """
    def __init__(self, minlat, maxlat, minlon, maxlon, dlat, dlon, v0=3.0, errV=0.3):
        self.Nlat=int((maxlat-minlat)/dlat)+1
        self.dlat=dlat
        self.Nlon=int((maxlon-minlon)/dlon)+1
        self.dlon=dlon
        self.lon=np.arange(self.Nlon)*self.dlon+minlon
        self.lat=np.arange(self.Nlat)*self.dlat+minlat
        self.minlat=minlat
        self.maxlat=self.lat.max()
        self.maxlon=self.lon.max()
        self.minlon=minlon
        if maxlon!=self.maxlon or maxlat!=self.maxlat:
            warnings.warn('maxlat='+ str(self.maxlat)+' maxlon='+str(self.maxlon), UserWarning, stacklevel=1)
        self.lonArr, self.latArr = np.meshgrid(self.lon, self.lat)
        self.vArr=np.ones(self.lonArr.shape)*v0
        self.errV=np.ones(self.lonArr.shape)*errV
        # cushion boundary is applied
        self.cushion_vArr=np.ones((self.Nlat+2, self.Nlon+2))*v0
        self.cushion_errV=np.ones((self.Nlat+2, self.Nlon+2))*errV
        return
    
    def _change_cushion(self):
        self.cushion_vArr[1:-1, 1:-1]=self.vArr
        return
    
    def _change_v(self):
        self.vArr=self.cushion_vArr[1:-1, 1:-1]
        return
    
    def read_cv(self, infname):
        """
        Read txt velocity model file
        """
        InArr=np.loadtxt(infname)
        inlon=InArr[:,0]
        inlat=InArr[:,1]
        inV=InArr[:,2]
        for i in xrange(inlon.size):
            lon=inlon[i]
            lat=inlat[i]
            index = np.where((self.lonArr==lon)*(self.latArr==lat))
            self.vArr[index[0], index[1]]=inV[i]
        self._change_cushion()
            # print lon, lat, self.lonArr[index[0], index[1]], self.latArr[index[0], index[1]]
        return
    
    def smooth(self, sigma):
        v_filtered=self.cushion_vArr.copy()
        for iteration in xrange(int(sigma)):
            for i in np.arange(1,self.lat.size+1):
                for j in np.arange(1,self.lon.size+1):
                    v_filtered[i,j]=(self.vArr[i,j]+self.vArr[i+1,j]+self.vArr[i-1,j]+self.vArr[i,j+1]+self.vArr[i,j-1])/5.0
        self.cushion_vArr=v_filtered
        self._change_v()
        return
    
    def plot(self, projection='lambert'):
        fig=plt.figure(num=None, figsize=(12, 12), dpi=80, facecolor='w', edgecolor='k')
        lat_centre = (self.maxlat+self.minlat)/2.0
        lon_centre = (self.maxlon+self.minlon)/2.0
        if projection=='merc':
            m=Basemap(projection='merc', llcrnrlat=self.minlat, urcrnrlat=self.minlat, llcrnrlon=self.minlon,
                      urcrnrlon=self.minlon, lat_ts=20, resolution=resolution)
            m.drawparallels(np.arange(self.lat_min,self.lat_max,self.d_lon), labels=[1,0,0,1])
            m.drawmeridians(np.arange(self.lon_min,self.lon_max,self.d_lat), labels=[1,0,0,1])
        
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
                                self.maxlat+1.7, self.minlon) # distance is in m
            m = Basemap(width=distEW, height=distNS,
            rsphere=(6378137.00,6356752.3142),\
            resolution='l', projection='lcc',\
            lat_1=self.minlat, lat_2=self.maxlat, lon_0=lon_centre+1.2, lat_0=lat_centre,)
            m.drawparallels(np.arange(-80.0,80.0,10.0), linewidth=2, dashes=[2,2], labels=[1,0,0,0], fontsize=15)
            m.drawmeridians(np.arange(-170.0,170.0,10.0), linewidth=2, dashes=[2,2], labels=[0,0,1,1], fontsize=15)
        m.drawcoastlines(linewidth=1.0)
        m.drawcountries()
        # m.drawmapboundary(fill_color=[1.0,1.0,1.0])
        m.fillcontinents(lake_color='#99ffff',zorder=0.2)
        # m.drawlsmask(land_color='0.8', ocean_color='#99ffff')
        m.drawmapboundary(fill_color="white")
        cmap = colors.get_colormap('tomo_80_perc_linear_lightness')
        x, y=m(self.lonArr, self.latArr)
        m.pcolormesh(x, y, self.vArr, cmap=cmap, shading='gouraud')
        m.colorbar()
        plt.show()
        return
    
    def write(self, outdir):
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        outfname=outdir+'/gridi.vtx'
        # outvArr=self.vArr[::-1, :]
        with open(outfname, 'w') as f:
            f.writelines('%12.d%12.d\n' %(self.lat.size, self.lon.size))
            f.writelines('%14.8f%14.8f\n' %(self.lat.min(), self.lon.min()))
            f.writelines('%14.8f%14.8f\n\n' %(-self.dlat, self.dlon))
            for ilon in xrange(self.Nlon+2):
                for ilat in xrange(self.Nlat+2):
                    f.writelines('%12.8f%12.8f\n' %(self.cushion_vArr[ilat, ilon], self.cushion_errV[ilat, ilon]))
                f.writelines('\n')
        return
#     
#     def BlockHomoAnomaly(self, Xmin, Xmax, Zmin, Zmax, va, dv=None):
#         """
#         Inplement block anomaly in the model for Vs
#         =============================================================================
#         Input Parameters:
#         Xmin, Xmax, Zmin, Zmax    - defines the bound
#         va                                           - anomalous velocity
#         dv                                           - velocity anomaly in percentage( default is None, which means use va )
#         =============================================================================
#         """
#         Xindex=(self.XArr>=Xmin)*(self.XArr<=Xmax)
#         Zindex=(self.ZArr>=Zmin)*(self.ZArr<=Zmax)
#         Index=Xindex*Zindex
#         if dv!=None:
#             self.VsArr[Index]=self.VsArr[Index]*(1+dv)
#         else:
#             self.VsArr[Index]=va
#         if self.plotflag==True:
#             Xindex=(self.XArrPlot>=Xmin)*(self.XArrPlot<=Xmax)
#             Zindex=(self.ZArrPlot>=Zmin)*(self.ZArrPlot<=Zmax)
#             Index=Xindex*Zindex
#             if dv!=None:
#                 self.VsArrPlot[Index]=self.VsArrPlot[Index]*(1+dv)
#             else:
#                 self.VsArrPlot[Index]=va
#         return
#     
#     def CircleHomoAnomaly(self, Xc, Zc, R, va, dv=None):
#         """
#         Inplement circle anomaly in the model for Vs
#         =============================================================================
#         Input Parameters:
#         Xc, Zc     - center of the circle
#         R             - radius
#         va           - anomalous velocity
#         dv           - velocity anomaly in percentage( default is None, which means use va )
#         =============================================================================
#         """
#         print 'Adding homo circle anomaly Xc=', Xc,' Zc=', Zc, ' R=',R
#         dArr = np.sqrt( (self.XArr-Xc)**2 + (self.ZArr-Zc)**2)
#         Index = dArr <= R
#         if dv!=None:
#             self.VsArr[Index]=self.VsArr[Index]*(1+dv)
#         else:
#             self.VsArr[Index]=va
#         if self.plotflag==True:
#             dArr = np.sqrt( (self.XArrPlot-Xc)**2 + (self.ZArrPlot-Zc)**2)
#             Index = dArr <= R
#             if dv!=None:
#                 self.VsArrPlot[Index] = self.VsArrPlot[Index]*(1+dv)
#             else:
#                 self.VsArrPlot[Index]=va
#         return
#     
#     def CircleLinearAnomaly(self, Xc, Zc, R, va, dv=None):
#         """
#         Inplement circle anomaly with linear change towards center in the model for Vs
#         Assuming the background Vs is homogeneous
#         =============================================================================
#         Input Parameters:
#         Xc, Zc     - center of the circle
#         R             - radius
#         va           - anomalous velocity
#         dv           - velocity anomaly in percentage( default is None, which means use va )
#         =============================================================================
#         """
#         dArr = np.sqrt( (self.XArr-Xc)**2 + (self.ZArr-Zc)**2)
#         if dv==None:
#             dva = va - self.Vs
#         else:
#             dva = self.Vs*dv
#         delD = R - dArr
#         IndexIn = (delD >=0)
#         # self.VsArr = 0.5 * (np.sign(delD) + 1) * delD/R * dva + self.VsArr
#         self.VsArr = IndexIn * delD/R * dva + self.VsArr
#         if self.plotflag==True:
#             dArr = np.sqrt( (self.XArrPlot-Xc)**2 + (self.ZArrPlot-Zc)**2)
#             delD = R - dArr
#             IndexIn = delD >=0
#             # self.VsArrPlot = 0.5 * (np.sign(delD) + 1) * delD/R * dva + self.VsArrPlot
#             self.VsArrPlot = IndexIn * delD/R * dva + self.VsArrPlot
#         return
#     
#     def CircleCosineAnomaly(self, Xc, Zc, R, va=None, dv=None):
#         """
#         Inplement circle anomaly with cosine change towards center in the model for Vs
#         Assuming the background Vs is homogeneous
#         =============================================================================
#         Input Parameters:
#         Xc, Zc     - center of the circle
#         R             - radius
#         va           - anomalous velocity
#         dv           - velocity anomaly in percentage( default is None, which means use va )
#         =============================================================================
#         """
#         dArr = np.sqrt( (self.XArr-Xc)**2 + (self.ZArr-Zc)**2)
#         if va !=None:
#             dva = va - self.Vs
#         else:
#             dva = self.Vs*dv
#         delD = R - dArr
#         IndexIn = (delD >=0)
#         # self.VsArr = 0.5 * (np.sign(delD) + 1) * ( 1+np.cos( np.pi* dArr / R ) ) * dva + self.VsArr
#         self.VsArr = IndexIn * ( 1+np.cos( np.pi* dArr / R ) )/2. * dva + self.VsArr
#         if self.plotflag==True:
#             dArr = np.sqrt( (self.XArrPlot-Xc)**2 + (self.ZArrPlot-Zc)**2)
#             delD = R - dArr
#             IndexIn = (delD >=0)
#             # self.VsArrPlot = 0.5 * (np.sign(delD) + 1) * ( 1+np.cos( np.pi* dArr / R ) ) * dva + self.VsArrPlot
#             self.VsArrPlot = IndexIn * ( 1+np.cos( np.pi* dArr / R ) )/2. * dva + self.VsArrPlot
#         return
#     
#     def RingHomoAnomaly(self, Xc, Zc, Rmax, Rmin, va, dv=None):
#         """
#         Inplement ring anomaly in the model for Vs
#         =============================================================================
#         Input Parameters:
#         Xc, Zc           - center of the circle
#         Rmax, Rmin - radius max/min
#         va                 - anomalous velocity
#         dv                 - velocity anomaly in percentage( default is None, which means use va )
#         =============================================================================
#         """
#         if Rmin < self.dx:
#             self.CircleHomoAnomaly(Xc=Xc, Zc=Zc, R=Rmax, va=va, dv=dv)
#             return
#         print 'Adding homo ring anomaly Xc=', Xc,' Zc=', Zc, ' Rmax=',Rmax, ' Rmin=', Rmin
#         dArr = np.sqrt( (self.XArr-Xc)**2 + (self.ZArr-Zc)**2)
#         Index = (dArr <= Rmax) * (dArr > Rmin) 
#         if dv!=None:
#             self.VsArr[Index]=self.VsArr[Index]*(1+dv)
#         else:
#             self.VsArr[Index]=va
#         if self.plotflag==True:
#             dArr = np.sqrt( (self.XArrPlot-Xc)**2 + (self.ZArrPlot-Zc)**2)
#             Index = (dArr <= Rmax) * (dArr > Rmin) 
#             if dv!=None:
#                 self.VsArrPlot[Index] = self.VsArrPlot[Index]*(1+dv)
#             else:
#                 self.VsArrPlot[Index]=va
#         return
#     
#     
#     def ASDFmodel(self, infname, per=10., phgr=1, verbose=True):
#         """
#         Read ASDF model
#         =============================================================================
#         Input Parameters:
#         infname        - input file name
#         per                - period
#         phgr              - use phase(1) or group(2) velocity
#         =============================================================================
#         """
#         dbase = pyasdf.ASDFDataSet(infname)
#         if verbose==True:
#             print dbase.auxiliary_data.Disp
#         perArr = dbase.auxiliary_data.Disp.VP000.data.value[0, :]
#         vArr=dbase.auxiliary_data.Disp.VP000.data.value[phgr, :]
#         if not np.any(perArr==per):
#             raise ValueError('Period ', per,' sec not in the theoretical dispersion curve !' )
#         Vs0 = vArr[perArr==per]* 1000.
#         self.VsArr[:] = Vs0 
#         if self.plotflag == True:
#             self.VsArrPlot[:]=Vs0 
#         for auxid in dbase.auxiliary_data.Disp.list()[1:]:
#             perArr = dbase.auxiliary_data.Disp[auxid].data.value[0, :]
#             vArr=dbase.auxiliary_data.Disp[auxid].data.value[phgr, :]
#             Rmax=dbase.auxiliary_data.Disp[auxid].parameters['Rmax']
#             Rmin=dbase.auxiliary_data.Disp[auxid].parameters['Rmin']
#             x=dbase.auxiliary_data.Disp[auxid].parameters['x']
#             y=dbase.auxiliary_data.Disp[auxid].parameters['y']
#             Vs = vArr[perArr==per] * 1000.
#             self.RingHomoAnomaly(Xc=x, Zc=y, Rmax=Rmax, Rmin=Rmin, va=Vs)
#         return
#             
#     
#     
#     def get_GLL(self):
#         """
#         Set Gauss-Lobatto-Legendre(GLL) points for a given Lagrange polynomial degree.
#         To construct a polynomial of degree n passing through n+1 data points. 
#         """
#         if self.lpd == 2:
#             knots = np.array([-1.0, 0.0, 1.0])
#         elif self.lpd == 3:
#             knots = np.array([-1.0, -0.4472135954999579, 0.4472135954999579, 1.0])
#         elif self.lpd == 4:
#             knots = np.array([-1.0, -0.6546536707079772, 0.0, 0.6546536707079772, 1.0])
#             # knots = np.array([-1.0, -1./7.*np.sqrt(21.), 0.0, 1./7.*np.sqrt(21.), 1.0])
#         elif self.lpd == 5:
#             knots = np.array([-1.0, -0.7650553239294647, -0.2852315164806451, 0.2852315164806451, 0.7650553239294647, 1.0])
#         elif self.lpd == 6:
#             knots = np.array([-1.0, -0.8302238962785670, -0.4688487934707142, 0.0, 0.4688487934707142, 0.8302238962785670, 1.0])
#         elif self.lpd == 7:
#             knots = np.array([-1.0, -0.8717401485096066, -0.5917001814331423,\
#                 -0.2092992179024789, 0.2092992179024789, 0.5917001814331423, 0.8717401485096066, 1.0])
#         self.knots=knots
#         return
#     
#     def plot(self, ds=1000, unit='km', vmin=2.5, vmax=3.5):
#         """Plot velocity model
#         =============================================================================
#         Input Parameters:
#         ds                          - grid spacing
#         unit                       - unit
#         vmin, vmax          - vmin,vmax for colorbar
#         =============================================================================
#         """
#         
#         if self.plotflag==False:
#             raise ValueError('No plot array!')
#         plt.figure(figsize=(16,13))
#         if self.regular==True:
#             plt.pcolormesh(self.XArrPlot/ds, self.ZArrPlot/ds, self.VsArrPlot/ds, cmap='seismic_r', vmin=vmin, vmax=vmax)
#         else:
#             xi = np.linspace(self.xmin, self.xmax, self.Nx*10)
#             zi = np.linspace(self.zmin, self.zmax, self.Nz*10)
#             self.xi, self.zi = np.meshgrid(xi, zi)
#             #-- Interpolating at the points in xi, yi
#             self.vi = griddata(self.XArr, self.ZArr, self.VsArr, self.xi, self.zi, 'linear')
#             plt.pcolormesh(self.xi/ds, self.zi/ds, ma.getdata(self.vi)/ds, cmap='seismic_r', vmin=vmin, vmax=vmax)
#         ##########################################
#         # plt.plot( 320, 320 , 'y*', markersize=30)
#         ##########################################
#         plt.xlabel('x('+unit+')', fontsize=30)
#         plt.ylabel('z('+unit+')', fontsize=30)
#         plt.colorbar()
#         plt.axis([self.xmin/ds, self.xmax/ds, self.zmin/ds, self.zmax/ds])
#         # plt.axis('scaled')
#         plt.yticks(fontsize=20)
#         plt.xticks(fontsize=20)
#         plt.show()
#         return
#     
#     def GetMinMaxV(self):
#         """
#         Get minimum/maximum vs 
#         """
#         vmin=self.VsArr.min()
#         vmax=self.VsArr.max()
#         return vmin, vmax
#     
# class InputChecker(object):
#     """
#     An object to check stability condition given input parameters.
#     =============================================================================
#     Parameters:
#     dt                   - time step
#     dx, dz            - element spacing 
#     fc                    - central frequency
#     lpd                 - Lagrange polynomial degree
#     vmin, vmax   - minimum/maximum velocity
#     =============================================================================
#     """
#     def __init__(self, dt, dx, dz, fc, lpd, vmin, vmax):
#         self.dt=dt
#         self.dx=dx
#         self.dz=dz
#         self.fc=fc
#         self.lpd=lpd
#         self.get_GLL()
#         self.vmin=vmin
#         self.vmax=vmax
#         return
#     
#     def get_GLL(self):
#         """
#         Set Gauss-Lobatto-Legendre(GLL) points for a given Lagrange polynomial degree.
#         To construct a polynomial of degree n passing through n+1 data points. 
#         """
#         if self.lpd == 2:
#             knots = np.array([-1.0, 0.0, 1.0])
#         elif self.lpd == 3:
#             knots = np.array([-1.0, -0.4472135954999579, 0.4472135954999579, 1.0])
#         elif self.lpd == 4:
#             knots = np.array([-1.0, -0.6546536707079772, 0.0, 0.6546536707079772, 1.0])
#         elif self.lpd == 5:
#             knots = np.array([-1.0, -0.7650553239294647, -0.2852315164806451, 0.2852315164806451, 0.7650553239294647, 1.0])
#         elif self.lpd == 6:
#             knots = np.array([-1.0, -0.8302238962785670, -0.4688487934707142, 0.0, 0.4688487934707142, 0.8302238962785670, 1.0])
#         elif self.lpd == 7:
#             knots = np.array([-1.0, -0.8717401485096066, -0.5917001814331423,\
#                 -0.2092992179024789, 0.2092992179024789, 0.5917001814331423, 0.8717401485096066, 1.0])
#         self.knots=knots
#         return
#     
#     def CheckMinLambda(self, freqfactor=2.5):
#         """
#         Check grid spacing with wavelength minimum wavelength.
#         ==============================================
#         Input Parameters:
#         freqfactor       - fmax = freqfactor*fc
#         ==============================================
#         """
#         lambdamin=self.vmin/self.fc/freqfactor
#         dxArr=self.dx*np.diff( (0.5+0.5*(self.knots)) )
#         dzArr=self.dz*np.diff( (0.5+0.5*(self.knots)) )
#         dsmax=max(dxArr.max(), dzArr.max())
#         # dsmax=max(self.dx, self.dz)
#         # Need checking! (in manual: threshold value is around 4.5 points
#         # per wavelength in elastic media and 5.5 in acoustic media), 4.5 grid points OR 4.5 element points
#         if dsmax * 4.5 > lambdamin:
#             raise ValueError('Grid spacing is too large: ', str(dsmax),' for ',lambdamin, ' m')
#         else:
#             print 'Grid spacing:', str(dsmax),'m for',lambdamin, 'm'
#         return
#     
#     def CheckCFLCondition(self, C=0.35):
#         """
#         Check Courant-Frieddrichs-Lewy stability condition
#         ===========================================
#         Input Parameters:
#         C - Courant number (default = 0.35, normally 0.3~0.4)
#         ===========================================
#         """
#         dxArr=self.dx*np.diff( (0.5+0.5*(self.knots)) )
#         dzArr=self.dz*np.diff( (0.5+0.5*(self.knots)) )
#         dsmin=min(dxArr.min(), dzArr.min())
#         dtCFL=C*dsmin/self.vmax
#         if self.dt > dtCFL:
#             raise ValueError('Time step violates Courant-Frieddrichs-Lewy Condition: ', dt, dtCFL)
#         else:
#             print 'Time Step: ',self.dt,' s Required Time Step: ', dtCFL, 's'
#         return
#     
#     def Check(self, freqfactor=2.5, C=0.35):
#         """
#         Check minimum wavelenght and Courant conditions
#         """
#         print '=========== Checking stability conditions ==========='
#         self.CheckMinLambda(freqfactor=freqfactor)
#         self.CheckCFLCondition(C=C)
#         print '===========  Stability conditions checked  ==========='
#         return 
