#!/usr/bin/env python

import numpy as np
import obspy
import warnings
from lasif import colors
import os
import matplotlib.pyplot as plt
from mpl_toolkits.basemap import Basemap, shiftgrid

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
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
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
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
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
                # print lines
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
    An object to handle input model for FMST
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
            if inV[i]==0:
                continue
            lon=inlon[i]
            lat=inlat[i]
            index = np.where((self.lonArr==lon)*(self.latArr==lat))
            self.vArr[index[0], index[1]]=inV[i]
        self._change_cushion()
            # print lon, lat, self.lonArr[index[0], index[1]], self.latArr[index[0], index[1]]
        return
    
    def write4field2d(self, outfname):
        OutArr=np.append(self.lonArr, self.latArr)
        OutArr=np.append(OutArr, self.vArr)
        OutArr=OutArr.reshape(3, self.lonArr.size)
        OutArr=OutArr.T
        np.savetxt(outfname, OutArr, fmt='%g')
        return
    
    def smooth(self, sigma):
        v_filtered=self.cushion_vArr.copy()
        for iteration in xrange(int(sigma)):
            for i in np.arange(1,self.lat.size+1):
                for j in np.arange(1,self.lon.size+1):
                    v_filtered[i,j]=(self.cushion_vArr[i,j]+self.cushion_vArr[i+1,j]+self.cushion_vArr[i-1,j]+self.cushion_vArr[i,j+1]+self.cushion_vArr[i,j-1])/5.0
        self.cushion_vArr=v_filtered
        self._change_v()
        return
    
    def plot(self, projection='lambert', geopolygons=None):
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
                                self.maxlat+2, self.minlon) # distance is in m
            m = Basemap(width=distEW, height=distNS, rsphere=(6378137.00,6356752.3142), resolution='l', projection='lcc',\
                lat_1=self.minlat, lat_2=self.maxlat, lon_0=lon_centre, lat_0=lat_centre+1.5)
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
        m.pcolormesh(x, y, self.vArr, cmap=cmap, shading='gouraud', vmin=2.9, vmax=3.4)
        try:
            geopolygons.PlotPolygon(inbasemap=m)
        except:
            pass
        m.colorbar()
        plt.show()
        return
    
    def write(self, outdir):
        """
        """
        if not os.path.isdir(outdir):
            os.makedirs(outdir)
        outfname=outdir+'/gridi.vtx'
        cushion_vArr=self.cushion_vArr[::-1, :]
        cushion_errV=self.cushion_errV[::-1, :]
        # outvArr=self.vArr[::-1, :]
        with open(outfname, 'w') as f:
            f.writelines('%12.d%12.d\n' %(self.lat.size, self.lon.size))
            f.writelines('%14.8f%14.8f\n' %(self.lat.max(), self.lon.min()))
            f.writelines('%14.8f%14.8f\n\n' %(self.dlat, self.dlon))
            for ilon in xrange(self.Nlon+2):
                for ilat in xrange(self.Nlat+2):
                    f.writelines('%12.8f%12.8f\n' %(cushion_vArr[ilat, ilon], cushion_errV[ilat, ilon]))
                f.writelines('\n')
        return

