import numpy as np
import pyasdf
from pyproj import Geod


class rayASDF(pyasdf.ASDFDataSet):

    def read_raydat(self, rayfile, period=10., factor=100):
        raylons=np.array([])
        raylats=np.array([])
        with open(rayfile, 'r') as f:
            inlines=f.readlines()
            for line in inlines:
                if line.split()[0]=='>':
                    if raylons.size!=0:
                        data=np.append(raylons, raylats)
                        data=data.reshape(2, raylons.size)
                        parameters={'lon':raylons[0], 'lat':raylats[0]}
                        path='L'+str(int(raylons[0]*100))+'L'+str(int(raylats[0]*100))
                        self.add_auxiliary_data(data=data,
                            data_type='Raypath', path=path, parameters=parameters)
                    raylons=np.array([])
                    raylats=np.array([])
                    # print lonArr.size
                    continue
                lArr=line.split()
                raylons=np.append(raylons, float(lArr[0]))
                raylats=np.append(raylats, float(lArr[1]))
            
    
    def get_pre_amp(self, field2d, outfname, verbose=True):
        
        lonArr=np.array([])
        latArr=np.array([])
        ampArr=np.array([])
        # lonArr=np.append(lonArr, raylons[0])
        # latArr=np.append(latArr, raylats[0])
        for dname in self.auxiliary_data.Raypath.list():
            subdset=self.auxiliary_data.Raypath[dname]
            lon=subdset.parameters['lon']
            lat=subdset.parameters['lat']
            data=subdset.data.value
            raylons=data[0,:]
            raylats=data[1,:]
            if verbose:
                print 'Get pre amp for:',lon,lat
            Amp=self._get_pre_amp(field2d, raylons, raylats)
            # return
            lonArr=np.append(lonArr, lon)
            latArr=np.append(latArr, lat)
            ampArr=np.append(ampArr, Amp)
        OutArr=np.append(lonArr, latArr)
        OutArr=np.append(OutArr, ampArr)
        OutArr=OutArr.reshape(3, lonArr.size)
        OutArr=OutArr.T
        np.savetxt(outfname, OutArr)
            
    
    def _get_pre_amp(self, field2d, raylons, raylats):
        dlon=field2d.dlon
        dlat=field2d.dlat
        midlons=(raylons[1:]+raylons[:-1])/2
        midlats=(raylats[1:]+raylats[:-1])/2
        lons0=raylons[:-1]
        lonsf=raylons[1:]
        lats0=raylats[:-1]
        latsf=raylats[1:]
        g = Geod(ellps='WGS84')
        az, baz, dsArr = g.inv(lons0, lats0, lonsf, latsf)
        dsArr=dsArr/1000.
        numbr=midlons.size
        lplcArr=np.array([])
        appV=np.array([])
        for i in xrange(numbr):
            lon=midlons[i]
            lat=midlats[i]
            # ilon=int(lon/dlon)
            # ilat=int(lat/dlat)
            ilon=(abs(field2d.lon-lon)).argmin()
            ilat=(abs(field2d.lat-lat)).argmin()
            # clon=ilon*dlon
            # clat=ilat*dlat
            # index=np.where((field2d.lonArr==clon)*(field2d.latArr==clat))
            lplcArr=np.append(lplcArr, field2d.lplc[ilat, ilon])
            appV=np.append(appV, field2d.appV[ilat, ilon])
            # print index, lon, clon, lat, clat, (abs(field2d.lon-lon)).argmin(), (abs(field2d.lat-lat)).argmin()
            # if index[0].size==0:
            #     lonArr=np.ones(field2d.lonArr.shape)*lon
            #     latArr=np.ones(field2d.lonArr.shape)*lat
            #     az, baz, distArr = g.inv(field2d.lonArr, field2d.latArr, lonArr, latArr)
            #     print distArr.argmin()
                
            # dlonArr=abs(field2d.lonArr-lon)
            # dlatArr=abs(field2d.latArr-lat)
            # index=np.where((dlonArr<dlon)*(dlatArr<dlat))
        corr=np.sum(lplcArr*dsArr*appV)
        Amp=np.exp(-0.5*corr)
        # print appV
        return Amp
            
            
        
        
        
            
    