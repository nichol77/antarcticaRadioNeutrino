## **************************************************************************
## anteom
##  This is a very hacky version of the AnitaGeomTool code from https://github.com/anitaNeutrino/eventReaderRoot/tree/master
## **************************************************************************


import numpy as np
from pyproj import Transformer

class AntGeom:

    def __init__(self):
        # This transformer takes us from WGS84 latitiude longitude altitude (https://epsg.io/4979) to WGS84 equatorial x y z (https://epsg.io/4978)
        self.trans_GPS_to_XYZ = Transformer.from_crs(4979, 4978, always_xy=True)
        self.trans_GPS_to_APS = Transformer.from_crs(4979,3031)
        
        #Anita 4 location (from Ben Strutt on Slack)
        #And now with an elog note https://www.phys.hawaii.edu/elog/anita_notes/699
        self.LATITUDE_WAIS_A4 = (-79.468116) #///< Latitude of WAIS divide pulser
        self.LONGITUDE_WAIS_A4 = -(112.059258) #///< Longitude of WAIS divide pulser
        self.ALTITUDE_WAIS_A4 = 1779.80 #///< Altitude of WAIS divide pulser

        self.waisX,self.waisY,self.waisZ=self.latlonaltToXYZ(self.LATITUDE_WAIS_A4,self.LONGITUDE_WAIS_A4,self.ALTITUDE_WAIS_A4)



        #///< Map from antenna to SURF. Both polarizations from an antenna go to the same SURF.
        self.antToSurfMap=np.array([11,5,10,4,11,4,10,5,11,5,10,4,11,4,10,5,
                    				9,3,8,2,8,3,9,2,9,3,8,2,8,3,9,2,
                    				6,0,7,1,6,1,7,0,6,0,7,1,6,1,7,0])
  
        #///< Map for VPOL channel of antenna to channel on SURF. (VPOL channels are 0-3)  
        self.vAntToChan=np.array([3,1,3,5,1,3,1,3,2,0,2,0,0,2,0,2,
                    			1,3,1,3,3,1,3,1,0,2,0,2,2,0,2,0,
                				3,1,3,1,1,3,1,3,2,0,2,0,0,2,0,2])
				  
        #///< Map for HPOL channel of antenna to channel on SURF. (HPOL channels are 4-7)
        self.hAntToChan=np.array([7,5,7,1,5,7,5,7,6,4,6,4,4,6,4,6,
                                			 5,7,5,7,7,5,7,5,4,6,4,6,6,4,6,4,
                            				 7,5,7,5,5,7,5,7,6,4,6,4,4,6,4,6])
  
        self.surfChanToAnt=np.zeros((12,9))
        for ant in range(len(self.antToSurfMap)):
            surf=self.antToSurfMap[ant]
            v=self.vAntToChan[ant]
            h=self.hAntToChan[ant]
            self.surfChanToAnt[surf][v]=-1*(ant+1)
            self.surfChanToAnt[surf][h]=+1*(ant+1)


    def latlonToAntarctica(self,lat,lon):
        return self.trans_GPS_to_APS.transform(lat,lon)
        
    def latlonaltToXYZ(self,lat,lon,alt):
        return self.trans_GPS_to_XYZ.transform(lon,lat,alt)

    def getDistToWais(self,lat,lon,alt):
        x,y,z=self.latlonaltToXYZ(lat,lon,alt)
        return np.sqrt( (x-self.waisX)**2 + (y-self.waisY)**2 + (z-self.waisZ)**2)

    
