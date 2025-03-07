## **************************************************************************
## anteom
##  This is a very hacky version of the AnitaGeomTool code from https://github.com/anitaNeutrino/eventReaderRoot/tree/master
## **************************************************************************


import numpy as np
from pyproj import Transformer
import math
#import scipy.spatial.transform.Rotation

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

    # get an array of rotation matrices from two arrays of vectors
    def get_rot_matrices(self,A, B):
        assert A.shape == B.shape
        v = np.cross(A, B)
        #print("v.shape",v.shape)
        s = np.linalg.norm(v,axis=1)
        #print("s.shape",s.shape)
        c = np.dot(A, B[0])
        #print("c.shape",c.shape)
        vx = np.zeros( (s.shape[0],3,3))
        #print("vx.shape",vx.shape)
        vx[:,0,1]=-v[:,2]
        vx[:,0,2]=v[:,1]
        vx[:,1,0]=v[:,2]
        vx[:,1,2]=-v[:,0]
        vx[:,2,0]=-v[:,1]
        vx[:,2,1]=v[:,0]
        #vx = np.array([[0, -v[2], v[1]], [v[2], 0, -v[0]], [-v[1], v[0], 0]]) 
        #vxdot=np.einsum('...j,...j',vx,vx)
        
        #print("vx.shape",vx.shape)
        vxdot=np.einsum('ijl,ilk->ijk',vx,vx)
        #print("vxdot.shape",vxdot.shape)
        #r = np.eye(3) + vx + (np.dot(vx, vx) * (1-c)/(s**2))
        r = np.eye(3) + vx + (np.einsum('ijk,i->ijk',vxdot,(1-c)/(s**2)))
        return r

    def getThetaPhiWaveWais(self,lat,lon,alt,heading):
        #First up we get the balloon posiiton in cartesian coordinatres
        x,y,z=self.latlonaltToXYZ(lat,lon,alt)
        x2,y2,z2=self.latlonaltToXYZ(lat,lon,alt+100)
    
        #Next we determine the vector of local up at the baloon position
        mag2=np.sqrt((x2-x)**2 + (y2-y)**2 + (z2-z)**2)
        upx=(x2-x)/mag2
        upy=(y2-y)/mag2
        upz=(z2-z)/mag2
        up=np.array([upx,upy,upz]).T
        #print("up.shape",up.shape)
    
        #Now we get the vector in the direction of WAIS
        waisDirx=self.waisX-x
        waisDiry=self.waisY-y
        waisDirz=self.waisZ-z
        waisDir=np.array([waisDirx,waisDiry,waisDirz]).T
        #print(waisDir.shape)
        #Now we tile the zhat unit vector so we can do some awful matrix multiplication
        B=np.tile([0,0,1],(up.shape[0],1))
        #print("B.shape",B.shape)
        #print("B[0]",B[0])
    
        #Then we get the array of rotation matrices from local up to [0,0,1]
        rotArray=self.get_rot_matrices(up,B)
    
        #Now we apply the rotation to the direction to WAIS
        waisDirPrime=np.einsum('ijk,ik->ij',rotArray,waisDir)
        #waisDirPrime=np.array([r.dot(w) for r,w in zip(rotArray,waisDir)])
        #print("waisDirPrime.shape",waisDirPrime.shape)
    
        #Last step is to determine theta and phi wave
        r=np.linalg.norm(waisDirPrime,axis=1)
        rxy=np.sqrt(waisDirPrime[:,0]**2 + waisDirPrime[:,1]**2)
        theta=np.arccos(waisDirPrime[:,2]/r)
        phi=np.sign(waisDirPrime[:,1])*np.arccos(waisDirPrime[:,0]/rxy)
    
        #The final step would be to add or subract the heading (after converting it from degrees to radians)
        phi= (phi + np.deg2rad(heading) )  #Can't remember if we add or subtract
        return theta,np.mod(phi,2*np.pi)

    
