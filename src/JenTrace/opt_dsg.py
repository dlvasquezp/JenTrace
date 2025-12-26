# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 21:21:14 2020
@author: David Vasquez
Class OpDesign
"""
try: import JenTrace
except ModuleNotFoundError: 
    import os, sys
    sys.path.insert(0,os.path.dirname(os.getcwd()))

from JenTrace.plt_fnc import plot_system, plot_rayTrace
from JenTrace.opt_sys import OpSysData
from JenTrace.ray_src import RaySource,PointSource,InfinitySource
from JenTrace.ray_trc import trace,print_report
from JenTrace.mrt_fnc import LMN_apertureStop,XYZ_apertureStop,XYZ_image
from JenTrace.spt_dgm import spot_diagram
from scipy.optimize import minimize, brute, fmin
import matplotlib.pyplot as plt
import numpy as np
import random

class OpDesign:
    '''
    dsn_src: point source / infinity source 
    opt_sys: optical system
    
    Falta: Documentacion
    '''
    def __init__(self,usrSrc,optSys,aprRad=1.0,aprInd=1,systemType='default'):
        #Attributes
        assert isinstance(usrSrc,(PointSource,InfinitySource)), 'Ray source should be either PointSource or InfinitSource'
        self.usrSrc  = usrSrc
        self.optSys  = optSys
        self.aprRad  = 1.0 
        self.aprInd  = int(1)
        self.pupRad  = 1.0
        self.pupPos  = 1.0
        self.designType = systemType
        #Assign aperture attributes
        self.change_aperture_radius(aprRad)
        self.change_aperture_index(aprInd)
        #design attributes
        self.dsgPtoSrc = PointSource([0,0,0],self.usrSrc.Wavelength) 
        self.dsgInfSrc = InfinitySource([0,0,1],self.usrSrc.Wavelength)
        #solve attributes
        self.dsgSolved = False
        self.dsgError  = []
        self.tolError  = 0.005
        
        #Initial Ray estimation
        self.initRayTrace = self.initial_ray_estimation()
        #solve design
        self.solve_dsg()
        
    def change_aperture_radius(self,aprRad):
        assert isinstance(aprRad,(int,float)), 'Invalid aperture radius (aperRad) data type'
        assert aprRad > 0, 'Invalid aperture (aperRad) radius value' 
        self.aprRad=aprRad
        
    def change_aperture_index(self,aprInd):
        surf_len = len(self.optSys.SurfaceData)
        assert isinstance(aprInd,int),'Invalid aperture index (aprInd) data type'
        assert aprInd > 0 and aprInd < (surf_len-1), 'Invalid aperture index (aprInd) value'
        self.aprInd=aprInd
    
    def initial_ray_estimation(self):
        samSrcWvln= self.usrSrc.Wavelength
        
        #Define directions from the point source to the first surface rand
        surf1Position = self.optSys.SurfaceData[0][0]
        surf1SemiDia  = self.optSys.SurfaceData[0][4]
        
        noRings=9
        centers=[]
        spacing= 1/(noRings+2)
        for q in range(-noRings, noRings + 1):
            r1 = max(-noRings, -q - noRings)
            r2 = min(noRings, -q + noRings)
            for r in range(r1, r2 + 1):
                # axial to cartesian (pointy-top orientation)
                x = spacing * (np.sqrt(3) * q + np.sqrt(3)/2 * r)
                y = spacing * (3/2 * r)
                centers.append([x, y])  
        centers=np.array(centers)
        
        mask = np.sqrt(centers[:, 0]**2 + centers[:, 1]**2) <= 1
        pts=centers[mask]
        pts = np.multiply(pts,surf1SemiDia)
        
        raysPupXYZ = []
        for px, py in pts:
            raysPupXYZ.append(np.array([px,py,surf1Position]))
        
        if isinstance(self.usrSrc,PointSource):
            XYZ= self.usrSrc.Position 

            raysImg2Pup = []
            for ray in raysPupXYZ:
                raysImg2Pup.append(ray-XYZ)
            
            LMN     = RaySource.calc_direcCos(raysImg2Pup[0])
            samSrc  = RaySource (XYZ,LMN,samSrcWvln)
            
            for ray in raysImg2Pup [1:]:
                LMN     = RaySource.calc_direcCos(ray)
                samSrc.new_ray(XYZ,LMN,samSrcWvln)
            
        if isinstance(self.usrSrc,InfinitySource): 
            LMN= self.usrSrc.DirecCos 

            raysImg2Pup = []
            for ray in raysPupXYZ:
                raysImg2Pup.append(list(ray-(np.multiply(LMN,surf1Position/LMN[2]))))
            
            samSrc  = RaySource (raysImg2Pup[0],LMN,samSrcWvln)
            
            for XYZ in raysImg2Pup [1:]:
                samSrc.new_ray(XYZ,LMN,samSrcWvln)
            
        # Make ray trace
        initRayTrace = trace(samSrc.RayList, self.optSys.SurfaceData)

        return initRayTrace 

        
    def solve_dsg(self, caller=None):
        self.dsgSolved = False
        #self.initRayTrace = self.initial_ray_estimation()
        self.propagate_essential_rays(caller)
        self.trace_optical_design()
        if max(self.dsgError)< self.tolError:
            self.dsgSolved = True
            self.calculate_entrance_pupil()
        else:
            print(self.dsgError)
            raise Warning('numerical error out of bounds, check aperture index or radius')
            '''
            if self.dsgError[2] > self.tolError:
                if self.designType =='telecentric':
                    self.dsgSolved   = True
                    self.dsgInfSrc   = np.nan
                    self.dsgInfTrace = np.nan
                    self.dsgError[2] = np.nan
                else:
                    print(self.dsgError)
                    raise Warning('numerical error out of bounds, to disable this warning due to telecentricity, pass designType="telecentric"')
            else:
                print(self.dsgError)
                raise Warning('numerical error out of bounds, check aperture index or radius')
            '''
    def trace_optical_design(self):
        self.raySrcTrace = trace(self.usrSrc.RayList   ,self.optSys.SurfaceData)
        self.dsgPtoTrace = trace(self.dsgPtoSrc.RayList,self.optSys.SurfaceData)
        #self.dsgInfTrace = trace(self.dsgInfSrc.RayList,self.optSys.SurfaceData)
        
    def propagate_essential_rays(self, caller=None):
        self.dsgError  = []
        usrSrcError=[]
        dsgPtoSrcError=[]
        #dsgInfSrcError=[]
        for rayIndex in range(5):
            #User ray source
            if isinstance(self.usrSrc,PointSource):
                # Calculate initial optimization direction
                if caller == 'minimization':
                    cosDirZ = self.dsgPtoTrace[rayIndex,21,0]
                    x0  = [self.dsgPtoTrace[rayIndex,19,0]/cosDirZ,self.dsgPtoTrace[rayIndex,20,0]/cosDirZ]
                else:
                    nearIdx = self.nearest_aperture_ray(rayIndex)
                    cosDirZ = self.initRayTrace[nearIdx,21,0]
                    x0  = [self.initRayTrace[nearIdx,19,0]/cosDirZ,self.initRayTrace[nearIdx,20,0]/cosDirZ]
                LMN, rayError = self.propagate_ray (self.usrSrc   , rayIndex, x0)
                self.usrSrc.change_LMN(LMN,rayIndex)
                usrSrcError.append(rayError)
            if isinstance(self.usrSrc,InfinitySource):
                # Calculate intial optimization position
                d   = self.optSys.SurfaceData[0][0]
                m   = self.usrSrc.DirecCos
                x0  = [-d*m[0],-d*m[1]]
                XYZ, rayError = self.propagate_ray (self.usrSrc   , rayIndex, x0)
                self.usrSrc.change_XYZ(XYZ,rayIndex)
                usrSrcError.append(rayError)
            #Design point source
            LMN, rayError = self.propagate_ray (self.dsgPtoSrc, rayIndex, [0,0])
            self.dsgPtoSrc.change_LMN(LMN,rayIndex)
            dsgPtoSrcError.append(rayError)
            #Design source at infinity
            #XYZ,rayError = self.propagate_ray (self.dsgInfSrc, rayIndex, [0,0])
            #self.dsgInfSrc.change_XYZ(XYZ,rayIndex)
            #dsgInfSrcError.append(rayError)
        
        self.dsgError.append(sum(usrSrcError))
        self.dsgError.append(sum(dsgPtoSrcError))
        #self.dsgError.append(sum(dsgInfSrcError))
        
    def nearest_aperture_ray(self,rayIndex)->int:
        xPosList = self.initRayTrace[:,8,self.aprInd]
        yPosList = self.initRayTrace[:,9,self.aprInd]
        
        if rayIndex == 0:
            nearIdx = np.argmin((np.power(xPosList,2)+np.power(yPosList,2)))
        
        if rayIndex == 1:
            distDif = np.array(yPosList)-self.aprRad
            nearIdx = np.argmin((np.power(xPosList,2)+np.power(distDif,2)))

        if rayIndex == 2:
            distDif = np.array(yPosList)+self.aprRad
            nearIdx = np.argmin((np.power(xPosList,2)+np.power(distDif,2)))
            
        if rayIndex == 3:
            distDif = np.array(xPosList)-self.aprRad
            nearIdx = np.argmin((np.power(distDif,2)+np.power(yPosList,2)))
            
        if rayIndex == 4:
            distDif = np.array(xPosList)+self.aprRad
            nearIdx = np.argmin((np.power(distDif,2)+np.power(yPosList,2)))
        
        return nearIdx
            
    
    def propagate_ray (self,ptoSrc,rayIndex,x0):
        
        if isinstance(ptoSrc,PointSource):
            res = minimize(LMN_apertureStop, x0, args=(self,ptoSrc,rayIndex), method='Nelder-Mead')
            # Get result
            rayError = res.fun
            x1  = res.x
            # If error is out of boundary, try method='Newton-CG' x1
            if (rayError > self.tolError/3):
                print('Powell PS')
                res = minimize(LMN_apertureStop, x0, args=(self,ptoSrc,rayIndex), method='Powell')
                # Get result
                rayError = res.fun
                x1  = res.x
            # If error is out of boundary, try brute algorithm near x1
            if (rayError > self.tolError/3):
                print('brute PS, prev. error {}'.format(rayError))
                rranges = (slice(x1[0]-0.05, x1[0]+0.05, 0.01), slice(x1[1]-0.05, x1[1]+0.05, 0.01))
                resbrute = brute(LMN_apertureStop, rranges,args=(self,ptoSrc,rayIndex), full_output=True,finish=fmin)
                rayError = resbrute[1]
                x1  = resbrute[0]
                # Get result
                #rayError = res.fun
                #x1  = res.x
                print(rayError)
            LMN = RaySource.calc_direcCos([x1[0],x1[1],1])
            return LMN, rayError
        
        if isinstance(ptoSrc,InfinitySource):
            res = minimize(XYZ_apertureStop, x0, args=(self,ptoSrc,rayIndex), method='Nelder-Mead')
            # Get result
            rayError = res.fun
            x1  = res.x
            # If error is out of boundary, try brute algorithm near x1
            if (rayError > self.tolError/3):
                print('brute IS, prev. error {}'.format(rayError))
                rranges = (slice(x1[0]-0.05, x1[0]+0.05, 0.01), slice(x1[1]-0.05, x1[1]+0.05, 0.01))
                resbrute = brute(XYZ_apertureStop, rranges,args=(self,ptoSrc,rayIndex), full_output=True,finish=fmin)
                rayError = resbrute[1]
                x1  = resbrute[0]
                # Get result 
                #rayError = res.fun
                print(rayError)
            XYZ = [x1[0],x1[1],0]
            return XYZ, rayError
            
    def autofocus(self):
        #Perform optimization
        x0 = self.optSys.SurfaceData[-2][0]
        sptSrc,sptTrace = spot_diagram(self,noRays=100)
        res= minimize(XYZ_image, x0,args=(self,sptSrc),method='Nelder-Mead')
        #Replace value
        x1 = res.x
        self.optSys.SurfaceData[-2][0]=x1[0]
        #Actualize trace
        self.solve_dsg()
        
    def plot(self):
        print(self.optSys)
        fig, ax = plt.subplots()
        fig, ax = plot_system(self,fig,ax)
        fig, ax = plot_rayTrace(self.raySrcTrace,fig=fig,ax=ax, color='b')
        fig, ax = plot_rayTrace(self.dsgPtoTrace,fig=fig,ax=ax, color='g')
        #fig, ax = plot_rayTrace(self.dsgInfTrace,fig=fig,ax=ax, color='r')
        ax.set_xlabel('z[mm]')
        ax.set_ylabel('y[mm]')
        ax.axis('equal')
        
    def calculate_entrance_pupil(self):
        if self.dsgSolved:
            #Find pupil Z position
            chiefRay = self.usrSrc.RayList[0]
            chiefRayXYZ = np.array(chiefRay[0])
            chiefRayLMN = np.array(chiefRay[1])
            #Find chiefRay optical axis intersection
            fun = lambda z: np.sum(np.power(chiefRayXYZ+np.multiply(chiefRayLMN, z/chiefRayLMN[2]),2)[0:2])
            x0 = [self.optSys.SurfaceData[0][0]]
            res = minimize(fun, x0, method='Nelder-Mead')
            pupilZ = res.x
            self.pupPos = pupilZ[0]
            
            #Find mean pupil radius
            pupilSize=[]
            for q in range(1,5):
                marginalRay= self.usrSrc.RayList[q]
                marginalRayXYZ = np.array(marginalRay[0])
                marginalRayLMN = np.array(marginalRay[1])
                #Propagate marginal rays to the entrance pupil
                pupilRand = marginalRayXYZ + np.multiply(marginalRayLMN,pupilZ/marginalRayLMN[2])
                pupilSize.append(pupilRand[0:2])
            pupilRadius = np.mean([pupilSize[0][1],-pupilSize[1][1],pupilSize[2][0],-pupilSize[3][0]])
            self.pupRad = pupilRadius
                                

    
if __name__=='__main__':
    
    # Instantiate optical system
    syst1 = OpSysData()
    syst1.change_surface(20 ,0 ,1 ,surfIndex=0)
    syst1.add_surface(    1 ,+0.20,1.42)
    syst1.add_surface(    2 ,-0.1 ,1.5 )
    syst1.add_surface(    5 ,-0.2 ,1.0 )
    clearSemDia=[1,2,1.8,1.8,1.0]
    syst1.change_clearSemDia(clearSemDia)
    syst1.plot()
    
    # Point source test
    pto1  = PointSource([0,3,0],635)
    design1  = OpDesign(pto1,syst1,aprRad=1.0,aprInd=3)
    design1.plot()
    design1.autofocus() 
    design1.plot()
    
    #design1.calculate_ray_source()
    '''
    pto2  = InfinitySource(RaySource.calc_direcCos([+0.0,-0.4,1.0]), 635)
    
    #design2  = OpDesign(pto2,syst1,aprRad=1.6,aprInd=2)
    
    
    #RayTrace = design2.raySrcTrace
    #Pto      = design2.usrSrc
    
    #Plot design
    #clearSemDia=[1,2,1.8,1.8]
    design1.plot()
    pto1.print_report()
    #design2.plot_design(clearSemDia)
    pto2.print_report()
    
    #Documentacion
    spot_diagram(design1,show=True)
    '''

