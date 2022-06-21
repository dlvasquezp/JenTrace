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
from JenTrace.ray_trc import trace
from JenTrace.mrt_fnc import LMN_apertureStop,XYZ_apertureStop,XYZ_image
from JenTrace.spt_dgm import spot_diagram
from scipy.optimize import minimize, brute, fmin
import matplotlib.pyplot as plt
import numpy as np

class OpDesign:
    '''
    dsn_src: point source 
    opt_sys: optical system
    
    Falta: Documentacion
    '''
    def __init__(self,usrSrc,optSys,aprRad=1.0,aprInd=1,systemType='default'):
        #Attributes
        self.usrSrc  = usrSrc
        self.optSys  = optSys
        self.aprRad  = 1.0 
        self.aprInd  = 1
        self.designType = systemType
        #Assign aperture attributes
        self.change_aperture_radius(aprRad)
        self.change_aperture_index(aprInd)
        #design attributes
        self.dsgPtoSrc = PointSource([0,0,0],self.usrSrc.Wavelength) 
        self.dsgInfSrc = InfinitySource([0,0,1],self.usrSrc.Wavelength)
        #Solve attriutes
        self.dsgSolved = False
        self.dsgError  = []
        self.tolError  = 0.005
        
        # solve design
        self.solve_dsg()
        
    def change_aperture_radius(self,aprRad):
        assert isinstance(aprRad,(int,float)), 'Invalid aperture radius (aperRad) data type'
        assert aprRad > 0, 'Invalid aperture (aperRad) radius value' 
        self.aprRad=aprRad
        
    def change_aperture_index(self,aprInd):
        surf_len = len(self.optSys.SurfaceData)
        assert isinstance(aprInd,int),'Invalid aperture index (aprInd) data type'
        assert aprInd > 0 and aprInd < surf_len, 'Invalid aperture index (aprInd) value'
        self.aprInd=aprInd
        
    def solve_dsg(self):
        self.dsgSolved = False
        self.propagate_essential_rays()
        self.trace_optical_design()
        if max(self.dsgError)< self.tolError:
            self.dsgSolved = True
            #print(self.dsgError)
        else:
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
        
    def trace_optical_design(self):
        self.raySrcTrace = trace(self.usrSrc.RayList   ,self.optSys.SurfaceData)
        self.dsgPtoTrace = trace(self.dsgPtoSrc.RayList,self.optSys.SurfaceData)
        self.dsgInfTrace = trace(self.dsgInfSrc.RayList,self.optSys.SurfaceData)
        
    def propagate_essential_rays(self):
        self.dsgError  = []
        usrSrcError=[]
        dsgPtoSrcError=[]
        dsgInfSrcError=[]
        
        for rayIndex in range(5):
            #User ray source
            if isinstance(self.usrSrc,PointSource):
                LMN, rayError = self.propagate_ray (self.usrSrc   , rayIndex)
                self.usrSrc.change_LMN(LMN,rayIndex)
                usrSrcError.append(rayError)
            if isinstance(self.usrSrc,InfinitySource):
                XYZ, rayError = self.propagate_ray (self.usrSrc   , rayIndex)
                self.usrSrc.change_XYZ(XYZ,rayIndex)
                usrSrcError.append(rayError)
            #Design point source
            LMN, rayError = self.propagate_ray (self.dsgPtoSrc, rayIndex)
            self.dsgPtoSrc.change_LMN(LMN,rayIndex)
            dsgPtoSrcError.append(rayError)
            #Design source at infinity
            XYZ,rayError = self.propagate_ray (self.dsgInfSrc, rayIndex)
            self.dsgInfSrc.change_XYZ(XYZ,rayIndex)
            dsgInfSrcError.append(rayError)
        
        self.dsgError.append(sum(usrSrcError))
        self.dsgError.append(sum(dsgPtoSrcError))
        self.dsgError.append(sum(dsgInfSrcError))
            
    
    def propagate_ray (self,ptoSrc,rayIndex):
        
        if isinstance(ptoSrc,PointSource):
            #Perform optimization
            x0  = [0,0]
            res = minimize(LMN_apertureStop, x0, args=(self,ptoSrc,rayIndex), method='Nelder-Mead')
            # Get result
            rayError = res.fun
            x1  = res.x
            # If error is out of boundary, try brute algorithm near x1
            if (rayError > self.tolError):
                rranges = (slice(x1[0]-0.05, x1[0]+0.05, 0.01), slice(x1[1]-0.05, x1[1]+0.05, 0.01))
                resbrute = brute(LMN_apertureStop, rranges,args=(self,ptoSrc,rayIndex), full_output=True,finish=fmin)
                rayError = resbrute[1]
                x1  = resbrute[0]
            LMN = RaySource.calc_direcCos([x1[0],x1[1],1])
            return LMN, rayError
        
        if isinstance(ptoSrc,InfinitySource):
            #Initial point is calculated with image distance and inclination
            d   = self.optSys.SurfaceData[0][0]
            m   = ptoSrc.DirecCos
            #Perform optimization
            x0  = [-d*m[0],-d*m[1]]
            res = minimize(XYZ_apertureStop, x0, args=(self,ptoSrc,rayIndex), method='Nelder-Mead')
            # Get result
            rayError = res.fun
            x1  = res.x
            # If error is out of boundary, try brute algorithm near x1
            if (rayError > self.tolError):
                rranges = (slice(x1[0]-0.05, x1[0]+0.05, 0.01), slice(x1[1]-0.05, x1[1]+0.05, 0.01))
                resbrute = brute(XYZ_apertureStop, rranges,args=(self,ptoSrc,rayIndex), full_output=True,finish=fmin)
                rayError = resbrute[1]
                x1  = resbrute[0]
            XYZ = [x1[0],x1[1],0]
            return XYZ, rayError
            
    def autofocus(self):
        #Perform optimization
        #rayIndex = 1
        x0 = self.optSys.SurfaceData[-2][0]
        sptSrc,sptTrace = spot_diagram(self,noRays=1000)
        #res= minimize(XYZ_image, x0,args=(self,rayIndex),method='Nelder-Mead')
        res= minimize(XYZ_image, x0,args=(self,sptSrc),method='Nelder-Mead')
        #Replace value
        x1 = res.x
        self.optSys.SurfaceData[-2][0]=x1[0]
        #Actualize trace
        self.solve_dsg()
        
    def plot_design(self,clearSemDia=[]):
        fig, ax = plt.subplots()
        fig, ax = plot_system(self,fig,ax,clearSemDia=clearSemDia)
        fig, ax = plot_rayTrace(self.raySrcTrace,fig=fig,ax=ax, color='b')
        fig, ax = plot_rayTrace(self.dsgPtoTrace,fig=fig,ax=ax, color='g')
        ax.set_xlabel('z[mm]')
        ax.set_ylabel('y[mm]')

    
if __name__=='__main__':
    
    # Instantiate optical system
    syst1 = OpSysData()
    syst1.change_surface(20 ,0 ,1 ,surfIndex=0)
    syst1.add_surface(    1 ,+0.20,1.42)
    syst1.add_surface(    2 ,-0.1 ,1.5 )
    syst1.add_surface(    5 ,-0.2 ,1.0 )
    clearSemDia=[1,2,1.8,1.8,1.0]
    syst1.plot_optical_system(clearSemDia)
    
    # Instantiate point source
    pto1  = PointSource([0,0.9,0],635)
    pto2  = InfinitySource(RaySource.calc_direcCos([+0.0,-0.4,1.0]), 635)
    
    design1  = OpDesign(pto1,syst1,aprRad=1.0,aprInd=1)
    #design2  = OpDesign(pto2,syst1,aprRad=1.6,aprInd=2)
    design1.autofocus() 
    
    #RayTrace = design2.raySrcTrace
    #Pto      = design2.usrSrc
    
    #Plot design
    #clearSemDia=[1,2,1.8,1.8]
    design1.plot_design(clearSemDia)
    pto1.print_report()
    #design2.plot_design(clearSemDia)
    #pto2.print_report()
    
    #Documentacion
    spot_diagram(design1,show=True)

