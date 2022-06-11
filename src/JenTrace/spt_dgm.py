# -*- coding: utf-8 -*-
"""
Created on Sat May  2 17:23:54 2020
@author: David Vasquez
Function: spot_diagram
"""
try: import JenTrace
except ModuleNotFoundError: 
    import os, sys
    sys.path.insert(0,os.path.dirname(os.getcwd()))
    
import numpy as np
import random
import matplotlib.pyplot as plt
from JenTrace.ray_trc import trace,print_report
from JenTrace.ray_src import RaySource,PointSource,InfinitySource

def spot_diagram(optDsg, noRays=1000, show=False, plotType ='posXYZ', surfIndex=-1, color='b'):
    '''
    spot_diagram creates a bunch of rays from the source and porpagate it throw the system. 
    The rays are filtereed (and discarted) using the ray position in the aperture plane.
    
    This function return a RayTrace (sptTrace) where the rays information can be extracted.
    
    optDsg: optical Design object
    noRays: Initial number of rays to porpagate. Note that due to the filtering, the final rays are less
    show: Boolean to show the plot
    plotType: Two types are posible, either the ray position (posXYZ) or ray cosine direction (cosDir) 
    surfIndex: Surface index of the selected plane. By dafault is the image plane selcted
    color: plot color, see Matplotlib.
    '''
    
    if optDsg.dsgSolved == True:
        optSys = optDsg.optSys
        random.seed(102629122021)
        if isinstance(optDsg.usrSrc,PointSource):
            usrSrc = optDsg.usrSrc
            xRad = np.array([usrSrc.RayList[rayIndex][1][0]for rayIndex in range(1,5)])
            yRad = np.array([usrSrc.RayList[rayIndex][1][1]for rayIndex in range(1,5)])
            xlim = [xRad.min(),xRad.max()]
            ylim = [yRad.min(),yRad.max()]
            
            # Instantiate sampling (sam) source
            samSrcXYZ = usrSrc.Position
            samSrcWvln= usrSrc.Wavelength
            samSrcLM  =[[random.uniform(xlim[0],xlim[1]),random.uniform(ylim[0],ylim[1])] for q in range(noRays)]
            
            LMN     = RaySource.calc_direcCos([samSrcLM[0][0],samSrcLM[0][1],1])
            samSrc  = RaySource (samSrcXYZ,LMN,samSrcWvln)
            
            for LM in samSrcLM [1:]:
                LMN     = RaySource.calc_direcCos([LM[0],LM[1],1])
                samSrc.new_ray(samSrcXYZ,LMN,samSrcWvln)
                
        if isinstance(optDsg.usrSrc,InfinitySource): 
            usrSrc = optDsg.usrSrc
            xRad = np.array([usrSrc.RayList[rayIndex][0][0]for rayIndex in range(1,5)])
            yRad = np.array([usrSrc.RayList[rayIndex][0][1]for rayIndex in range(1,5)])
            xlim = [xRad.min(),xRad.max()]
            ylim = [yRad.min(),yRad.max()]
            
            # Instantiate sampling (sam) source
            samSrcLMN = usrSrc.DirecCos
            samSrcWvln= usrSrc.Wavelength
            samSrcXYZ =[[random.uniform(xlim[0],xlim[1]),random.uniform(ylim[0],ylim[1]),0] for q in range(noRays)]
            
            samSrc  = RaySource (samSrcXYZ[0],samSrcLMN,samSrcWvln)
            
            for XYZ in samSrcXYZ [1:]:
                samSrc.new_ray(XYZ,samSrcLMN,samSrcWvln)
            
        # Make  Trace
        samTrace = trace(samSrc.RayList,optSys.SurfaceData)
        
        # Filter the rays that hit outside the aperture
        aprRad = optDsg.aprRad
        aprInd = optDsg.aprInd
        delInd = []
        for rayInd in range(len(samSrc.RayList)):
            # Ray position at aperture plane
            x = samTrace[rayInd][8][aprInd]
            y = samTrace[rayInd][9][aprInd]
            if x**2 + y**2 > aprRad**2:
                delInd.append(rayInd)
        
        sptTrace = np.delete(samTrace,delInd,0)
        for index in sorted(delInd, reverse=True):
            del samSrc.RayList[index]

        if show == True:
            if plotType =='posXYZ':
                # Chiefray is used as the reference point
                xCent = optDsg.raySrcTrace[0][8][surfIndex]
                yCent = optDsg.raySrcTrace[0][9][surfIndex]
                # XY positions are stored in x/yPtos
                xPtos = np.array([sptTrace[q][8][surfIndex] for q in range(len(sptTrace))])
                yPtos = np.array([sptTrace[q][9][surfIndex] for q in range(len(sptTrace))])
                # Plot limits
                xDelt = abs(xPtos.max()-xPtos.min())
                yDelt = abs(yPtos.max()-yPtos.min())
                pltlim= np.max([xDelt,yDelt])
                
                plt.figure()
                plt.plot(xPtos,yPtos,'o',markersize=1,color=color)
                plt.plot(xCent,yCent,'kx')
                plt.xlim(xCent-pltlim,xCent+pltlim)
                plt.ylim(yCent-pltlim,yCent+pltlim)
                plt.axis('equal')
                plt.grid('on')
                plt.xlabel('x[mm]')
                plt.ylabel('y[mm]')
                plt.show()
            
            if plotType =='cosDir':
                # Chiefray is used as the reference point
                xCent = optDsg.raySrcTrace[0][19][surfIndex]
                yCent = optDsg.raySrcTrace[0][20][surfIndex]
                # XY positions are stored in x/yPtos
                xPtos = np.array([sptTrace[q][19][surfIndex] for q in range(len(sptTrace))])
                yPtos = np.array([sptTrace[q][20][surfIndex] for q in range(len(sptTrace))])
                # Plot limits
                xDelt = abs(xPtos.max()-xPtos.min())
                yDelt = abs(yPtos.max()-yPtos.min())
                pltlim= np.max([xDelt,yDelt])
                
                plt.figure()
                plt.plot(xPtos,yPtos,'o',markersize=1,color=color)
                plt.plot(xCent,yCent,'kx')
                plt.xlim(xCent-pltlim,xCent+pltlim)
                plt.ylim(yCent-pltlim,yCent+pltlim)
                plt.axis('equal')
                plt.grid('on')
                plt.xlabel('cosDir_x[rad]')
                plt.ylabel('cosDir_y[rad]')
                plt.show()
        
        
        return samSrc,sptTrace
        
if __name__ == '__main__':
    from opt_sys import OpSysData
    from opt_dsg import OpDesign
    
    # Instantiate optical system
    syst1 = OpSysData()
    syst1.change_surface(60     ,0         ,1      ,surfIndex=0) 
    syst1.add_surface   (3.50   ,1/15.37   ,'N-BK7')
    syst1.add_surface   (1.50   ,1/-11.10  ,'N-SF5')
    syst1.add_surface   (10     ,1/-31.47  ,1      )
    #syst1.changeAperture(1,surfIndex=2)
    clearSemDia=[1,5.0,5.0,5.0,1]
    syst1.plot_optical_system(clearSemDia)
    
    # Instantiate point source
    pto1  = PointSource([0,2.5,0],635)
    
    # Instantiate optical design
    design1  = OpDesign(pto1,syst1,aprRad=2)
    
    # autofocus
    design1.autofocus()
    
    # plot design
    design1.plot_design(clearSemDia)
    samSrc1,samTrace1 = spot_diagram(design1,show=True)
    
    # Instantiate infinity source
    pto2    = InfinitySource(RaySource.calc_direcCos([+0.0,-0.4,1.0]), 635)
    design2 = OpDesign(pto2, syst1)
    design2.autofocus()
    design2.plot_design(clearSemDia)
    samSrc2,samTrace2 = spot_diagram(design2,show=True)
    
    print_report(samTrace2)
    
    #print_report(samTrace2, 'prop',index=8)
    #print_report(samTrace2, 'prop',index=9)