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

def spot_diagram(optDsg, noRings=9, noRays=1000, show=False, plotType ='posXYZ', surfIndex=-1, color='b'):
    '''
    spot_diagram creates a bunch of rays from the source and porpagate it throw the system. 
    The rays are filtereed (and discarted) using the ray position in the aperture plane.
    
    This function return a RayTrace (sptTrace) where the rays information can be extracted.
    
    optDsg: optical Design object
    noRays: Initial number of rays to porpagate. Note that due to the filtering, the final rays are less
    show: Boolean to show the plot
    plotType: Two types are posible, either the ray position (posXYZ) or ray cosine direction (cosDir) 
    surfIndex: Surface index of the selected plane. By dafault is the image plane selected
    color: plot color, see Matplotlib.
    '''
    
    if optDsg.dsgSolved == True:
        optSys = optDsg.optSys
        
        # Generate ray source with ray pattern
        samSrc=ray_pattern(optDsg,'hexapolar', noRings=noRings)
             
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
    
def ray_pattern(optDsg, pattern:str, noRings:int, noRays=100, show=False, color='b'):
    pupilRad = optDsg.pupRad
    pupilPos = optDsg.pupPos
    usrSrcWvln= optDsg.usrSrc.Wavelength
    centers = []
    
    if pattern == 'hexapolar':
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
    
    if pattern == 'sagital':
        x = np.linspace(-1,1,(noRings*2+1))
        y = np.zeros(len(x))
        centers=[[q,w] for q,w in zip (x, y)]  
        centers=np.array(centers)
        
    if pattern == 'tangential':
        y = np.linspace(-1,1,(noRings*2+1))
        x = np.zeros(len(y))
        centers=[[q,w] for q,w in zip (x, y)]  
        centers=np.array(centers)
        
    mask = np.sqrt(centers[:, 0]**2 + centers[:, 1]**2) <= 1
    pts=centers[mask]
    if show:
        fig, ax = plt.subplots(figsize=(6, 6))
        ax.scatter(pts[:, 0], pts[:, 1], s=40)
        circle = plt.Circle((0, 0), 1, fill=False)
        ax.add_patch(circle)
    
    #Resize pupil radious from 1 to design value
    pts = np.multiply(pts,pupilRad)
    
    raysPupXYZ = []
    for px, py in pts:
        raysPupXYZ.append(np.array([px,py,pupilPos]))
    
    if isinstance(optDsg.usrSrc,PointSource):
        XYZ= optDsg.usrSrc.Position 

        raysImg2Pup = []
        for ray in raysPupXYZ:
            raysImg2Pup.append(ray-XYZ)
    
        LMN     = RaySource.calc_direcCos(raysImg2Pup[0])
        samSrc  = RaySource (XYZ,LMN,usrSrcWvln)
    
        for ray in raysImg2Pup [1:]:
            LMN     = RaySource.calc_direcCos(ray)
            samSrc.new_ray(XYZ,LMN,usrSrcWvln)
    
    if isinstance(optDsg.usrSrc,InfinitySource): 
        LMN= optDsg.usrSrc.DirecCos 

        raysImg2Pup = []
        for ray in raysPupXYZ:
            raysImg2Pup.append(list(ray-(np.multiply(LMN,pupilPos/LMN[2]))))
        
        samSrc  = RaySource (raysImg2Pup[0],LMN,usrSrcWvln)
        
        for XYZ in raysImg2Pup [1:]:
            samSrc.new_ray(XYZ,LMN,usrSrcWvln)
    
    return samSrc
        
if __name__ == '__main__':
    import time
    from opt_sys import OpSysData
    from opt_dsg import OpDesign
    
    start = time.time()
    # Instantiate optical system
    syst1 = OpSysData()
    syst1.change_surface(30     ,0         ,1      ,surfIndex=0) 
    syst1.add_surface   (3.50   ,1/12.37   ,'N-BK7')
    syst1.add_surface   (1.50   ,1/-11.10  ,'N-SF5')
    syst1.add_surface   (5     ,1/-25.47  ,1.4     )
    syst1.add_surface   (5     ,0         ,1      )
    #syst1.changeAperture(1,surfIndex=2)
    clearSemDia=[1,5.0,5.0,5.0,5.0,1]
    syst1.plot(clearSemDia)
    
    # Instantiate point source
    pto1  = PointSource([4,4,0],635)
    
    # Instantiate optical design
    design1  = OpDesign(pto1,syst1,aprRad=2,aprInd=4)
    
    # autofocus
    design1.autofocus()
    
    # plot design
    design1.plot()
    samSrc1,samTrace1 = spot_diagram(design1,show=True)
    
    # Instantiate infinity source
    pto2    = InfinitySource(RaySource.calc_direcCos([+0.0,-0.4,1.0]), 635)
    design2 = OpDesign(pto2, syst1, aprInd=4)
    design2.autofocus()
    design2.plot()
    samSrc2,samTrace2 = spot_diagram(design2,show=True)
    
    print_report(samTrace2)
    
    #print_report(samTrace2, 'prop',index=8)
    #print_report(samTrace2, 'prop',index=9)
    t1=ray_pattern(design1,'hexapolar', noRings=5)
    t2=ray_pattern(design1,'sagital', noRings=9)
    t3=ray_pattern(design2,'tangential', noRings=9)
    end = time.time()
    print('-------------------->>>>>>>>>>>>>',end - start)
    
