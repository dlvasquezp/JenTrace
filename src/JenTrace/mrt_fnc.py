# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 06:47:46 2020
@author: David Vasquez
functions:  -LMN_apertureStop
            -XYZ_apertureStop
            -ap_error_calc
            -XYZ_image
"""
try: import JenTrace
except ModuleNotFoundError: 
    import os, sys
    sys.path.insert(0,os.path.dirname(os.getcwd()))
    
from JenTrace.ray_trc import trace
import numpy as np

def LMN_apertureStop (x0,*arg):
    '''
    Merit function to calculate the error produced by the rays propaged in the 
    direction [x0[0],x0[1],1] and the aperture stop
    
    x0:   (list) [vector_compX, vector_comp_y] 
    arg: (list) [object Optical desing, object ray source, int indexRay]
    
    arg[0] must be Optical design
    arg[1] must be point source (pointSource)
    arg[2] must be int [0,1,2,3,4]
    '''   
    assert arg[0].__class__.__name__=='OpDesign' ,'arg[0] must be an OpDesign object' 
    assert arg[1].__class__.__name__=='PointSource' ,'arg[1] must be a PointSource object' 
    assert arg[2] >=0 and arg[2] <=4, 'Invalid ray index (indexRay)'
    #reference the values
    vector        = x0 
    SurfaceData   = arg[0].optSys.SurfaceData
    ApertureRadio = arg[0].aprRad
    ApertureIndex = arg[0].aprInd
    RayList       = arg[1].RayList
    indexRay      = arg[2]
    
    #Calculate direction cosines
    LMN = arg[1].calc_direcCos([vector[0],vector[1],1])  
    #replace cosine director    
    arg[1].change_LMN(LMN,indexRay)
    #Make Raytrace
    RayTrace  = trace(RayList,SurfaceData)
    #Calculate error
    error = ap_error_calc(RayTrace,ApertureRadio,indexRay,ApertureIndex)
    
    return error

def XYZ_apertureStop (x0,*arg):
    '''
    Merit function to calculate the error produced by the rays propaged in the 
    position [x0[0],x0[1],0] and the aperture stop
    
    x0:   (list) [position X, position y] 
    arg: (list) [object Optical desing, object ray source, int indexRay]
    
    arg[0] must be Optical design
    arg[1] must be source at infinity (infinitySource)
    arg[2] must be int [0,1,2,3,4]
    '''
    assert arg[0].__class__.__name__=='OpDesign' ,'arg[0] must be an OpDesign object' 
    assert arg[1].__class__.__name__=='InfinitySource' ,'arg[1] must be an InfinitySource object' 
    assert arg[2] >=0 and arg[2] <=4, 'Invalid ray index (indexRay)'
    
    #reference the values
    vector        = x0 
    SurfaceData   = arg[0].optSys.SurfaceData
    ApertureRadio = arg[0].aprRad
    ApertureIndex = arg[0].aprInd
    RayList       = arg[1].RayList
    indexRay      = arg[2]    
    
    #Rewrite position vector
    XYZ = [vector[0],vector[1],0]   
    #replace cosine director    
    arg[1].change_XYZ(XYZ,indexRay)
    #Make Raytrace
    RayTrace  = trace(RayList,SurfaceData)
    #Calculate error
    error = ap_error_calc(RayTrace,ApertureRadio,indexRay,ApertureIndex)
    
    return error

def ap_error_calc(RayTrace,ApertureRadio,indexRay,ApertureIndex):
    """
    if indexRay == 0:
        error    = (abs(RayTrace[indexRay, 9 ,ApertureIndex])
                   +abs(RayTrace[indexRay, 8 ,ApertureIndex]))
        
    if indexRay == 1:
        error    = (abs(+ApertureRadio - RayTrace[indexRay, 9 ,ApertureIndex])
                   +abs(RayTrace[indexRay, 8 ,ApertureIndex]))
                    
    if indexRay == 2:
        error    = (abs(-ApertureRadio - RayTrace[indexRay, 9 ,ApertureIndex])
                   +abs(RayTrace[indexRay, 8 ,ApertureIndex]))
                   
    if indexRay == 3:
        error    = (abs(+ApertureRadio - RayTrace[indexRay, 8 ,ApertureIndex])
                   +abs(RayTrace[indexRay, 9 ,ApertureIndex]))
        
    if indexRay == 4:
        error    = (abs(-ApertureRadio - RayTrace[indexRay, 8 ,ApertureIndex])
                   +abs(RayTrace[indexRay, 9 ,ApertureIndex]))
    """
    if indexRay == 0:
        error    = (abs(RayTrace[indexRay, 4 ,ApertureIndex])
                   +abs(RayTrace[indexRay, 3 ,ApertureIndex]))
        
    if indexRay == 1:
        error    = (abs(+ApertureRadio - RayTrace[indexRay, 4 ,ApertureIndex])
                   +abs(RayTrace[indexRay, 3 ,ApertureIndex]))
                    
    if indexRay == 2:
        error    = (abs(-ApertureRadio - RayTrace[indexRay, 4 ,ApertureIndex])
                   +abs(RayTrace[indexRay, 3 ,ApertureIndex]))
                   
    if indexRay == 3:
        error    = (abs(+ApertureRadio - RayTrace[indexRay, 3 ,ApertureIndex])
                   +abs(RayTrace[indexRay, 4 ,ApertureIndex]))
        
    if indexRay == 4:
        error    = (abs(-ApertureRadio - RayTrace[indexRay, 3 ,ApertureIndex])
                   +abs(RayTrace[indexRay, 4 ,ApertureIndex]))
    
    return error
    

def XYZ_image (x0,*arg):
    '''
    Merit function to calculate the distance from a ray and the last 
    surface (image) origin [0,0,0]. If the point source is in the optical axis, 
    the rays with index from 1 to 4 describe the focus error distance.
    
    x0:   (list[float]) distance
    *arg: (list) [object optical design, int indexRay]
    
    # x0 must be float
    # arg[0] must be Optical design
    # arg[1] must be int [0,1,2,3,4]
    
    '''
    #rename values
    dist        = x0[0] 
    #RayList     = arg[0].dsgPtoSrc.RayList
    #RayList     = arg[0].usrSrc.RayList
    SurfaceData = arg[0].optSys.SurfaceData
    #indexRay    = arg[1]
    sptSrc      = arg[1]
    
    #replace surface distance
    surf_len = len(SurfaceData)
    surf_idx = (surf_len-2) 
    arg[0].optSys.change_surface(dist,SurfaceData[surf_idx][1],SurfaceData[surf_idx][2],surfIndex=surf_idx)
    
    #Make Raytrace
    #RayTrace  = trace(RayList,SurfaceData)
    RayTrace  = trace(sptSrc.RayList,SurfaceData)
    
    #Fist moment of inertia
    xCoor     = RayTrace[:,8,-1]
    yCoor     = RayTrace[:,9,-1]
    N         = len(xCoor)
    centroidX = sum(xCoor) / N 
    centroidY = sum(yCoor) / N
    
    #Root mean square radius
    ms =np.sum(np.power(xCoor-centroidX,2)+np.power(yCoor-centroidY,2))
    mrs=np.sqrt(ms)
    
    error = mrs
    
    '''
    #Calculate error
    error = (abs(RayTrace[indexRay, 7 ,-1])
            +abs(RayTrace[indexRay, 8 ,-1])
            +abs(RayTrace[indexRay, 9 ,-1]))
    
    error = (abs(RayTrace[0, 7 ,-1] - RayTrace[indexRay, 7 ,-1])
            +abs(RayTrace[0, 8 ,-1] - RayTrace[indexRay, 8 ,-1])
            +abs(RayTrace[0, 9 ,-1] - RayTrace[indexRay, 9 ,-1]))
    '''
    return error



if __name__=='__main__':
    from ray_src import PointSource
    from opt_sys import OpSysData
    from opt_dsg import OpDesign
    from plt_fnc import plot_system,plot_rayTrace 

    import matplotlib.pyplot as plt
    
    pto1  = PointSource([0,1,0],635)
    syst1 = OpSysData()
    syst1.add_surface(2,0.05,1.7)
    syst1.add_surface(10,-0.5,1.4)
    #syst1.changeAperture(1,surfIndex = 1)
    
    design1  = OpDesign(pto1,syst1)
    design1.autofocus()
    
    fig, ax = plt.subplots()
    fig, ax = plot_system(design1, fig=fig, ax=ax)
    fig, ax = plot_rayTrace(design1.raySrcTrace,fig=fig,ax=ax)
    fig, ax = plot_rayTrace(design1.dsgPtoTrace,fig=fig,ax=ax)
    fig, ax = plot_rayTrace(design1.dsgInfTrace,fig=fig,ax=ax)
    
    
    