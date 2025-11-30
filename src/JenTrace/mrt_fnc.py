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
    

def XYZ_image (x0:float,*arg):
    '''
    Merit function to calculate the distance from a ray and the last 
    surface (image) origin [0,0,0]. If the point source is in the optical axis, 
    the rays with index from 1 to 4 describe the focus error distance.
    
    x0:   (list[float]) distance
    *arg: (list) [object optical design, object point source]
    
    # x0 must be float
    # arg[0] must be Optical design
    # arg[1] must be point source
    
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
    #arg[0].optSys.change_surface(dist,SurfaceData[surf_idx][1],SurfaceData[surf_idx][2],surfIndex=surf_idx)
    arg[0].optSys.SurfaceData[surf_idx][0]=dist
    
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
    
    return error

def score_function (x0:list,*arg) -> float:
    '''
    Merit wizard to generalize the construction of merit functions
    x0:   (list[floats]) variable values
    *arg: (list) [0] optical design object,
                 [1] list of ray parameters,
                 [2] list of error functions,
                 [3] list of function parameters,
                 [4] List of weights
    
    '''
    #Rename arguments
    opticalDesign  = arg[0]
    raysParam      = arg[1]
    errorFun       = arg[2]
    funParam       = arg[3]
    weights        = arg[4]
    
    #Replace x0 in optical design
    opticalDesign.optSys.set_varValues(x0)
    opticalDesign.solve_dsg()
    
    #Generate RaySources
    raySourceList=[]
    for ray in raysParam:
        raySourceList.append(opticalDesign.calcRaySource(ray))
        
    #Make Raytrace
    rayTraceList = []
    for source in raySourceList:
        rayTraceList.append(trace(source.RayList,opticalDesign.SurfaceData))
        
    #Evaluate error functions
    errorList = []
    for fun, rayTrace, param , w in zip(errorFun, rayTraceList, funParam,weights):
        errorList.append(fun(rayTrace)*w)
    
    error = np.sum(errorList)
    
    return error

def spotSqrt (RayTrace,*arg)->[float]:
    '''
    Merit function to calculate the mean sqrt of the Raytrace in the specified surface.
    
    RayTrace:   (list[float]) distance
    *arg: (list) [surf_idx]
    
    # surf_idx:[int] -> Surface index where the mean sqrt value of the spot is calculated
    
    '''
    #rename values
    surf_idx      = arg[0]
    
    #Fist moment of inertia
    xCoor     = RayTrace[:,8,surf_idx]
    yCoor     = RayTrace[:,9,surf_idx]
    N         = len(xCoor)
    centroidX = sum(xCoor) / N 
    centroidY = sum(yCoor) / N
    
    #Root mean square radius
    ms =np.sum(np.power(xCoor-centroidX,2)+np.power(yCoor-centroidY,2))
    mrs=np.sqrt(ms)
    
    error = mrs
    
    return error



if __name__=='__main__':
    from JenTrace.ray_src import PointSource
    from JenTrace.opt_sys import OpSysData
    from JenTrace.opt_dsg import OpDesign
    #from JenTrace.plt_fnc import plot_system,plot_rayTrace 

    #import matplotlib.pyplot as plt
    
    pto1  = PointSource([0,1,0],635)
    syst1 = OpSysData()
    syst1.add_surface(2,0.05,1.7)
    syst1.add_surface(-2,-0.5,1.4)
    #syst1.add_surface(10,0,1)
    syst1.plot([1.5])
    #syst1.changeAperture(1,surfIndex = 1)
    
    design1  = OpDesign(pto1,syst1,aprRad=0.66,aprInd=2)
    design1.plot()
    design1.autofocus()
    #print(design1.optSys)
    design1.plot()
    #fig, ax = plt.subplots()
    #fig, ax = plot_system(design1, fig=fig, ax=ax)
    #fig, ax = plot_rayTrace(design1.raySrcTrace,fig=fig,ax=ax)
    #fig, ax = plot_rayTrace(design1.dsgPtoTrace,fig=fig,ax=ax)
    #fig, ax = plot_rayTrace(design1.dsgInfTrace,fig=fig,ax=ax)
    
    
    