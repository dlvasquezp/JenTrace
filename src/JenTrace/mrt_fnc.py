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
from JenTrace.spt_dgm import ray_pattern, spot_diagram
from JenTrace.ray_src import RaySource#calc_direcCos
import numpy as np
import copy

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
    assert arg[0].__class__.__name__=='OpDesign' ,'arg[0] must be an OpDesign object, instead of {} '.format(type(arg[0])) 
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
    print(vector)
    LMN = RaySource.calc_direcCos([vector[0],vector[1],1])  
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
    assert arg[0].__class__.__name__=='OpDesign' ,'arg[0] must be an OpDesign object, instead of {} '.format(type(arg[0]))  
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
    print('########## Merit x0: {}#####'.format(x0))
    print(opticalDesign.aprInd,opticalDesign.aprRad)
    print(opticalDesign.optSys)
    opticalDesign.optSys.set_varValues(x0)
    print(opticalDesign.optSys)
    opticalDesign.solve_dsg(caller = 'minimization')
    print('x0:', x0)
    opticalDesign.plot
    #Generate RaySources
    raySourceList=[]
    for param in raysParam:
        #print([opticalDesign,*param])
        raySourceList.append(ray_pattern(opticalDesign,*param))
    #print(raySourceList[0])    
    #Make Raytrace
    rayTraceList = []
    for source in raySourceList:
        rayTraceList.append(trace(source.RayList,opticalDesign.optSys.SurfaceData))
        
    #Evaluate error functions
    errorList = []
    for fun, rayTrace, param , w in zip(errorFun, rayTraceList, funParam,weights):
        #print(rayTrace)
        errorList.append(fun(rayTrace,*funParam)*w)
    
    error = np.sum(errorList)
    
    return error

def spot_sqrt (RayTrace,*arg)->[float]:
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

def dfdx(fun,x0,dx=[1e-6], *arg):
    '''
    Finite difference differentiation
    Parameters
    ----------
    fun : TYPE
        DESCRIPTION.
    x0 : TYPE
        DESCRIPTION.
    dx : list(floats)
        Initial step. (broadcasteable)
    *arg : TYPE
        DESCRIPTION.

    Returns
    -------
    int
        DESCRIPTION.

    '''
    assert len(x0)==len(dx), 'Missmatch in x0 and dx'
    assert any(dx), 'dx equal 0'
    
    dim= len(x0)
    yp = []
    y0 = fun(x0,*arg)
    
    for i, dx_i in enumerate(dx):
        h = np.zeros(dim)
        h[i] = dx_i
        
        y1 = fun(x0+h,*arg)
        
        
        if np.isnan(y1) or np.isnan(y0):
            return False
        else:
            yp_i = (y1 - y0) / np.sum(h)
            yp.append(yp_i)
            
    return yp, y0

def gradient_descent (fun,x0,fvalue=5e-3,eps=1e-5, eta=2e-3,*arg):
    
    # dx will be the change in the solution -- we'll iterate until this
    # is small
    dx = [1e-6,1e-6]
    xp_old = x0.copy()

    #y = fun(xp_old,*arg)
    grad, y = dfdx(fun,xp_old,dx,*arg)

    dx_norm=1
    while (dx_norm > eps and y >fvalue):
        print(xp_old,eta,grad)
        xp = xp_old - np.multiply(eta,grad)
        
        #y = fun(xp,*arg)
        grad, y = dfdx(fun,xp,dx, *arg)
        
        dx_norm = np.linalg.norm(xp - xp_old)
        xp_old = xp.copy()
    
    return xp

def calc_direcCos_2(vector):
    assert (all([isinstance(q,(int,float)) for q in vector]) 
            and len(vector) == 3),'Invalid vector'
    #Calculate direction cosines
    norm     = np.sqrt(np.sum(np.power(vector,2))) 
    cosDirX  = vector[0]/norm
    cosDirY  = vector[1]/norm
    cosDirZ  = vector[2]/norm 
    
    return [cosDirX,cosDirY,cosDirZ]
    

def LMN_apertureStop_2 (x0,*arg):
    '''
    Merit function to calculate the error produced by the rays propaged in the 
    direction [x0[0],x0[1],1] and the aperture stop
    
    x0:   (list) [vector_compX, vector_comp_y] 
    arg: (list) [object Optical desing, object ray source, int indexRay]
    
    arg[0] must be Optical design
    arg[1] must be point source (pointSource)
    arg[2] must be int [0,1,2,3,4]
    '''   
    assert arg[0].__class__.__name__=='OpDesign' ,'arg[0] must be an OpDesign object, instead of {} '.format(type(arg[0])) 
    assert arg[1].__class__.__name__=='PointSource' ,'arg[1] must be a PointSource object' 
    assert arg[2] >=0 and arg[2] <=4, 'Invalid ray index (indexRay)'
    #reference the values
    #vectorList    = np.asarray(x0) 
    
    #print("x.shape =", np.shape(x0))
    m, batch = x0.shape[0], x0.shape[1:]
    #print('m, batch: ',m,batch)
    x = np.reshape(x0, (m, -1)) 
    
    
    SurfaceData   = arg[0].optSys.SurfaceData
    ApertureRadio = arg[0].aprRad
    ApertureIndex = arg[0].aprInd
    ptoSrc        = copy.deepcopy(arg[1])
    #RayList       = (arg[1].RayList).copy
    indexRay      = arg[2]
    
    #print(x[0],x[1])
    errorList=[]
    for vecX, vecY in zip(x[0],x[1]):
        #print(vecX,vecY)
        #Calculate direction cosines
        #print(vectorList)
        LMN = ptoSrc.calc_direcCos([vecX,vecY,1])  
        #replace cosine director    
        ptoSrc.change_LMN(LMN,indexRay)
        #Make Raytrace
        RayTrace  = trace(ptoSrc.RayList,SurfaceData)
        #Calculate error
        error = ap_error_calc(RayTrace,ApertureRadio,indexRay,ApertureIndex)
        errorList.append([error])
            
    #res = np.array(errorList).reshape(len(vectorList),)
    res = np.array(errorList)
    #return res

    #n = res.shape[0]
    #print(res) 
    #print( '#########',(1,) + batch)
    res2 = np.reshape(res, (1,) + batch) # return shape (2, ...)
    #print(res2) 
    #print("y.shape =", np.shape(res2))
    
    return res2
    
    '''
    #print(x0, vectorList.ndim)
    print("x.shape =", x0.shape)
    if vectorList.ndim == 1:
        vector = vectorList
        #Calculate direction cosines
        #print(vector)
        #print(vector[0], type[vector[0]])
        #print(vector[1], type[vector[1]])
        
        LMN = ptoSrc.calc_direcCos([vector[0],vector[1],1])  
        #replace cosine director    
        ptoSrc.change_LMN(LMN,indexRay)
        #Make Raytrace
        RayTrace  = trace(ptoSrc.RayList,SurfaceData)
        #Calculate error
        error = ap_error_calc(RayTrace,ApertureRadio,indexRay,ApertureIndex)
        res = np.array(error).reshape(1,)
        #print(np.shape(x0),np.shape(res))
        print("y.shape =", np.shape(res))
        return res
    
    if vectorList.ndim == 2:
        #Calculate direction cosines
        #print(vectorList)
        errorList=[]
        for vector in vectorList:
            LMN = ptoSrc.calc_direcCos([vector[0],vector[1],1])  
            #replace cosine director    
            ptoSrc.change_LMN(LMN,indexRay)
            #Make Raytrace
            RayTrace  = trace(ptoSrc.RayList,SurfaceData)
            #Calculate error
            error = ap_error_calc(RayTrace,ApertureRadio,indexRay,ApertureIndex)
            errorList.append([error])
        res = np.array(errorList).reshape(len(vectorList),)
        print(res) 
        print("y.shape =", np.shape(res))
        return res
    '''


if __name__=='__main__':
    import time
    from scipy.optimize import minimize
    from JenTrace.ray_src import PointSource
    from JenTrace.opt_sys import OpSysData
    from opt_dsg import OpDesign
    
    from scipy.optimize import approx_fprime
    
    from JenTrace.plt_fnc import plot_system,plot_rayTrace 
    from JenTrace.ray_trc import print_report
    from scipy.differentiate import jacobian
    from scipy.optimize import lsq_linear
    #%%
    print('######### Dsg 0 ##########')
    pto1  = PointSource([0,1,0],635)
    syst1 = OpSysData()
    syst1.change_surface(2,0,1,surfIndex=0)
    syst1.add_surface(12.28,0.05,1.7,varProp='100')
    syst1.add_surface(5,-0.05,1.4)
    syst1.add_surface(9.910, -0.386,1,varProp='110')
    syst1.plot([3])
    
    design1  = OpDesign(pto1,syst1,aprRad=1,aprInd=2)
    qwe = np.copy(design1.usrSrc)
    arg=(design1,design1.usrSrc,3)
    x0=[0,0]
    LMN_apertureStop(x0,*arg)
    dfdx(LMN_apertureStop,x0,[1e-6,1e-6], *arg)
    
    start_time = time.time()
    res=gradient_descent (LMN_apertureStop,x0,5e-3,1e-5,2e-3,*arg)
    print("--- %s seconds ---" % (time.time() - start_time))
    
    print(qwe)
    print(res)
    
    
    f_wrapped = lambda x: LMN_apertureStop_2(x,*arg)
    
    x0=np.array([0.1,0.2])
    print(x0.shape,'>>>>',f_wrapped(x0))
    
    #x0=np.array([[0,0]])
    #print(x0.shape,'>>>>',f_wrapped(np.array(x0)))
    
    #x0=np.array([[0,0],[0.1,0.1],[0,0]])
    #print(x0.shape,'>>>>',f_wrapped(np.array(x0)))
    
    start_time = time.time()
    print(x0)
    for _ in range(10):
        res = jacobian(f_wrapped, x0)
        #print(res.df,res.df.shape)    
        J = res.df
        A = np.matmul(J.T,J)
        b = np.matmul(-J.T,f_wrapped(x0))
        #b = np.array([ 0, 0])
        res2 = lsq_linear(A,b)
        x0 = res2.x
        print(x0)
    print("--- %s seconds ---" % (time.time() - start_time))
    
    
    #x0=np.array([[0,0],[0.1,0.1],[0,0]])
    
    
    

    
    #%%
    #import matplotlib.pyplot as plt
    print('######### Dsg 1 ##########')
    pto1  = PointSource([0,1,0],635)
    syst1 = OpSysData()
    syst1.change_surface(2,0,1,surfIndex=0)
    syst1.add_surface(12.28,0.05,1.7,varProp='100')
    syst1.add_surface(5,-0.05,1.4)
    syst1.add_surface(9.910, -0.386,1,varProp='110')
    syst1.plot([3])
    
    design1  = OpDesign(pto1,syst1,aprRad=1,aprInd=2)
    design1.plot()
    design1.autofocus()
    print(design1.optSys)
    #design1.plot()
    #spot_diagram(design1,noRings=9,show=True)
    
    
    #fig, ax = plt.subplots()
    #fig, ax = plot_system(design1, fig=fig, ax=ax)
    #fig, ax = plot_rayTrace(design1.raySrcTrace,fig=fig,ax=ax)
    #fig, ax = plot_rayTrace(design1.dsgPtoTrace,fig=fig,ax=ax)
    #fig, ax = plot_rayTrace(design1.dsgInfTrace,fig=fig,ax=ax)
    #%%   
    import matplotlib.pyplot as plt 
    
    print('######### Dsg 2 ##########')
    pto2  = PointSource([0,1,0],635)
    syst2 = OpSysData()
    syst2.change_surface(2,0,1,surfIndex=0)
    syst2.add_surface(47,0.05,1.7,varProp='100')
    syst2.add_surface(5,-0.05,1.4)
    syst2.add_surface(20.586,-0.162,1,varProp='110')
    syst2.plot([3])
    
    
    design2  = OpDesign(pto2,syst2,aprRad=1,aprInd=2)
    design2.plot()
    spot_diagram(design2,noRings=9,show=True)
    #print_report(design2.initRayTrace, index=0)
    #print(design2.usrSrc)
    #print(design2.dsgPtoSrc)
    #print(design2.dsgInfSrc)
    
    fig, ax = plt.subplots()
    fig, ax = plot_system(design2,fig,ax)
    fig, ax = plot_rayTrace(design2.initRayTrace,fig=fig,ax=ax, color='r')
    
    print_report(design2.initRayTrace)
    
    #%%   
    print('######### Dsg 3 ##########')
    pto3  = PointSource([0,1,0],635)
    syst3 = OpSysData()
    syst3.change_surface(2,0,1,surfIndex=0)
    syst3.add_surface(40,0.05,1.7,varProp='100')
    syst3.add_surface(5,-0.05,1.4)
    syst3.add_surface(19.950,-0.2,1,varProp='110')
    syst3.plot([3])
    #syst1.changeAperture(1,surfIndex = 1)
    
    
    #Merit function
    #opticalDesign  = arg[0]
    #raysParam      = arg[1]
    #errorFun       = arg[2]
    #funParam       = arg[3]
    #weights        = arg[4]
    print('-------------------->>>>>>>>>>>>>')
    start = time.time()
    design3  = OpDesign(pto3,syst3,aprRad=1,aprInd=2)
    #design3.autofocus()
    design3.plot()
    args = (design3,
           [['hexapolar',9]],
           [spot_sqrt],
           [[-1]],
           [1])
    x0 = design3.optSys.get_varValues()
    
    #eps = np.sqrt(np.finfo(float).eps)
    #fprime = approx_fprime(x0, lambda u:score_function(u, *args), eps)
    
    def gradient(x, *args):
        eps = np.sqrt(np.finfo(float).eps)
        # wrapper so approx_fprime sees a function f(x)
        return approx_fprime(x, lambda u:score_function(u, *args), eps)
    
    bounds=[(0,100),(0,100),(-0.3,0.3)]
    res = minimize(score_function, x0, args=args, method='Nelder-Mead',bounds=bounds)
    #res = minimize(score_function, x0, jac=gradient, args=args, method='Newton-CG',bounds=bounds)
    #res = minimize(score_function, x0, args=args, method='L-BFGS-B',bounds=bounds)
    #res = minimize(score_function, x0, args=args, method='Powell',bounds=bounds)
    
    design3.solve_dsg()
    print(res)  
    design3.plot()
    spot_diagram(design3,noRings=9,show=True)
    end = time.time()
    print('-------------------->>>>>>>>>>>>>',end - start)
    
    
    
    
    