# -*- coding: utf-8 -*-
"""
Created on Sat Apr 18 09:57:31 2020
@author: David Vasquez
Function:   -seidel_coef
            -plot_seidel
"""
try: import JenTrace
except ModuleNotFoundError: 
    import os, sys
    sys.path.insert(0,os.path.dirname(os.getcwd()))

import numpy 
import matplotlib.pyplot as plt

def seidel_coef(OptDsg):

    
    i   = 1                             # Principal Ray index
    i_b = 0                             # Chief Ray index
    j   = 5                             # 5 Aberrations
    k   = len(OptDsg.optSys.SurfaceData)# Surface k
    
    '''
    _p : sufix for prime
    _b : sufix for bar
    _pb: sufix for prime bar
    '''
    
    SeidelCoef = numpy.zeros((j,k))                            
    
    
    for k in range (1, k):             # k -> surface indice (start in 1)
        # Principal ray - surface n-1
        n     =  OptDsg.dsgPtoTrace[i,2 ,k-1]
        U     =  OptDsg.dsgPtoTrace[i,20,k-1] 
        h     =  OptDsg.dsgPtoTrace[i,9,k]     #h = h_p
        
        # Principal ray - surface n
        c     =  OptDsg.dsgPtoTrace[i,1,k]
        n_p   =  OptDsg.dsgPtoTrace[i,2 ,k]
        U_p   =  OptDsg.dsgPtoTrace[i,20,k]
        #I_p   =  math.acos(OptDsg.dsgPtoTrace[i,17,k])
        
        # Chief ray - surface n-1
        U_b   =  OptDsg.raySrcTrace[i_b,20,k-1] 
        h_b   =  OptDsg.raySrcTrace[i_b,9,k]  #h_b = h_pb
        
        # Chief ray - surface n
        #n_pb  =  OptDsg.raySrcTrace[i_b,2 ,k]
        #I_pb  =  math.acos(OptDsg.raySrcTrace[i_b,17,k])
        
        # Abreviations
        A     =  n*(h  *c+U  )
        A_b   =  n*(h_b*c+U_b)  
        H     =  n*(U*h_b-U_b*h) #Lagrange invariant      
            
        #Spheric aberration
        SeidelCoef[0,k] = -( (A**2)  * h * ((U_p/n_p)- (U/n) ))
        #Coma
        SeidelCoef[1,k] = -( (A*A_b) * h * ((U_p/n_p)- (U/n) ))
        #Astigmatism
        SeidelCoef[2,k] = -( (A_b**2)* h * ((U_p/n_p)- (U/n) ))
        #Petzval curvature
        SeidelCoef[3,k] = -( (H**2)  * c * (  (1/n_p)- (1/n) ))
        #Distortion
        SeidelCoef[4,k] = -( ((A_b**3)/A)     * h * ((U_p/n_p)- (U/n)  )  
                           + (A_b/A) * (H**2) * c * (  (1/n_p)- (1/n)  ))
        
    return SeidelCoef

def plot_seidel(SeidelCoef):
    (i,NoSurf)=SeidelCoef.shape
    SeidelSum =numpy.sum(SeidelCoef,1)
    # Bar width and position
    width = 0.1
    posRel = numpy.linspace(-2.5*width,+2.5*width,6)
    # Surfaces label
    labels = [str(q) for q in range(1,NoSurf-1)]
    labels.append('SUM')
    # Aberration label
    aber_label = ['Spherical','Coma','Astigmatism','Petzval curvature','Distortion']
    
    plt.figure()
    plt.title('Seidel Coefficients')
    for seidel in range(i):
        plt.bar(range (NoSurf-1)+posRel[seidel], [*SeidelCoef[seidel,1:-1],SeidelSum[seidel]], width=width,label=aber_label[seidel])
    for pos    in range(NoSurf-2):
        plt.axvline(x=pos+0.5, color='black',ls='--')
    plt.xticks(range (NoSurf-1), labels) 
    plt.grid(True,axis='y')  
    plt.legend(bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)
    plt.show()
        
    

if __name__=='__main__':
    from opt_sys import OpSysData
    from ray_src import PointSource
    from opt_dsg import OpDesign

    # Instantiate optical system
    syst1 = OpSysData()
    syst1.change_surface(10    ,0         ,1      ,surfIndex=0) 
    syst1.add_surface   (3.50  ,1/15.37   ,'N-BK7')
    syst1.add_surface   (1.50  ,1/-11.10  ,'N-SF5')
    #syst1.add_surface   (5     ,1/-31.47  ,1      )
    #syst1.add_surface   (5     ,1/-3.47   ,'N-SF5')
    syst1.add_surface   (5     ,0         ,1      )
    #syst1.changeAperture(1,surfIndex=2)
    clearSemDia=[1,5.0,5.0,5.0,2]#,2]
    syst1.plot_optical_system(clearSemDia)
    
    # Instantiate point source
    pto1  = PointSource([0,1.0,0],635)
    
    # Instantiate optical design
    design1  = OpDesign(pto1,syst1,aprInd=3)
    
    SCoef =  seidel_coef(design1)
    print(SCoef)
    
    # plot design
    design1.plot_design(clearSemDia)
    plot_seidel(SCoef)