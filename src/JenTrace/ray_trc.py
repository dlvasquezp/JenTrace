# -*- coding: utf-8 -*-
"""
Created on Sat Mar 21 05:44:55 2020
@author: David Vasquez
functions:  -trace
            -fill_first
            -makeTrace           
            -print_report
            -makeTrace_paraxial *
            -makeTrace_grin *
* implemented in future versions
"""
try: import JenTrace
except ModuleNotFoundError: 
    import os, sys
    sys.path.insert(0,os.path.dirname(os.getcwd()))

from JenTrace.catalog import OPT_GLASS, sellmeierDispForm
import math
import numpy

def trace(RayList,SurfaceData):
    '''
    RayList    : list of lists, ([x,y,z],[cosX,cosY,cosZ], lambda)
    SurfaceData: list of lists, (d,C,n,surfType)
    RayTrace: numpy.ndarray  [d,c,n]            distance,curvature,refraction index
                             [X0,Y0]            last position X,Y
                             [F,G,Delta]        -
                             [X,Y,Z]            new position X,Y,Z
                             [check1]           -> Boolean
                             [alfa,beta,gamma]  last direction cosines 
                             [lambda]           Wavelength
                             [cosI,CosIp]       -
                             [K]                -
                             [L,M,N]            new direction cosines
                             [Check2]           -> Boolean
                             [surfType]         -> 1:Standard, 2:Paraxial
                         
              The position is represented by the index.
              In the array, i represent the ray, j the data and k the surface
            
    '''
    i = len(RayList)                          # Ray i
    j = 24                                    # [d,c,n][X0,Y0][F,G,Delta][X,Y,Z][check1][alfa,beta,gamma][lambda][cosI,CosIp][K][L,M,N][Check2][surfType]
    k = len(SurfaceData)                      # Surface k
    RayTrace = numpy.zeros([i,j,k])
    RayTrace = fill_first(RayList,SurfaceData,RayTrace,i,k)
    RayTrace = makeTrace(RayTrace,i,k)
    
    return RayTrace

def fill_first(RayList,SurfaceData,RayTrace,i,k):
    for q in range (0,i):                                       # i -> ray indice
        RayTrace[q,8 ,0] = RayList[q][0][0]                     # X coordinate from object
        RayTrace[q,9 ,0] = RayList[q][0][1]                     # Y coordinate from object
        RayTrace[q,10,0] = RayList[q][0][2]                     # Z coordinate from object
        
        RayTrace[q,11,0] = 1                                    # Chech 1 -> true
        RayTrace[q,15,:] = RayList[q][2]                        # Wavelength in nm
        
        RayTrace[q,19,0] = RayList[q][1][0]                     # X cosine director
        RayTrace[q,20,0] = RayList[q][1][1]                     # Y cosine director
        RayTrace[q,21,0] = RayList[q][1][2]                     # Z cosine director 
        
        RayTrace[q,22,0] = 1                                    # Chech 2 -> true
        RayTrace[q,23,0] = 1                                    # Surftype -> 1: Standard
        
        for w in range (0,k):                                  # k -> surface indice
            RayTrace[q,0,w] = SurfaceData[w][0]                # distance in mm
            RayTrace[q,1,w] = SurfaceData[w][1]                # curvature in 1/mm
            # refraction indice 
            if isinstance(SurfaceData[w][2],str):
                try:
                    RayTrace[q,2,w] = sellmeierDispForm (OPT_GLASS[SurfaceData[w][2]]
                                                        ,RayList[q][2])
                except KeyError:
                    raise ValueError('%s, non existent glass material' % SurfaceData[w][2])
            else:
                RayTrace[q,2,w] = SurfaceData[w][2] 
            # Surface type 
            if SurfaceData[w][3]=="standard":               
                RayTrace[q,23,w]=1
            if SurfaceData[w][3]=="paraxial":
                RayTrace[q,23,w]=2
            
    return RayTrace

def makeTrace(RayTrace,i,k):
    for q in range (0,i):                                             #i -> ray indice
        for w in range (1,k):                                         #k -> surface indice, start in 1
                
            #Transfer part1
            dmin1 = RayTrace[q,0 ,w-1]
            
            Xmin1 = RayTrace[q,8 ,w-1]
            Ymin1 = RayTrace[q,9 ,w-1]
            Zmin1 = RayTrace[q,10,w-1]
            
            Lmin1 = RayTrace[q,19,w-1]
            Mmin1 = RayTrace[q,20,w-1]
            Nmin1 = RayTrace[q,21,w-1]
            
            X0 = Xmin1 + (Lmin1/Nmin1)*(dmin1-Zmin1)
            Y0 = Ymin1 + (Mmin1/Nmin1)*(dmin1-Zmin1)
            
            #Transfer part2
            c = RayTrace[q,1,w]
            
            F = c*(X0**2+Y0**2)
            G = Nmin1 - c*(Lmin1*X0+Mmin1*Y0)
            
            try:
                Delta = F / (G + math.sqrt(G**2-c*F) )
            except ValueError:
                #Surface not found
                RayTrace[q,3:11,w:]  = numpy.NaN
                RayTrace[q,12:15,w:] = numpy.NaN
                RayTrace[q,16:22,w:] = numpy.NaN
                break 
            
            X = X0 + Lmin1*Delta
            Y = Y0 + Mmin1*Delta
            Z = Nmin1*Delta
            
            # The ray i meet the surface k
            check1 = 1
            
            #Refraction
            n   = RayTrace[q,2,w-1]
            np  = RayTrace[q,2,w]
            
            alfa = -c*X
            beta = -c*Y
            gamma= 1 - c*Z
            
            
            CosI  = math.sqrt(G**2 - c*F)
            
            try:
                CosIp = (1/np)*math.sqrt(np**2 - (n**2)*(1-CosI**2))
            except ValueError:
                #ray reflected instead of refracted
                RayTrace[q,3:11,w:]  = numpy.NaN
                RayTrace[q,12:15,w:] = numpy.NaN
                RayTrace[q,16:22,w:] = numpy.NaN
                break
            
            K = c*(np*CosIp - n*CosI) 
            
            L = (1/np)*(n*Lmin1 - K*X )
            M = (1/np)*(n*Mmin1 - K*Y )
            N = (1/np)*(n*Nmin1 - K*Z + np*CosIp - n*CosI)
            
            #Check1 and Check for direction cosines
            check2 = int(check1 and numpy.isclose((L**2+M**2+N**2),1))
            
            # Replace
            
            #RayTrace[q,0 ,w] = d                        # From SysData in function "fill_first"
            #RayTrace[q,1 ,w] = c                        # From SysData in function "fill_first"
            #RayTrace[q,2 ,w] = np                       # From SysData in function "fill_first"
            RayTrace[q,3 ,w] = X0
            RayTrace[q,4 ,w] = Y0
            RayTrace[q,5 ,w] = F
            RayTrace[q,6 ,w] = G
            RayTrace[q,7 ,w] = Delta
            RayTrace[q,8 ,w] = X
            RayTrace[q,9 ,w] = Y
            RayTrace[q,10,w] = Z
            RayTrace[q,11,w] = check1
            RayTrace[q,12,w] = alfa
            RayTrace[q,13,w] = beta
            RayTrace[q,14,w] = gamma
            #RayTrace[q,15,w] = Lambda                   # From sysData in function "fill_first"
            RayTrace[q,16,w] = CosI
            RayTrace[q,17,w] = CosIp
            RayTrace[q,18,w] = K
            RayTrace[q,19,w] = L
            RayTrace[q,20,w] = M
            RayTrace[q,21,w] = N
            RayTrace[q,22,w] = check2
            #RayTrace[q,23 ,w] = SurfType                # From SysData in function "fill_first"    
    
    return RayTrace

def print_report(RayTrace,info='ray',index=0):
    (i,j,k) = RayTrace.shape
    if j==24:
        headers = ["d","c","n","X0","Y0","F","G","Delta","X","Y","Z"
                      ,"check1","alfa","beta","gamma","lambda","cosI","CosIp"
                      ,"K","L","M","N","Check2","surfType"]
        if info == 'surf':
            surfInfo= RayTrace[:,:,index].T
            print("\nRAYTRACE REPORT - SURFACE #{:}\n".format(index))
            str_row = "{: >12} ".format("Prop\Ray#")
            for rayNo in range(i):
                str_row += "{: >12} ".format(rayNo)
            print(str_row)
            for prop in range(j):
                str_val = "{: >12} ".format(headers[prop])
                row = surfInfo[prop,:]
                for data in row:
                    str_val +=  "{:12f} ".format(data)   
                print(str_val)
                
        if info == 'ray':
            surfInfo= RayTrace[index,:,:]
            print("\nRAYTRACE REPORT - RAY #{:}\n".format(index))
            str_row = "{: >12} ".format("Prop\Surf#")
            str_row += "{: >12} ".format("Object")
            for surfNo in range(1,k-1):
                str_row += "{: >12} ".format(surfNo)
            str_row += "{: >12} ".format("Image")
            print(str_row)
            for prop in range(j):
                str_val = "{: >12} ".format(headers[prop])
                row = surfInfo[prop,:]
                for data in row:
                    str_val +=  "{:12f} ".format(data)   
                print(str_val)
                
        if info == 'prop':
            if index < 24:
                surfInfo= RayTrace[:,index,:]
                print("\nRAYTRACE REPORT - PROPERTY #{:}\n".format(headers[index]))
                str_row  = "{: >12} ".format("Ray#\Surf#")
                str_row += "{: >12} ".format("Object")
                for surfNo in range(1,k-1):
                    str_row += "{: >12} ".format(surfNo)
                str_row += "{: >12} ".format("Image")
                print(str_row)
                for ray in range(i):
                    str_val = "{: >12} ".format(ray)
                    row = surfInfo[ray,:]
                    for data in row:
                        str_val +=  "{:12f} ".format(data)   
                    print(str_val)
            else:
                raise Warning('Incorrect property index (0-23).')
    else:
        raise Warning('The RayTrace has an incorrect shape.')
        
 
if __name__ == '__main__':
    from ray_src import PointSource
    from opt_sys import OpSysData
    
    pto1  = PointSource([0,0.1,0],1040)
    syst1 = OpSysData()
    syst1.add_surface(10,0.005,'N-BK7')
    syst1.add_surface(10,5,1.42)
    syst1.add_surface(10,0.005,'N-K5')
    
    trace1 = trace(pto1.RayList,syst1.SurfaceData)
    
    print_report(trace1, 'prop',index=8)
    print_report(trace1, 'prop',index=9)
    
    print_report(trace1)
    print_report(trace1, 'ray',index=4)
    