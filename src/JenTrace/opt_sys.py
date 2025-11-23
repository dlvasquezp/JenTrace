# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 19:29:17 2020
@author: David Vasquez
class OpSysData
"""
try: import JenTrace
except ModuleNotFoundError: 
    import os, sys
    sys.path.insert(0,os.path.dirname(os.getcwd()))
    
import matplotlib.pyplot as plt
from JenTrace.plt_fnc import plot_system
from JenTrace.catalog import OPT_GLASS

class OpSysData:
    '''
    Optical system data storing the physical data of the optical surfaces.
    Attributes:
        surfIndex: Surface position in the sequential optical system
        SurfaceData: list of lists, (d,C,n,surfType,clearSemDia,varProp)
        d: distance (mm)
        C: curvature (1/mm)
        n: refraction index (-)
        surfType: surface type ('standard', 'paraxial', 'grin', 'aspheric')
        surfPara: surfaces parameter, need for the complete description of the surface such as gradient profile or aspheric curvatures (Pending). 
        clearSemDia: clear semi-diameter (mm)
        varProp: String indicating the property (d,C,n) to be considered variable in the merit function. 
                d,C,n = ###,
                Example 1: '100' d is variable and C and n are fixed. 
                Example 2: '110' d and C are variable and n is fixed.
                Example 3: '000' d, C and n are fixed.
    '''
    
    surfaceTypes = {'standard', 'paraxial', 'grin', 'aspheric'}
    
    def __init__(self):                      
        self.SurfaceData=[]
        self.add_surface(10,0,1,'standard',surfIndex=0) #object
        self.add_surface( 0,0,1,'standard',surfIndex=1) #image
    
    def add_surface(self, d:(int,float), C:(int,float), n:(int, float, str), surfType:str='standard', clearSemDia:float=10.0, surfPara:list=[], varProp:str='000', surfIndex:int=-1):
        self.check_arguments(d=d,C=C,n=n,surfType=surfType)
        surf_len= len(self.SurfaceData)
        
        if surfIndex >= 0 and surfIndex <= surf_len:
            self.SurfaceData.insert(surfIndex,[d,C,n,surfType,clearSemDia,surfPara,varProp]) 
        else:
            if surfIndex == -1:
                #if surfIndex not defined, append to the last position before 'image'
                self.SurfaceData.insert((surf_len-1),[d,C,n,surfType,clearSemDia,surfPara,varProp]) 
            else:
                raise ValueError('addSurface: invalid surface index')
                
                
    def change_surface(self, d:(int,float), C:(int,float), n:(int, float, str), surfType='standard', clearSemDia:float=10.0, surfPara:list=[], varProp:str='000', surfIndex:int=-1):
        self.check_arguments(d=d,C=C,n=n,surfType=surfType)
        surf_len= len(self.SurfaceData)
        
        if surfIndex >= 0 and surfIndex <= surf_len:
            self.SurfaceData[surfIndex]=[d,C,n,surfType,clearSemDia,surfPara,varProp] 
        else:
            raise ValueError('changeSurface: invalid surface index ')
            
            
    def delete_surface(self,surfIndex=-1):
        surf_len= len(self.SurfaceData)
        
        if surfIndex >= 0 and surfIndex < surf_len:
            del self.SurfaceData[surfIndex]
        else:
            if surfIndex == -1:
                del self.SurfaceData[(surf_len-2)]    
            else:
                raise ValueError('deleteSurface: invalid surface index')
            
            
    def invert_surface_order(self, surf1, surf2):                   
        surf_len= len(self.SurfaceData)
        
        assert surf1 >= 0       ,     'invertSurface: surf1 invalid'
        assert surf2 <  surf_len,     'invertSurface: surf2 invalid'
        assert surf1 <  surf2   ,     'invertSurface: surf1 > surf2'
        
        surfCopy_inv = self.SurfaceData[surf1:surf2+1][::-1]
        dist = [d[0] for d in surfCopy_inv]
        curv = [C[1] for C in surfCopy_inv]
        refI = [n[2] for n in surfCopy_inv]
        len_surf = len(surfCopy_inv)
        
        for q in range(len_surf):
            surfCopy_inv[q][0]= +dist[(q+1)%len_surf]
            surfCopy_inv[q][1]= -curv[q]
            surfCopy_inv[q][2]= refI[(q+1)%len_surf]
        
        self.SurfaceData[surf1:surf2+1] = surfCopy_inv
        
    def check_arguments(self,d='nan',C='nan',n='nan',surfType='nan'):
        if d!='nan':
            assert isinstance(d,(int,float)), 'distance [d] must be either int of float'
        if C!='nan':
            assert isinstance(C,(int,float)), 'curvature [C] must be either int of float'
        if n!='nan':
            assert isinstance(n,(int,float,str)), 'refraction index [n] must be either int,float or str'
            if isinstance (n,str):
                assert n in OPT_GLASS, '{} is not found in the optical glass catalog. The available materials are: {} '.format(n, [k for k in OPT_GLASS.keys()])
        if surfType!='nan':
            assert surfType in self.surfaceTypes, 'Surface type [surfType] not supported'
 
    def __str__(self):
        headers = ['#','Distance','Curvature','Material','Type', 'SemiDiameter','Parameters']
        print("\nSURFACE LIST")
        print("Class:" + self.__class__.__name__)
        print("{: >5}{: >12}{: >12}{: >12}{: >12}{: >13}  {:12}".format(*headers))
        for idx, surface in enumerate(self.SurfaceData):
            print("{: 5d}".format(idx)+
                  "{:12.3f}{:12.3f}{:>12}{:>12}{:>13.2f}".format(*surface[:5])+
                  "  {:12}".format(str(surface[5]))
                  )
        return str('#####')
    
    def change_clearSemDia(self,clearSemDia_usr:list[float]):
        if   len(clearSemDia_usr)== 1:
            assert clearSemDia_usr[0] > 0,'Invalid value of clearSemDia: {}'.format(clearSemDia_usr[0])
            for idx, surface in enumerate (self.SurfaceData):
                surface[4] = clearSemDia_usr[0]
                self.SurfaceData[idx]=surface
        elif len(clearSemDia_usr)== len(self.SurfaceData):
            assert (all([isinstance(q,(int,float)) for q in clearSemDia_usr])),'Invalid data type in clearSemDia: {}'.format(clearSemDia_usr)
            assert (all([ q>0 for q in clearSemDia_usr])),'Negative clearSemDia value: {}'.format(clearSemDia_usr)
            
            for idx, surface in enumerate (self.SurfaceData):
                surface[4] = clearSemDia_usr[idx]
                self.SurfaceData[idx]=surface
        else:
            print('Unmatch clearSemDia and surface dimension. # Surfaces: {}, # Semi Diameters: {}'.format(len(self.SurfaceData),len(clearSemDia_usr)))
            
    def plot(self,clearSemDia_usr=[]):
        if clearSemDia_usr != []:
            self.change_clearSemDia(clearSemDia_usr)
        fig, ax = plt.subplots()
        print(self)
        plot_system(self,fig=fig,ax=ax)
               
if __name__=='__main__':
    #Create system
    syst1 = OpSysData()
    print('\nDefault Optical system')
    #syst1.print_report()
    print(syst1)
    
    syst1.change_surface(2,0,1.1,surfIndex=0)
    syst1.add_surface(2,1.0,2)
    syst1.add_surface(2,1.0,2,surfType='paraxial',surfPara=[1,2,3,4,5])
    syst1.add_surface(10,1/2.0,2)
    syst1.add_surface(3,1/3.0,'N-BK7')
    syst1.add_surface(4,1/4.0,1)
    syst1.add_surface(4,1/4.0,1.5,surfPara=['test'])
    print('\nOptical system data')
    #syst1.print_report()
    print(syst1)
    
    #Test changeSurface
    syst1.change_surface(9,1/9,4,surfIndex=3)
    print('\nSurface 3 changed')
    print(syst1)
    
    #Test deleteSurface
    syst1.delete_surface(surfIndex=3)
    print('\nSurface 3 deleted')
    print(syst1)
    #Plot system
    syst1.change_clearSemDia(clearSemDia_usr=[4])
    syst1.plot()
    
    #Test invert surface
    syst1.invert_surface_order(0,5)
    print('\nSurfaces 0 to 5 inverted')
    print(syst1)
    #Plot system
    syst1.plot()
    syst1.plot(clearSemDia_usr=[1.5])
    syst1.plot(clearSemDia_usr=[3,3,2,1.5,2,2,3])
    
    #Test invalid clearSemDia
    print('---Invalid data type---')
    try:
        syst1.plot(clearSemDia_usr=[3,'3',2,1.5,2,2,3])
    except Exception as e: print(e)
    
    print('---Negative value---')
    try:
        syst1.plot(clearSemDia_usr=[3,3,2,-1.5,2,2,3])
    except Exception as e: print(e)
        
    print('---Unmatch dimensions---')
    try:
        syst1.plot(clearSemDia_usr=[3,5,3,2,4,1.5,2,2,3])
    except Exception as e: print(e)
    

    
