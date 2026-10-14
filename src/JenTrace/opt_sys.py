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
import JenTrace.opt_sfc as optElement

class OpSysData:
    
    surfaceTypes = {'standard', 'paraxial', 'grin', 'aspheric'}
    
    def __init__(self):                      
        self.surfaceData=[]
        self.add_surface(optElement.standard(thickness= 10, curvature=0,material=1),surfIndex=0) #object
        self.add_surface(optElement.standard(thickness=  0, curvature=0,material=1),surfIndex=1) #image
                       
    def add_surface(self, surface, surfIndex:int=-1):
        surf_len= len(self.surfaceData)
        if not isinstance(surface, optElement.OpSurface):
            raise TypeError(
                "surface must be an OpticalSurface"
            )
        if surfIndex >= 0 and surfIndex <= surf_len:
            self.surfaceData.insert(surfIndex,surface) 
        else:
            if surfIndex == -1:
                #if surfIndex not defined, append to the last position before 'image'
                self.surfaceData.insert((surf_len-1),surface) 
            else:
                raise ValueError('addSurface: invalid surface index')               
                
    def change_surface(self, surface, surfIndex:int=-1):
        surf_len= len(self.surfaceData)
        
        if not isinstance(surface, optElement.OpSurface):
            raise TypeError(
                "surface must be an OpticalSurface"
            )
        
        if surfIndex == 0:
            surface.surfType  ='Standard'
            surface.curvature = 0
            self.surfaceData[surfIndex] = surface
        elif surfIndex >= 0 and surfIndex <= surf_len:
            self.surfaceData[surfIndex] = surface
        else:
            raise ValueError('changeSurface: invalid surface index ')   
            
    def delete_surface(self,surfIndex=-1):
        surf_len= len(self.surfaceData)
        
        if surfIndex >= 0 and surfIndex < surf_len:
            del self.surfaceData[surfIndex]
        else:
            if surfIndex == -1:
                del self.surfaceData[(surf_len-2)]    
            else:
                raise ValueError('deleteSurface: invalid surface index')
            
            
    def invert_surface_order(self, surf1, surf2):                   
        surf_len= len(self.surfaceData)
        
        assert surf1 > 0        ,    'invertSurface: surf1 invalid'
        assert surf2 <  surf_len,    'invertSurface: surf2 invalid'
        assert surf1 <  surf2   ,    'invertSurface: surf1 > surf2'
        
        surfCopy_inv = self.surfaceData[surf1:surf2+1][::-1]
        dist = [d.thickness for d in surfCopy_inv]
        curv =[]
        for surface in surfCopy_inv:
            if isinstance(surface, optElement.standard):
                curv.append(surface.curvature)
            elif isinstance(surface, optElement.paraxial):
                curv.append(surface.focalLength)
            elif isinstance(surface, optElement.grin):
                raise TypeError ("Grin surface not supported in surface order invertion")
            else:
                raise TypeError ("invert_surface_order: Surface type not supported")
                
        refI = [n.material for n in surfCopy_inv]
        len_surf = len(surfCopy_inv)
        
        for q in range(len_surf):
            surfCopy_inv[q].thickness = +dist[(q+1)%len_surf]
            if isinstance(surface, optElement.standard):
                surfCopy_inv[q].curvature= -curv[q]
            elif isinstance(surface, optElement.paraxial):
                surfCopy_inv[q].focalLength= +curv[q]
            else:
                raise TypeError ("invert_surface_order: Surface type not supported")
            surfCopy_inv[q].material = refI[(q+1)%len_surf]
        
        self.surfaceData[surf1:surf2+1] = surfCopy_inv
        
    def check_arguments(self,d='nan',C='nan',n='nan',surfType='nan'):
        if d!='nan':
            assert isinstance(d,(int,float)), 'thickness [d] must be either int of float'
        if C!='nan':
            assert isinstance(C,(int,float)), 'curvature [C] must be either int of float'
        if n!='nan':
            assert isinstance(n,(int,float,str)), 'refraction index [n] must be either int,float or str'
            if isinstance (n,str):
                assert n in OPT_GLASS, '{} is not found in the optical glass catalog. The available materials are: {} '.format(n, [k for k in OPT_GLASS.keys()])
        if surfType!='nan':
            assert surfType in self.surfaceTypes, 'Surface type [surfType] not supported'
 
    def fmt(self, value, name, varProp,spec):
        text = format(value, spec)
        if name in varProp:
            text = text[1:]+'*'
        return  text
 
    def __str__(self):
        headerStandard = ['#','Type','Thickness','Curvature'  ,'Material', 'SemiDiameter']
        
        print("\nSURFACE LIST")
        print("Class:" + self.__class__.__name__)
        print("{: >5}{: >12}{: >12}{: >12}{: >12}{: >13}".format(*headerStandard))
        actualHeader ='standard'
        for idx, surface in enumerate(self.surfaceData):
            if isinstance(surface,optElement.standard):
                if not actualHeader == 'standard':
                    headerStandard = ['','','','','', '']
                    print("{: >5}{: >12}{: >12}{: >12}{: >12}{: >13}".format(*headerStandard))
                
                print(
                    f"{idx:5d}"
                    f"{surface.surfType:>12}"
                    f"{self.fmt(surface.thickness, 'thickness', surface.varProp ,'12.3f')}"
                    f"{self.fmt(surface.curvature, 'curvature', surface.varProp ,'12.3f')}"
                    f"{self.fmt(surface.material, 'material', surface.varProp ,'>12')}"
                    f"{surface.clearSemDia:13.2f}"
                    )
                
                actualHeader ='standard'
                
            if isinstance(surface,optElement.paraxial):
                if not actualHeader == 'paraxial':
                    headerParaxial = ['','','','','','','FocalLength']
                    print("{: >5}{: >12}{: >12}{: >12}{: >12}{: >13}{: >12}".format(*headerParaxial))
                
                print(
                    f"{idx:5d}"
                    f"{surface.surfType:>12}"
                    f"{self.fmt(surface.thickness, 'thickness', surface.varProp ,'12.3f')}"
                    f"{'-':>12}"
                    f"{self.fmt(surface.material, 'material', surface.varProp ,'>12')}"
                    f"{surface.clearSemDia:13.2f}"
                    f"{self.fmt(surface.focalLength, 'focalLength', surface.varProp ,'12.3f')}"
                    )
                
                actualHeader ='paraxial'
                
            if isinstance(surface,optElement.grin):
                if not actualHeader == 'grin':
                    headerGrin = ['','','','','', '','B','C','K0','K1','K2']
                    print("{: >5}{: >12}{: >12}{: >12}{: >12}{: >13}{: >6}{: >6}{: >6}{: >6}{: >6}".format(*headerGrin))
                
                print(
                    f"{idx:5d}"
                    f"{surface.surfType:>12}"
                    f"{self.fmt(surface.thickness, 'thickness', surface.varProp ,'12.3f')}"
                    f"{self.fmt(surface.curvature, 'curvature', surface.varProp ,'12.3f')}"
                    f"{'-':>12}"
                    f"{surface.clearSemDia:13.2f}"
                    f"{self.fmt(surface.B, 'B',  surface.varProp ,'6.2f')}"
                    f"{self.fmt(surface.C, 'C',  surface.varProp ,'6.2f')}"
                    f"{self.fmt(surface.K0,'K0', surface.varProp ,'6.2f')}"
                    f"{self.fmt(surface.K1,'K1', surface.varProp ,'6.2f')}"
                    f"{self.fmt(surface.K2,'K2', surface.varProp ,'6.2f')}"
                    
                    )
                
                actualHeader ='grin'
                
            if isinstance(surface,optElement.aspheric):
                if not actualHeader == 'aspheric':
                    headerAspheric = ['','','','','', '','K','A4','A6','A8','A10','A12','A14']
                    print("{: >5}{: >12}{: >12}{: >12}{: >12}{: >13}{: >6}{: >8}{: >8}{: >8}{: >8}{: >8}{: >8}".format(*headerAspheric))
                
                print(
                    f"{idx:5d}"
                    f"{surface.surfType:>12}"
                    f"{self.fmt(surface.thickness, 'thickness', surface.varProp ,'12.3f')}"
                    f"{self.fmt(surface.curvature, 'curvature', surface.varProp ,'12.3f')}"
                    f"{self.fmt(surface.material, 'material', surface.varProp ,'>12')}"
                    f"{surface.clearSemDia:13.2f}"
                    f"{self.fmt(surface.K  ,  'B', surface.varProp ,'6.2f')}"
                    f"{self.fmt(surface.A4 , 'A4', surface.varProp ,'8.1E')}"
                    f"{self.fmt(surface.A6 , 'A6', surface.varProp ,'8.1E')}"
                    f"{self.fmt(surface.A8 , 'A8', surface.varProp ,'8.1E')}"
                    f"{self.fmt(surface.A10,'A10', surface.varProp ,'8.1E')}"
                    f"{self.fmt(surface.A10,'A12', surface.varProp ,'8.1E')}"
                    f"{self.fmt(surface.A10,'A14', surface.varProp ,'8.1E')}"
                    )
                
                actualHeader ='aspheric'
                

        return str('#####')
    
    def change_clearSemDia(self,clearSemDia_usr:list[float]):
        if   len(clearSemDia_usr)== 1:
            assert clearSemDia_usr[0] > 0,'Invalid value of clearSemDia: {}'.format(clearSemDia_usr[0])
            for idx, surface in enumerate (self.surfaceData):
                surface.clearSemDia = clearSemDia_usr[0]
                self.surfaceData[idx]=surface
        elif len(clearSemDia_usr)== len(self.surfaceData):
            assert (all([isinstance(q,(int,float)) for q in clearSemDia_usr])),'Invalid data type in clearSemDia: {}'.format(clearSemDia_usr)
            assert (all([ q>0 for q in clearSemDia_usr])),'Negative clearSemDia value: {}'.format(clearSemDia_usr)
            
            for idx, surface in enumerate (self.surfaceData):
                surface.clearSemDia = clearSemDia_usr[idx]
                self.surfaceData[idx]=surface
        else:
            print('Unmatch clearSemDia and surface dimension. # Surfaces: {}, # Semi Diameters: {}'.format(len(self.surfaceData),len(clearSemDia_usr)))
            
    def plot(self):
        fig, ax = plt.subplots()
        print(self)
        plot_system(self,fig=fig,ax=ax)
        
    def get_varValues(self):
        x0=[]
        for surface in self.surfaceData:
            for prop in surface.varProp:
                x0.append(getattr(surface, prop))
        return x0
    
    def set_varValues(self, x0):
        varValues = self.get_varValues()
        assert len(x0)==len(varValues), 'x0 does not match the number of variables. # of variables:{}, x0 length:{}'.format(len(varValues),len(x0))
        assert all([isinstance(element,(int,float)) for element in x0]), 'Wrong type of data in x0: {}'.format(x0)
        
        idxX0 = int(0)
        for surface in self.surfaceData:
            for prop in surface.varProp:
                setattr(surface, prop,x0[idxX0])
                idxX0 += int(1)

    def mark_variable(self, varName=None, surfIndex=-2):
        if varName is not None :
            self.surfaceData[surfIndex].mark_variable(varName)

        
    def unmark_variable(self, varName=None, surfIndex=-2):
        if varName is not None :
            self.surfaceData[surfIndex].unmark_variable(varName)
        
if __name__=='__main__':
    #Create system
    syst1 = OpSysData()
    print('\nDefault Optical system')
    #syst1.print_report()
    print(syst1)

    syst1.change_surface(optElement.standard(thickness= 2, curvature=0,material=1.1),surfIndex=0)
    
    #optElement.standard(thickness= 2, curvature=1,material=2),surfIndex=0
    syst1.add_surface(optElement.standard(thickness= 2, curvature=1,material=2, clearSemDia=4))
    
    #optElement.paraxial(thickness=10, focalLength=5, clearSemDia=12)
    syst1.add_surface(optElement.paraxial(thickness=10, focalLength=5, clearSemDia=12))
    syst1.mark_variable('focalLength')
    
    #syst1.add_surface(10,1/2.0,2,varProp='010')
    syst1.add_surface(optElement.standard(thickness= 10, curvature=1/2,material=1.45))
    syst1.mark_variable(surfIndex=3,varName='curvature')
    
    #syst1.add_surface(3,1/3.0,'N-BK7',varProp='010')
    syst1.add_surface(optElement.standard(thickness= 13, curvature=1/3,material='N-BK7'))
    syst1.mark_variable(varName='curvature')
    
    #syst1.add_surface(4,1/4.0,1,varProp='110')
    syst1.add_surface(optElement.standard(thickness= 8, curvature=1/4,material=1,clearSemDia=5))
    syst1.mark_variable('thickness',surfIndex=5)
    syst1.mark_variable(surfIndex=5,varName='curvature')
    
    syst1.add_surface(optElement.grin(curvature=-1/4,clearSemDia=7))
    
    
    #syst1.add_surface(4,1/4.0,1.5,surfPara=['test'],varProp='100')
    syst1.add_surface(optElement.standard(thickness= 10, curvature=1/100,material=1.1,clearSemDia=4))
    syst1.mark_variable('thickness')
    
    syst1.change_surface(optElement.standard(thickness= 3.14, curvature=1/100,material=1.1,clearSemDia=5),surfIndex=0)


    syst1.add_surface(optElement.aspheric(thickness= 10, curvature=1/100,material=1.1,clearSemDia=4))
    syst1.add_surface(optElement.aspheric(thickness= 10, curvature=-1/100,material=1.1,clearSemDia=4))
    syst1.mark_variable('thickness')
    

    print('\nOptical system data')
    #syst1.print_report()
    print(syst1)
    
    #Test changeSurface
    #syst1.change_surface(9,1/9,4,surfIndex=3)
    syst1.change_surface(optElement.standard(thickness= 9, curvature=1/9,material=4),surfIndex=3)
    print('\nSurface 3 changed')
    print(syst1)
    
    #Test deleteSurface
    syst1.delete_surface(surfIndex=3)
    print('\nSurface 3 deleted')
    print(syst1)
    #Plot system
    #syst1.change_clearSemDia(clearSemDia_usr=[4])
    syst1.plot()
    
    #Test invert surface
    print('---Invalid data type---')
    try:
        syst1.invert_surface_order(1,5)
    except Exception as e: print(e)
    syst1.invert_surface_order(1,4)
    print('\nSurfaces 1 to 4 inverted')
    print(syst1)
    #Plot system
    syst1.plot()
    syst1.mark_variable('material')
    
    
    syst1.change_clearSemDia([3])
    syst1.plot()
    # syst1.plot(clearSemDia_usr=[3,3,2,1.5,2,2,3])
    
    # #Test invalid clearSemDia
    # print('---Invalid data type---')
    # try:
    #     syst1.plot(clearSemDia_usr=[3,'3',2,1.5,2,2,3])
    # except Exception as e: print(e)
    
    # print('---Negative value---')
    # try:
    #     syst1.plot(clearSemDia_usr=[3,3,2,-1.5,2,2,3])
    # except Exception as e: print(e)
        
    # print('---Unmatch dimensions---')
    # try:
    #     syst1.plot(clearSemDia_usr=[3,5,3,2,4,1.5,2,2,3])
    # except Exception as e: print(e)
    
    #Test get values for optimization
    print(syst1)
    x0 = syst1.get_varValues()
    print(x0)
    syst1.set_varValues([1.111,2.222,3.333,4.444,7.777,1.999,8.88])
    print(syst1)
    syst1.plot()
    

    
