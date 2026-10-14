# -*- coding: utf-8 -*-
"""
Created on Wed Sep 23 22:57:23 2026

@author: vasquezpinzondavid
"""

try: import JenTrace
except ModuleNotFoundError: 
    import os, sys
    sys.path.insert(0,os.path.dirname(os.getcwd()))
    
from JenTrace.catalog import OPT_GLASS


class OpSurface:
    varProp     = set()
    validOptimizationParam = []
    
    def check_argument(self,x='nan',varName=None):
        if x!='nan':
            assert isinstance(x,(int,float)), varName + ' must be either int or float.'
        if varName=='clearSemDia':
            assert x>0, varName + '({}) should be higher than zero.'.format(x)
        
    def check_material(self,material='nan'):
        if material!='nan':
            assert isinstance(material,(int,float,str)), 'refraction index [n] must be either int,float or str'
            if isinstance (material,str):
                assert material in OPT_GLASS, '{} is not found in the optical glass catalog. The available materials are: {} '.format(material, [k for k in OPT_GLASS.keys()])
                
    def mark_variable(self,varName):
        assert varName in self.validOptimizationParam , '{} is not a valid optimization parameter'.format(varName)
        if varName == 'material':
            assert not isinstance(getattr(self, varName),str) , 'material: {} can not be set as optimization parameter'.format(getattr(self, varName))
        self.varProp.add(varName)
        
    def unmark_variable(self,varName):
        assert varName in self.validOptimizationParam , '{} is not a valid optimization parameter'.format(varName)
        self.varProp.discard(varName)
        
class standard (OpSurface):

    validOptimizationParam = ['thickness','curvature','material']
    
    def __init__(self, thickness=10, curvature=10, material=1, clearSemDia=10): 
         
        self.check_argument(thickness,'thickness') 
        self.check_argument(curvature,'curvature') 
        self.check_material(material) 
        self.check_argument(clearSemDia,'clearSemDia')
        
        self.surfType  = 'Standard'
        self.thickness  = thickness 
        self.curvature = curvature    
        self.material  = material 
        self.clearSemDia = clearSemDia
        
        self.varProp     = set()
        
    def __repr__(self):
        return (
            f"standard("
            f"thickness={self.thickness}, "
            f"curvature={self.curvature}, "
            f"material={self.material}, "
            f"clearSemDia={self.clearSemDia})"
        )
        
class paraxial (OpSurface):

    validOptimizationParam = ['thickness','focalLength']
    
    def __init__(self, thickness=10, focalLength=10, material=1, clearSemDia=10): 
        
        assert focalLength != 0, 'focalLength must be different from 0'
        self.check_argument(thickness,'thickness') 
        self.check_argument(focalLength,'focalLength') 
        self.check_material(material) 
        self.check_argument(clearSemDia,'clearSemDia')
        
        self.surfType   = 'Paraxial'
        self.thickness  = thickness 
        self.curvature  = None
        self.focalLength = focalLength  
        self.material    = material 
        self.clearSemDia = clearSemDia
        
        self.varProp     = set()
        
    def __repr__(self):
        return (
            f"paraxial("
            f"thickness={self.thickness}, "
            f"focalLength={self.focalLength}, "
            f"material={self.material}, "
            f"clearSemDia={self.clearSemDia})"
        )
        
class grin (OpSurface):
    '''
    Radial-Gradient-Index Media
    
    n(r)=n0(1-((A*(r**2))/2))
    
    n0 = B + C/(lam**2)
    
    A = (K0 + K1/(lam**2) + K2/(lam**4))**2
    
    '''
    
    validOptimizationParam = ['thickness','curvature',
                              'B', 'C', 'K0', 'K1', 'K2']
    
    def __init__(self, thickness=10, curvature=0,clearSemDia=10 , B=1, C=0, K0=1,K1=0,K2=0): 
        
        self.check_argument(thickness,'thickness') 
        self.check_argument(curvature,'curvature')
        self.check_argument(B,'B')
        self.check_argument(C,'C')
        self.check_argument(K0,'K0')
        self.check_argument(K1,'K1')
        self.check_argument(K2,'K2')
        self.check_argument(clearSemDia,'clearSemDia')
        
        self.surfType  = 'grin'
        self.thickness  = thickness 
        self.curvature = curvature
        self.material  = None
        self.clearSemDia = clearSemDia
        self.B  = B
        self.C  = C
        self.K0 = K0
        self.K1 = K1
        self.K2 = K2
        
        self.varProp     = set()
        
    def __repr__(self):
        return (
            f"grin("
            f"thickness={self.thickness}, "
            f"curvature={self.curvature}, "
            f"clearSemDia={self.clearSemDia}, "
            f"B={self.B}, "
            f"C={self.C}, "
            f"K0={self.K0}, "
            f"K1={self.K1}, "
            f"K2={self.K2})"
        )
        
class aspheric (OpSurface):
    '''
    Rotationally symmetric aspheric surface
    
    z = (c*(r**2))/(1+sqrt(1-(1+K)(c**2)(r**2))) + (A4)(r**4) + (A6)(r**6) + (A8)(r**8) + ...
    
    '''

    validOptimizationParam = ['thickness','curvature','material',
                              'K','A4','A6','A8','A10','A12','A14']
    
    def __init__(self, thickness=10, curvature=10, material=1, clearSemDia=10 ,K=0, A4=0, A6=0, A8=0, A10=0, A12=0, A14=0): 
        

        self.check_argument(thickness,'thickness') 
        self.check_argument(curvature,'curvature')
        self.check_material(material) 
        self.check_argument(clearSemDia,'clearSemDia')
        self.check_argument(K,'K')
        self.check_argument(A4,'A4')
        self.check_argument(A6,'A6')
        self.check_argument(A8,'A8')
        self.check_argument(A10,'A10')
        self.check_argument(A12,'A12')
        self.check_argument(A14,'A14')
        
        self.surfType  = 'aspheric'
        self.thickness  = thickness 
        self.curvature = curvature
        self.material  = material
        self.clearSemDia = clearSemDia
        self.K   = K
        self.A4  = A4
        self.A6  = A6
        self.A8  = A8
        self.A10 = A10
        self.A12 = A12
        self.A14 = A14
        
        self.varProp     = set()
        
    def __repr__(self):
        return (
            f"aspheric("
            f"thickness={self.thickness}, "
            f"curvature={self.curvature}, "
            f"material={self.material}, "
            f"clearSemDia={self.clearSemDia}, "
            f"K={self.K}, "
            f"A4={self.A4}, "
            f"A6={self.A6}, "
            f"A8={self.A8}, "
            f"A10={self.A10}, "
            f"A12={self.A12}, "
            f"A14={self.A14})"
        )
        
if __name__=='__main__':

    q = standard(3,1/3.0,'N-BK7')
    try:
        q.mark_variable('material')
    except Exception as excep: 
        print(excep)
    
    w = standard(4,1/4.0,1)
    e = standard(4,1/4.0,1.5,15)
    e.mark_variable('thickness')
    e.mark_variable('curvature')
    e.unmark_variable('thickness')
    print(e.varProp)
    
    q = paraxial(10,1)
    print(q.surfType)
    
    r = grin()
    print(r)
    
    t = aspheric()
    print(t)
    