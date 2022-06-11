# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 20:43:53 2020
@author: David Vasquez
Classes RaySource, PointSource, InfinitySource
"""
import numpy as np

class RaySource:
    def __init__(self, XYZ=None,LMN=None,wvln=None): 
        '''
        Arguments
        XYZ : list [3xdouble], relative position in mm
        LMN : list [3xdouble], cosine directors in radians
        wvln: double, wavelength in nm
        
        Static method:
        calc_direcCos: Calculate direction cosine from a vector 
        ''' 
        self.RayList = [] 
        if XYZ != None or LMN != None or wvln != None:
            self.new_ray(XYZ,LMN,wvln) 

    def new_ray(self,XYZ,LMN,wvln):
        self.check_arguments(XYZ=XYZ, LMN=LMN, wvln=wvln)       
        self.RayList.append([XYZ,LMN,wvln])
    
    def change_XYZ(self,XYZ,Index):
        self.check_arguments(XYZ=XYZ)
        self.RayList[Index][0] = XYZ 
    
    def change_LMN(self, LMN, Index):
        self.check_arguments(LMN=LMN)
        self.RayList[Index][1] = LMN
        
    def change_wnlg(self,wvln,Index):
        self.check_arguments(wvln=wvln)
        self.RayList[Index][2] = wvln 
        
    def check_arguments(self,XYZ='nan',LMN='nan',wvln='nan'):
        if XYZ !='nan':
            assert (all([isinstance(q,(int,float)) for q in XYZ]) 
                    and len(XYZ) == 3), 'Invalid XYZ vector'
        if LMN !='nan':
            assert (all([isinstance(q,(int,float)) for q in LMN]) 
                    and len(LMN) == 3
                    and np.isclose(np.sum(np.power(LMN,2)),1)),'Invalid direction cosines'
        if wvln !='nan':
            assert isinstance(wvln,(int,float)),'Invalid wavelength'
            
    def print_report(self):
        headers = ['#','XPos','YPos','ZPos','XCosDir','YCosDir','ZCosDir','Wavelength']
        counter = 0
        print("\nRAY LIST")
        print("Class:" + self.__class__.__name__)
        print("{: >5} {: >12} {: >12} {: >12} {: >12} {: >12} {: >12} {: >12} ".format(*headers))
        for ray in self.RayList:
            print("{: 5d} ".format(counter)+
                  "{:12f} {:12f} {:12f} ".format(*ray[0])+
                  "{:12f} {:12f} {:12f} ".format(*ray[1])+
                  "{:12f} ".format( ray[2]))
            counter += 1

    
    @staticmethod        
    def calc_direcCos(Vector):
        assert (all([isinstance(q,(int,float)) for q in Vector]) 
                and len(Vector) == 3),'Invalid vector'
        #Calculate direction cosines
        norm     = np.sqrt(np.sum(np.power(Vector,2))) 
        cosDirX  = Vector[0]/norm
        cosDirY  = Vector[1]/norm
        cosDirZ  = Vector[2]/norm 
        
        return [cosDirX,cosDirY,cosDirZ]
        

class PointSource(RaySource):
    '''
    Inherit from RaySource. 
    A PointSource is a ray source with fixed position [XYZ] and wavelength (wvln).
    
    Attributes:
        Position  : list [3xdouble]
        Wavelength: double
        RayList   : list [rays] 
        
    Essential rays roll
    #           on axis                     off axis
    0           optical axis                chief ray
    1           y upper marginal ray        y upper coma ray          y=meridional
    2           y lower marginal ray        y lower coma ray
    3           x positive marginal ray     x positive coma ray       x=sagintal
    4           x negative marginal ray     x negative coma ray       positive/negative = pupil cartesian X sign 
    >= 5        screw ray                   screw ray
    '''
    
    def __init__(self, XYZ, wvln):                   #Vectors as entry
        # Initialize ray list
        self.Position   = XYZ
        self.Wavelength = wvln
        self.RayList    = []
        
        # Initialize essential rays
        for q in range (5):
            self.new_ray(self.Position,[0,0,1],self.Wavelength)
            
    # Override function to avoid XYZ and wnlg changes
    def change_XYZ(self,XYZ,Index):
        raise ReferenceError ('PointSource position [XYZ] not allowed to change') 
    def change_wnlg(self,wvln,Index):
        raise ReferenceError ('PointSource wavelength [wvln] not allowed to change')
  

class InfinitySource(RaySource):
    '''
    Inherit from RaySource. 
    A InfinitySource (source at infinity) is a ray source with parallel rays 
    (fixed direction cosine [LMN]) and fixed wavelength (wvln).
    
    Attributes:
        DirecCos  : list [3xdouble]
        Wavelength: double
        RayList   : list [rays] 
        
    Essential rays roll
    #           on axis                     off axis
    0           optical axis                chief ray
    1           y upper marginal ray        y upper coma ray          y=meridional
    2           y lower marginal ray        y lower coma ray
    3           x positive marginal ray     x positive coma ray       x=sagintal
    4           x negative marginal ray     x negative coma ray       positive/negative = pupil cartesian X sign 
    >= 5        screw ray                   screw ray
    '''
    
    def __init__(self, LMN, wvln):                   #Vectors as entry
        # Initialize ray list
        self.DirecCos = LMN
        self.Wavelength = wvln
        self.RayList = []
        
        # Initialize essential rays
        for q in range (5):
            self.new_ray([0,0,0],self.DirecCos,self.Wavelength)
            
    # Override function to avoid LMN and wnlg changes
    def change_LMN(self,XYZ,Index):
        raise ReferenceError ('InfinitySource direction [LMN] not allowed to change') 
    def change_wnlg(self,wvln,Index):
        raise ReferenceError ('InfinitySource wavelength [wvln] not allowed to change')
            
            
        
if __name__=='__main__':
    #Ray Source
    pto1=RaySource()
    #Add rays
    for Y in np.linspace(0,11,10):
        pto1.new_ray([0,Y,0],[0,1,0],365)  
    #Change wavelength
    pto1.change_wnlg(300,3)
    pto1.print_report()
    
    
    #Point Source
    pto2=PointSource([0,4,0],635)
    #change cosine director
    LMN=[0.48666426339228763, 0.3244428422615251, 0.8111071056538127]
    pto2.change_LMN(LMN,0)
    pto2.print_report()
    
    #InfinitySource
    pto3=InfinitySource([0,0,1],430)
    #Change position
    XYZ=[0,1,0]
    pto3.change_XYZ(XYZ,0)
    pto3.print_report()
    
    #Static function to calculate direction cosines from a vector
    LMN =RaySource.calc_direcCos([3,2,5])
    print(LMN)
    
    
    
