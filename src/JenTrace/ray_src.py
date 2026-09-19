# -*- coding: utf-8 -*-
"""
Created on Fri Mar 20 20:43:53 2020
@author: David Vasquez
Classes RaySource, PointSource, InfinitySource

Version 3.0: each ray (in RayList and HauptRay) is now stored as a
dictionary with keys 'XYZ', 'LMN' and 'wvln' instead of a plain list
[XYZ, LMN, wvln]. This makes the code more readable and less error-prone,
since ray attributes are accessed by name instead of by position.
"""
import numpy as np
import matplotlib.pyplot as plt


class RaySource:
    """
    Base class representing a generic optical ray source.

    A ray source stores a collection of rays. Each ray is represented as a
    dictionary with the following keys:
        'XYZ'  : list [3xdouble], relative position in mm
        'LMN'  : list [3xdouble], direction cosines (dimensionless)
        'wvln' : double, wavelength in nm

    Attributes:
        HauptRay : dict or None, "main rays" (chief/marginal rays) keyed by
                   name, each value being a ray dictionary as described above
        RayList  : list of ray dictionaries
    """

    def __init__(self, XYZ=None, LMN=None, wvln=None):
        """
        Arguments:
            XYZ  : list [3xdouble], relative position in mm
            LMN  : list [3xdouble], direction cosines (dimensionless)
            wvln : double, wavelength in nm

        If any argument is provided, a new ray is immediately created and
        appended to RayList via new_ray().

        Static method:
            calc_direcCos: Calculate direction cosines from a vector
        """
        self.HauptRay = None
        self.RayList = []
        if XYZ != None or LMN != None or wvln != None:
            self.new_ray(XYZ, LMN, wvln)

    def new_ray(self, XYZ, LMN, wvln):
        """
        Create a new ray and append it to RayList.

        Arguments:
            XYZ  : list [3xdouble], relative position in mm
            LMN  : list [3xdouble], direction cosines
            wvln : double, wavelength in nm
        """
        self.check_arguments(XYZ=XYZ, LMN=LMN, wvln=wvln)
        self.RayList.append({'XYZ': XYZ, 'LMN': LMN, 'wvln': wvln})

    def change_XYZ(self, XYZ, Index):
        """
        Change the position (XYZ) of the ray at position 'Index' in RayList.

        Arguments:
            XYZ   : list [3xdouble], new relative position in mm
            Index : int, index of the ray in RayList to modify
        """
        self.check_arguments(XYZ=XYZ)
        self.RayList[Index]['XYZ'] = XYZ

    def change_LMN(self, LMN, Index):
        """
        Change the direction cosines (LMN) of the ray at position 'Index'
        in RayList.

        Arguments:
            LMN   : list [3xdouble], new direction cosines
            Index : int, index of the ray in RayList to modify
        """
        self.check_arguments(LMN=LMN)
        self.RayList[Index]['LMN'] = LMN

    def change_wnlg(self, wvln, Index):
        """
        Change the wavelength (wvln) of the ray at position 'Index' in
        RayList.

        Arguments:
            wvln  : double, new wavelength in nm
            Index : int, index of the ray in RayList to modify
        """
        self.check_arguments(wvln=wvln)
        self.RayList[Index]['wvln'] = wvln

    def check_arguments(self, XYZ='nan', LMN='nan', wvln='nan'):
        """
        Validate ray arguments before they are stored.

        Arguments:
            XYZ  : list [3xdouble] or 'nan' (skip check)
            LMN  : list [3xdouble] or 'nan' (skip check)
            wvln : double or 'nan' (skip check)

        Raises:
            AssertionError if any provided argument is invalid.
        """
        if XYZ != 'nan':
            assert (all([isinstance(q, (int, float)) for q in XYZ])
                    and len(XYZ) == 3), 'Invalid XYZ vector'
        if LMN != 'nan':
            assert (all([isinstance(q, (int, float)) for q in LMN])
                    and len(LMN) == 3
                    and np.isclose(np.sum(np.power(LMN, 2)), 1)), \
                'Invalid direction cosines. Type {}'.format(type(LMN))
        if wvln != 'nan':
            assert isinstance(wvln, (int, float)), 'Invalid wavelength'

    def __str__(self):
        """
        Print a formatted table of HauptRay (if present) and RayList,
        and return a short marker string.
        """
        print("\nRAY LIST")
        print("Class:" + self.__class__.__name__)

        if self.HauptRay != None:
            Headers = ['HauptRay', 'XPos', 'YPos', 'ZPos',
                       'XCosDir', 'YCosDir', 'ZCosDir', 'Wavelength']
            print("{: >14} {: >12} {: >12} {: >12} {: >12} {: >12} {: >12} {: >12} ".format(*Headers))
            for RayKey, Ray in self.HauptRay.items():
                print("{: >14} ".format(RayKey) +
                      "{:12f} {:12f} {:12f} ".format(*Ray['XYZ']) +
                      "{:12f} {:12f} {:12f} ".format(*Ray['LMN']) +
                      "{:12f} ".format(Ray['wvln']))
            print(' ')

        Headers = ['Ray#', 'XPos', 'YPos', 'ZPos',
                   'XCosDir', 'YCosDir', 'ZCosDir', 'Wavelength']
        print("{: >5} {: >12} {: >12} {: >12} {: >12} {: >12} {: >12} {: >12} ".format(*Headers))
        Counter = 0
        for Ray in self.RayList:
            print("{: 5d} ".format(Counter) +
                  "{:12f} {:12f} {:12f} ".format(*Ray['XYZ']) +
                  "{:12f} {:12f} {:12f} ".format(*Ray['LMN']) +
                  "{:12f} ".format(Ray['wvln']))
            Counter += 1
        return str('#####')

    @staticmethod
    def calc_direcCos(Vector):
        """
        Calculate direction cosines from a 3D vector.

        Arguments:
            Vector : list [3xdouble], arbitrary (non-normalized) vector

        Returns:
            list [3xdouble], the normalized direction cosines [L, M, N]
        """
        assert (all([isinstance(q, (int, float, np.int64)) for q in Vector])
                and len(Vector) == 3), 'Invalid vector'
        # Calculate direction cosines
        Norm = np.sqrt(np.sum(np.power(Vector, 2)))
        CosDirX = Vector[0] / Norm
        CosDirY = Vector[1] / Norm
        CosDirZ = Vector[2] / Norm

        return [CosDirX, CosDirY, CosDirZ]
    
    def plot_ray_source(self, pupilPos=10, Ax=None, showHauptRay=False):
        """
        Plot RayList from RaySource in a 3D matplotlib plot.

        Draws:
            - For each ray in self.RayList, a line starting at the ray's
              position ('XYZ') and extending along its direction cosines
              ('LMN') towards positive z, up to the plane z = pupilPos
        """
        if Ax is None:
            Fig = plt.figure()
            Ax = Fig.add_subplot(111, projection='3d')

        Ax.set_title( self.__class__.__name__ )
        Ax.view_init(elev=30, azim=145, roll=-120)
            
        for Ray in self.RayList:
            XRay, YRay, ZRay = Ray['XYZ']
            L, M, N = Ray['LMN']
            T = (pupilPos - ZRay) / N
            XEnd, YEnd, ZEnd = XRay + T * L, YRay + T * M, ZRay + T * N
            Ax.plot([XRay, XEnd], [YRay, YEnd], [ZRay, ZEnd], color='blue')
        
        if self.HauptRay != None and showHauptRay:
            for Ray in self.HauptRay:
                XRay, YRay, ZRay = self.HauptRay[Ray]['XYZ']
                L, M, N = self.HauptRay[Ray]['LMN']
                T = (pupilPos - ZRay) / N
                XEnd, YEnd, ZEnd = XRay + T * L, YRay + T * M, ZRay + T * N
                Ax.plot([XRay, XEnd], [YRay, YEnd], [ZRay, ZEnd], color='black', linestyle='--')
            
        Xlim, Ylim, Zlim = Ax.get_xlim3d(), Ax.get_ylim3d(), Ax.get_zlim3d()
        Ax.set_axis_off()
        #Draw plain axis lines crossing through the origin instead
        Ax.plot([-5,5], [0, 0], [0, 0], color='black', lw=1)
        Ax.plot([0, 0], [-5,5], [0, 0], color='black', lw=1)
        Ax.plot([0, 0], [0, 0],[0,5], color='black', lw=1)
        Ax.text(5, 0, 0, 'X')
        Ax.text(0, 5, 0, 'Y')
        Ax.text(0, 0, 5, 'Z')
        Ax.set_xlim3d(Xlim)
        Ax.set_ylim3d(Ylim)
        Ax.set_zlim3d(Zlim)
        Ax.set_aspect('equal')
        
        return Ax
    
    def plot_ray_pupil(self, pts, pupilRad, pupilCenterXY=(0,0), Ax=None):
        """
        Plot RayList from RaySource in a 3D matplotlib plot.

        Draws:
            - For each ray in self.RayList, a line starting at the ray's
              position ('XYZ') and extending along its direction cosines
              ('LMN') towards positive z, up to the plane z = pupilPos
        """
        if Ax is None:
            Fig = plt.figure()
            Ax = Fig.add_subplot()
            
        Ax.scatter(pts[:, 0], pts[:, 1], s=40)
        circle = plt.Circle(pupilCenterXY, pupilRad, fill=False)
        Ax.add_patch(circle)
        Ax.set_aspect('equal')   # important: otherwise the circle looks like an ellipse
        # Labels, title, and a nice grid
        Ax.set_xlabel('X [mm]')
        Ax.set_ylabel('Y [mm]')
        Ax.set_title('Pupil plane')
        Ax.grid(True, linestyle='--', linewidth=0.5, alpha=0.7)
        Ax.set_axisbelow(True) 
        
        return Ax
    
    def ray_pattern(self, pupilRad, pupilPos, wvln, pattern:str, noRings:int=10, show=False, clearRayList=True):
        """
		Generate a fan of rays that sample the pupil plane in a chosen pattern,
		and add them to self.RayList.The pupil is treated as a disk of radius 
        'pupilRad' lying in the plane Z = pupilPos.
        How each ray is aimed through those pupil points depends on the type of 
        source:
			- PointSource   : all rays start at the fixed source position
							(self.Position) and are aimed at each pupil
							point.
			- InfinitySource: all rays share the same fixed direction
							(self.DirecCos).

		Arguments:
			pupilRad : double, radius of the pupil in mm
			pupilPos : double, Z position of the pupil plane in mm
			wvln     : double, wavelength assigned to every generated ray,
					in nm
			pattern  : str, sampling pattern to use. One of:
				'hexapolar'  : rays arranged in a hexagonal grid of
								'noRings' concentric rings, cropped to the
								unit disk (a classic pupil-sampling pattern
								for aberration analysis)
				'sagittal'    : rays arranged along a single line on the X
								axis (Y=0) -- i.e. the sagittal plane fan
				'tangential' : rays arranged along a single line on the Y
								axis (X=0) -- i.e. the tangential plane fan
			noRings  : int, for 'hexapolar' this is the number of hexagonal
					rings. For 'sagital'/'tangential', this instead
					controls the number of points along the line, which
					is (2*noRings + 1). Default 10.
			show     : bool, if True, display a side-by-side figure with a
					3D view of the rays and a 2D scatter plot of the
					sampled pupil points with the pupil boundary circle.
					Default False.

		"""
        assert isinstance(pupilRad,(int,float)), 'Invalid data type, pupilRad type is {}, and must be either int or float.'.format(type(pupilRad))
        assert isinstance(pupilPos,(int,float)), 'Invalid data type, pupilPos type is {}, and must be either int or float.'.format(type(pupilPos))
        assert isinstance(wvln,(int,float)), 'Invalid data type, wvln type is {}, and must be either int or float.'.format(type(wvln))
        assert wvln > 0, 'wvln must be positive'
        patternList=['hexapolar', 'sagittal','tangential']
        assert pattern in patternList, '{} is not a valid method. Available methods are {}'.format(pattern,patternList)
        assert isinstance(self,(PointSource,InfinitySource)),'Ray source must be either "PointSOurce" or "InfinitySource"'
        
        if clearRayList:
            self.RayList=[]
            
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
        
        elif pattern == 'sagittal':
            x = np.linspace(-1,1,(noRings*2+1))
            y = np.zeros(len(x))
            centers=[[q,w] for q,w in zip (x, y)]  
            centers=np.array(centers)
            
        elif pattern == 'tangential':
            y = np.linspace(-1,1,(noRings*2+1))
            x = np.zeros(len(y))
            centers=[[q,w] for q,w in zip (x, y)]  
            centers=np.array(centers)
        else:
            raise ValueError(f"Unknown pattern '{pattern}'")
            
        mask = np.sqrt(centers[:, 0]**2 + centers[:, 1]**2) <= 1
        pts=centers[mask]
        
        #Resize pupil radious from 1 to design value
        pts = np.multiply(pts,pupilRad)
        
        raysPupXYZ = []
        for px, py in pts:
            raysPupXYZ.append(np.array([px,py,pupilPos])) 
        
        if isinstance(self,PointSource):
            XYZ= self.Position 

            raysImg2Pup = []
            for ray in raysPupXYZ:
                raysImg2Pup.append(ray-XYZ)
        
            for ray in raysImg2Pup:
                LMN     = self.calc_direcCos(ray)
                self.new_ray(XYZ,LMN,wvln)
        
        elif isinstance(self,InfinitySource): 
            LMN=self.DirecCos 

            raysImg2Pup = []
            for ray in raysPupXYZ:
                raysImg2Pup.append(list(ray-(np.multiply(LMN,pupilPos/LMN[2]))))
            
            for XYZ in raysImg2Pup:
                self.new_ray(XYZ,LMN,wvln)
                
                
        if show:
            Fig = plt.figure(figsize=(12, 6))
            
            # 3D figure
            Ax0 = Fig.add_subplot(1, 2, 1, projection='3d')
            Ax0 = self.plot_ray_source(pupilPos, Ax=Ax0)   # pass Ax0 in
            
            # # Pupil plane
            Ax1 = Fig.add_subplot(1, 2, 2)
            Ax1 = self.plot_ray_pupil(pts, pupilRad, Ax=Ax1)
            
            plt.tight_layout()
            plt.show()
    

class PointSource(RaySource):
    """
    Inherit from RaySource.
    A PointSource is a ray source with fixed position [XYZ] and wavelength
    (wvln).

    Attributes:
        Position   : list [3xdouble]
        Wavelength : double
        RayList    : list of ray dictionaries ({'XYZ', 'LMN', 'wvln'})
        HauptRay   : dict of "main rays" keyed by name, each value a ray
                     dictionary ({'XYZ', 'LMN', 'wvln'})

    Haupt Rays naming convention:
    'chief'         on axis: optical axis / off axis: chief ray
    'marginal_py'   y upper marginal ray / y upper coma ray
    'marginal_ny'   y lower marginal ray / y lower coma ray
    'marginal_px'   x positive marginal ray / x positive coma ray
    'marginal_nx'   x negative marginal ray / x negative coma ray
    """

    HauptRayNames = ['chief', 'marginal_py', 'marginal_ny',
                     'marginal_px', 'marginal_nx']

    def __init__(self, XYZ, wvln):
        """
        Arguments:
            XYZ  : list [3xdouble], fixed position of the source in mm
            wvln : double, fixed wavelength in nm
        """
        # Initialize ray list
        self.Position = XYZ
        self.Wavelength = wvln
        self.RayList = []

        # Initialize Haupt rays
        self.HauptRay = {}
        for RayName in self.HauptRayNames:
            self.HauptRay[RayName] = {'XYZ': list(self.Position),
                                       'LMN': [0, 0, 1],
                                       'wvln': self.Wavelength}

    def change_LMN(self, LMN, RayKey):
        """
        Override: modify the direction cosines (LMN) of a Haupt ray,
        since PointSource rays are identified by name (RayKey) rather
        than index.

        Arguments:
            LMN    : list [3xdouble], new direction cosines
            RayKey : str, name of the Haupt ray to modify (see
                     HauptRayNames)
        """
        self.check_arguments(LMN=LMN)
        self.HauptRay[RayKey]['LMN'] = LMN

    def change_XYZ(self, XYZ, Index):
        """
        Override: disabled, since a PointSource's position is fixed.

        Raises:
            ReferenceError always.
        """
        raise ReferenceError('PointSource position [XYZ] not allowed to change')

    def change_wnlg(self, wvln, Index):
        """
        Override: disabled, since a PointSource's wavelength is fixed.

        Raises:
            ReferenceError always.
        """
        raise ReferenceError('PointSource wavelength [wvln] not allowed to change')


class InfinitySource(RaySource):
    """
    Inherit from RaySource.
    An InfinitySource (source at infinity) is a ray source with parallel
    rays (fixed direction cosines [LMN]) and fixed wavelength (wvln).

    Attributes:
        DirecCos   : list [3xdouble]
        Wavelength : double
        RayList    : list of ray dictionaries ({'XYZ', 'LMN', 'wvln'})
        HauptRay   : dict of "main rays" keyed by name, each value a ray
                     dictionary ({'XYZ', 'LMN', 'wvln'})

    Haupt Ray naming convention:
    'chief'         on axis: optical axis / off axis: chief ray
    'marginal_py'   y upper marginal ray / y upper coma ray
    'marginal_ny'   y lower marginal ray / y lower coma ray
    'marginal_px'   x positive marginal ray / x positive coma ray
    'marginal_nx'   x negative marginal ray / x negative coma ray
    """

    HauptRayNames = ['chief', 'marginal_py', 'marginal_ny',
                     'marginal_px', 'marginal_nx']

    def __init__(self, LMN, wvln):
        """
        Arguments:
            LMN  : list [3xdouble], fixed direction cosines of the source
            wvln : double, fixed wavelength in nm
        """
        # Initialize ray list
        self.DirecCos = LMN
        self.Wavelength = wvln
        self.RayList = []

        # Initialize Haupt rays
        self.HauptRay = {}
        for RayName in self.HauptRayNames:
            self.HauptRay[RayName] = {'XYZ': [0, 0, 0],
                                       'LMN': list(self.DirecCos),
                                       'wvln': self.Wavelength}

    def change_XYZ(self, XYZ, RayKey):
        """
        Override: modify the position (XYZ) of a Haupt ray, since
        InfinitySource rays are identified by name (RayKey) rather than
        index.

        Arguments:
            XYZ    : list [3xdouble], new position
            RayKey : str, name of the Haupt ray to modify (see
                     HauptRayNames)
        """
        self.check_arguments(XYZ=XYZ)
        self.HauptRay[RayKey]['XYZ'] = XYZ

    def change_LMN(self, XYZ, Index):
        """
        Override: disabled, since an InfinitySource's direction is fixed.

        Raises:
            ReferenceError always.
        """
        raise ReferenceError('InfinitySource direction [LMN] not allowed to change')

    def change_wnlg(self, wvln, Index):
        """
        Override: disabled, since an InfinitySource's wavelength is fixed.

        Raises:
            ReferenceError always.
        """
        raise ReferenceError('InfinitySource wavelength [wvln] not allowed to change')


if __name__ == '__main__':
    # Ray Source
    Pto1 = RaySource()
    # Add rays
    for Y in np.linspace(0, 11, 10):
        Pto1.new_ray([0, Y, 0], [0, 1, 0], 365)
    # Change wavelength
    Pto1.change_wnlg(300, 3)
    print(Pto1)

    # Point Source
    Pto2 = PointSource([0, 4, 0], 635)
    # Change cosine director
    LMN = [0.48666426339228763, 0.3244428422615251, 0.8111071056538127]
    Pto2.change_LMN(LMN, 'chief')
    print(Pto2)

    # InfinitySource
    Pto3 = InfinitySource([0, 0, 1], 430)
    # Change position
    XYZ = [0, 1, 0]
    Pto3.change_XYZ(XYZ, 'marginal_py')
    print(Pto3)

    # Static function to calculate direction cosines from a vector
    LMN = RaySource.calc_direcCos([3, 2, 5])
    print(LMN)
    
    #Raypattern
    t1=RaySource.ray_pattern(Pto2,10, 10, 300,'hexapolar', noRings=10, show=True)
    print(Pto2)
    t1=RaySource.ray_pattern(Pto2,10, 10, 300,'sagittal', noRings=5, show=True)
    print(Pto2)
    t1=RaySource.ray_pattern(Pto2,5, 10, 300,'tangential', noRings=5, show=True)
    print(Pto2)
    t1=RaySource.ray_pattern(Pto3,5, 13, 300,'hexapolar', noRings=5, show=True)
    print(Pto3)

