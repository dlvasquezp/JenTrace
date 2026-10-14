# -*- coding: utf-8 -*-
"""
Created on Tue Mar 31 20:33:47 2020
@author: David Vasquez
function:   -plot_system
            -plot_rayTrace
"""
import math
import matplotlib.patches
import matplotlib.pyplot as plt
from numpy import pi
#import JenTrace.opt_sfc as optElement

def mechanical_SemDia(z,surface):
    if surface.__class__.__name__ in ["standard","grin","aspheric"]:   
        if surface.curvature==0:
            return [z,surface.clearSemDia]            
        else:
            radius = 1.0/surface.curvature 
            if abs(+surface.clearSemDia/radius) > 1:    #Check size
                return [z+radius, +abs(radius)]
            else:
                #r**2=apRadius**2+b**2, zSag = (-)r+b or [+]r-b
                b = math.sqrt(radius**2-surface.clearSemDia**2)
                if radius < 0:                             #Decide sag direction
                    return [z+radius+b,+surface.clearSemDia]
                else:
                    return [z+radius-b,+surface.clearSemDia]
    else:
        return [z,surface.clearSemDia]   

def plot_system(optObj, fig=[], ax=[]):
    '''
    plot_system accets either OpSysData or OpDesign
    '''   
    # Check optObj type
    if optObj.__class__.__name__=='OpSysData':
        optSystem = optObj
        apRadius  = 1.0
        apIndex   = -1
    else:
        if optObj.__class__.__name__=='OpDesign':
            optSystem = optObj.optSys
            apRadius  = optObj.aprRad
            apIndex   = optObj.aprInd
        else:
            raise TypeError ('plot_system accets either OpSysData or OpDesign')
    
    #Check figure and axis instance
    assert fig.__class__.__name__=='Figure','Invalid figure [fig] instance'
    assert ax.__class__.__name__ =='Axes','Invalid axes [ax] instance'
    
    
    #fig, ax = plt.subplots()  
    
    #number of surfaces
    surf_len = len(optSystem.surfaceData)
    
    #list the Z coordinates
    z0=[]
    z0.append(0)
    for idx, surface in enumerate(optSystem.surfaceData):
        if idx<(surf_len-1):
            z0.append(z0[-1] + surface.thickness)
            
    # create arc patches and lines 
    surfPatch = []
    lines2D   = []
    arrowPatch= []
    
    # Optical axis line
    lines2D.append(matplotlib.lines.Line2D((z0[0],z0[-1]),(0,0)
                                           ,c='black',lw=1.0,ls='--'))
    
    for idx, surface in enumerate(optSystem.surfaceData):
        if idx == 0 or idx == (surf_len-1):
            lines2D.append(matplotlib.lines.Line2D((z0[idx],z0[idx])
                                                  ,(+surface.clearSemDia,-surface.clearSemDia)
                                                  ,c='black',lw=1.0,ls='-'))
            
        elif surface.__class__.__name__ in ["standard","grin","aspheric"]:
            if surface.curvature==0:
                lines2D.append(matplotlib.lines.Line2D((z0[idx],z0[idx])
                                                      ,(+surface.clearSemDia,-surface.clearSemDia)
                                                      ,c='black',lw=1.0,ls='-'))
            else:    
                radius = (1.0/surface.curvature) 
                zR = z0[idx] + radius 
                #radius = 1.0/optSystem.surfaceData[q][1] 
                diam = abs(2*radius)
                if abs(+surface.clearSemDia/radius) > 1:
                    theta = 90.0
                else:
                    theta = math.asin(abs(+surface.clearSemDia/radius))*(180/pi)
            
                if radius < 0:
                    angleRot = 0
                else:
                    angleRot = 180
                    
                surfPatch.append(matplotlib.patches.Arc((zR,0)
                                                        ,diam
                                                        ,diam
                                                        ,angle=angleRot
                                                        ,theta1=-theta
                                                        ,theta2=+theta
                                                        ,color='black',lw=1.0,ls='-'))
    
        elif surface.__class__.__name__=="paraxial":
            arrowPatch.append(matplotlib.patches.FancyArrowPatch((z0[idx], +surface.clearSemDia) 
                                                                 ,(z0[idx], -surface.clearSemDia) 
                                                                 ,arrowstyle='<->' 
                                                                 ,mutation_scale=20 
                                                                 ,color='black',lw=1.0,ls='-'))    
        
        if surface.material!=1.0 and idx < (surf_len-1):
            [x0,y0]=mechanical_SemDia(z0[idx],surface) 
            [x1,y1]=mechanical_SemDia(z0[idx+1],optSystem.surfaceData[idx+1]) 
            #print(surface.material!=1.0,surface.material,[x0,y0],[x1,y1])
            if y0>y1:
                lines2D.append(matplotlib.lines.Line2D((+x0,+x1,+x1)
                                                      ,(+y0,+y0,+y1)
                                                      ,c='black',lw=1.0,ls='-'))
                
                lines2D.append(matplotlib.lines.Line2D((+x0,+x1,+x1)
                                                      ,(-y0,-y0,-y1)
                                                      ,c='black',lw=1.0,ls='-')) 
            else:
                lines2D.append(matplotlib.lines.Line2D((+x0,+x0,+x1)
                                                      ,(+y0,+y1,+y1)
                                                      ,c='black',lw=1.0,ls='-'))
                
                lines2D.append(matplotlib.lines.Line2D((+x0,+x0,+x1)
                                                      ,(-y0,-y1,-y1)
                                                      ,c='black',lw=1.0,ls='-'))    
            
            # if surface.curvature==0:
            # else:
            #     radius = 1.0/surface.curvature 
            #     if abs(+surface.clearSemDia/radius) > 1:    #Check size
            #         yP.append([z0[q]+radius, +abs(radius)])
            #     else:
            #         #r**2=apRadius**2+b**2, zSag = (-)r+b or [+]r-b
            #         b = math.sqrt(radius**2-lensRadius[q]**2)
            #         if radius < 0:                             #Decide sag direction
            #             yP.append([z0[q]+radius+b,+lensRadius[q]])
            #         else:
            #             yP.append([z0[q]+radius-b,+lensRadius[q]])
            
    
    # for q in surfPatch:
    #     fig.gca().add_patch(q)
    # for w in lines2D:
    #     ax.add_line(w)
    # for e in arrowPatch:
    #     ax.add_patch(e)
        
    # xlim =[z0[0]-1,z0[-1]+1]
    # yPoints=[q.clearSemDia for q in optSystem.surfaceData]
    # maxY = max(yPoints)+1
    # ylim =[-maxY,+maxY]
    # #ylim =[-11,+11]
    
    # ax.set_xlim(xlim)
    # ax.set_ylim(ylim)
    # plt.show()
    
    # #Check if the lens radius are complete
    # lensRadius =[]
    # for q in range(surf_len):
    #     lensRadius.append(optSystem.surfaceData[q][4])
    
    
    # #list the surface center
    # zR=[]
    # for q in range(surf_len):
    #     if optSystem.surfaceData[q][1] == 0:
    #         zR.append(z0[q])
    #     else:
    #         zR.append(z0[q]+(1.0/optSystem.surfaceData[q][1]))
    
    # #List off-axis points // Only the upper point is plotted, the 2nd is mirrorred
    # yP=[]
    # for q in range(surf_len):
    #     if q == 0 or q==(surf_len-1):                          #Skipt the first and last surface
    #         yP.append([z0[q],+lensRadius[q]])
    #     else:
    #         if optSystem.surfaceData[q][1] == 0:               #Check curvature different from 0
    #             yP.append([z0[q],+lensRadius[q]])
    #         else:
    #             radius = 1.0/optSystem.surfaceData[q][1] 
    #             if abs(+lensRadius[q]/radius) > 1:                   #Check size
    #                 yP.append([z0[q]+radius, +abs(radius)])
    #             else:
    #                 #r**2=apRadius**2+b**2, zSag = (-)r+b or [+]r-b
    #                 b = math.sqrt(radius**2-lensRadius[q]**2)
    #                 if radius < 0:                             #Decide sag direction
    #                     yP.append([z0[q]+radius+b,+lensRadius[q]])
    #                 else:
    #                     yP.append([z0[q]+radius-b,+lensRadius[q]])

    ''' 
    List additional connection points // 'NA': Not applicable
                                      // [0,0]: Straight line
                                      // [z,y]: connection point coordinate
    '''
    # conP=[]
    # for q in range(surf_len):
    #     if q == (surf_len-1):                    #Skip the last surface
    #         conP.append('NA')
    #     else:
    #         if optSystem.surfaceData[q][2] == 1: #Refraction index equal to air
    #             conP.append('NA')
    #         else:
    #             if yP[q][1] == yP[q+1][1]:       #off axis points comparison
    #                 conP.append([0,0])
    #             else:
    #                 if yP[q][1] < yP[q+1][1]:    # Size comparizon
    #                     conP.append([yP[q][0],yP[q+1][1]])
    #                 else:
    #                     conP.append([yP[q+1][0],yP[q][1]])
   
    # # create arc patches and lines 
    # surfPatch= []
    # lines2D  = []
    
    # # Optical axis line
    # lines2D.append(matplotlib.lines.Line2D((z0[0],z0[-1]),(0,0)
    #                                        ,c='black',lw=1.0,ls='--'))
    
    # # Surfaces
    # for q in range(surf_len):
    #     #Plot first surface, last surface and curvature equal zero as planes
    #     if  q == 0 or q==(surf_len-1) or optSystem.surfaceData[q][1] == 0:
    #         lines2D.append(matplotlib.lines.Line2D((yP[q][0],yP[q][0])
    #                                               ,(yP[q][1],-yP[q][1])
    #                                               ,c='black',lw=1.0,ls='-'))
    #     else:
    #         radius = 1.0/optSystem.surfaceData[q][1] 
    #         diam = abs(2*radius)
    #         if abs(+lensRadius[q]/radius) > 1:
    #             theta = 90.0
    #         else:
    #             theta = math.asin(abs(+lensRadius[q]/radius))*(180/pi)
            
    #         if radius < 0:
    #             angleRot = 0
    #         else:
    #             angleRot = 180
                    
    #         surfPatch.append(matplotlib.patches.Arc((zR[q],0)
    #                                                     ,diam
    #                                                     ,diam
    #                                                     ,angle=angleRot
    #                                                     ,theta1=-theta
    #                                                     ,theta2=+theta
    #                                                     ,color='black',lw=1.0,ls='-'))
    #     #Connection lines
    #     if conP[q]!= 'NA' and q!= (surf_len-1):
    #         if conP[q]==[0,0]:
    #             lines2D.append(matplotlib.lines.Line2D((yP[q][0],yP[q+1][0])
    #                                                   ,(+yP[q][1],+yP[q+1][1])
    #                                                   ,c='black',lw=1.0,ls='-'))
                
    #             lines2D.append(matplotlib.lines.Line2D((yP[q][0],yP[q+1][0])
    #                                                   ,(-yP[q][1],-yP[q+1][1])
    #                                                   ,c='black',lw=1.0,ls='-')) 
    #         else:
    #             lines2D.append(matplotlib.lines.Line2D((yP[q][0],conP[q][0],yP[q+1][0])
    #                                                   ,(+yP[q][1],+conP[q][1],+yP[q+1][1])
    #                                                   ,c='black',lw=1.0,ls='-'))
                
    #             lines2D.append(matplotlib.lines.Line2D((yP[q][0],conP[q][0],yP[q+1][0])
    #                                                   ,(-yP[q][1],-conP[q][1],-yP[q+1][1])
    #                                                   ,c='black',lw=1.0,ls='-'))
    
    # Plot aperture stop
    HIGH=1.5
    if apIndex != -1:
        """ # The aperture is draw on the curvature of the surface and not simply on the reference plane
        lines2D.append(matplotlib.lines.Line2D((z0[apIndex],z0[apIndex]),(apRadius,apRadius*HIGH)
                                               ,c='black',lw=2.0,ls='-'))
        lines2D.append(matplotlib.lines.Line2D((z0[apIndex],z0[apIndex]),(-apRadius,-apRadius*HIGH)
                                               ,c='black',lw=2.0,ls='-'))
        """
        apX=[]
        if optSystem.surfaceData[apIndex][1] == 0:     #Check curvature different from 0
                apX.append(z0[apIndex])
        else:
            radius = 1/optSystem.surfaceData[apIndex][1]
            #r**2=apRadius**2+b**2, zSag = (-)r+b or [+]r-b
            b = math.sqrt(radius**2-apRadius**2)
            if radius < 0:                             #Decide sag direction
                apX.append(z0[apIndex]+radius+b)
            else:
                apX.append(z0[apIndex]+radius-b)
        
        
        lines2D.append(matplotlib.lines.Line2D((apX[0],apX[0]),(+apRadius,+apRadius*HIGH)
                                               ,c='black',lw=2.0,ls='-'))
        lines2D.append(matplotlib.lines.Line2D((apX[0],apX[0]),(-apRadius,-apRadius*HIGH)
                                               ,c='black',lw=2.0,ls='-'))
    
    
    #Append arcs and lines to the figigure (fig) and axes (ax)
    for q in surfPatch:
        fig.gca().add_patch(q)
    #fig, ax = plt.subplots()    
    for w in lines2D:
        ax.add_line(w)
    for e in arrowPatch:
         ax.add_patch(e)
    #plt.show()
    
    #Plot limits
    xlim =[z0[0]-1,z0[-1]+1]
    #yPoints=[q[1] for q in yP]
    yPoints=[q.clearSemDia for q in optSystem.surfaceData]
    maxY = max(yPoints)+1
    ylim =[-maxY,+maxY]
    
    ax.set_xlim(xlim)
    ax.set_ylim(ylim)
    
    return fig,ax

def plot_rayTrace(rayTrace, fig=[], ax=[], color='b'):
    #Check instances
    assert rayTrace.__class__.__name__=='ndarray','Invalid rayTrace' 
    assert fig.__class__.__name__=='Figure'      ,'Invalid figure [fig] instance'
    assert ax.__class__.__name__ =='Axes' ,'Invalid axes [ax] instance'
    #number of surfaces
    dim = rayTrace.shape
    lines2D  = []
    
    for q in range(dim[0]):
        #list the Z coordinates
        z0=[]
        for w in range(dim[2]):
            if w == 0:
                z0.append(0)
            else:
                z0.append(z0[w-1] + rayTrace[q][0][w-1])
    
        for w in range(1,dim[2]):
            lines2D.append(matplotlib.lines.Line2D((z0[w-1]+ rayTrace[q][7][w-1]
                                                   ,z0[w]  + rayTrace[q][7][w])
                                                   ,(rayTrace[q][9][w-1],rayTrace[q][9][w])
                                                   ,color=color))    
    
    for w in lines2D:
        ax.add_line(w)
    
    #Check Plot limits
    xlim = ax.get_xlim()
    ylim = ax.get_ylim()
    newXlim = [min(z0),max(z0)]
    newYlim = [rayTrace[:,9,:].min(),rayTrace[:,9,:].max()]
    
    if newXlim[0] < xlim[0]: xlim = [newXlim[0]-1,   xlim[1]  ]
    if newXlim[1] > xlim[1]: xlim = [   xlim[0]  ,newXlim[1]+1]
    if newYlim[0] < ylim[0]: ylim = [newYlim[0]-1,   ylim[1]  ]
    if newYlim[1] > ylim[1]: ylim = [   ylim[0]  ,newYlim[1]+1]

    ax.set_xlim(xlim)
    ax.set_ylim(ylim)    
    
    return fig,ax

if __name__ == '__main__':
    from ray_src import PointSource
    from opt_sys import OpSysData
    #from opt_dsg import OpDesign
    #Optical system
    syst1 = OpSysData()
    syst1.add_surface(2,0,1)
    syst1.add_surface(2,+0.2,1.5)
    syst1.add_surface(1,-0.2,1)
    syst1.add_surface(4,-0.2,1.1)
    #syst1.add_surface(0,0,1)
    syst1.add_surface(3,-0.3,1)
    syst1.add_surface(4,+0.4,1.8)
    print(syst1)
    clearSemDia=[2,1,3,2,2,1.5,1.5,1]
    syst1.plot(clearSemDia)
    #Ray source
    #pto1  = PointSource([0,0.5,0],635)
    #Plot Optical system
    #syst1.invert_surface_order(0,5)#;clearSemDia=[1,1,1,1,2,1]
    fig, ax = plt.subplots()
    plot_system(syst1,fig=fig,ax=ax)
    plt.show()
    #Plot Optical design
    #design1  = OpDesign(pto1,syst1,aprRad=1.5,aprInd=2)
    #fig, ax = plt.subplots()
    #fig, ax = plot_system(design1,clearSemDia=clearSemDia,fig=fig,ax=ax,show=True)
    #fig, ax = plot_rayTrace(design1.raySrcTrace,fig=fig,ax=ax,color='b')
    #fig, ax = plot_rayTrace(design1.dsgPtoTrace,fig=fig,ax=ax,color='g')
    #fig, ax = plot_rayTrace(design1.dsgInfTrace,fig=fig,ax=ax,color='r')
    
 