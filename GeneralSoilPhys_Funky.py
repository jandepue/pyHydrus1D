# -*- coding: utf-8 -*-
"""
Created on Wed Oct 28 13:11:50 2015

@author: Administrator
"""

#==============================================================================
# General Soil Physics functions
#==============================================================================

import numpy
import pylab

#==============================================================================
# Textuur
#==============================================================================

def DefineTexClassPolygons():
    
    #==============================================================================
    # define polygon
    #==============================================================================
    
    ## Draw polygons: Sand / Clay coordinates
    poly_U = numpy.array([[0.0  ,100.0],
                          [65.0 ,35.0 ],
                          [10.0 ,35.0 ],
                          [0.0  ,45.0 ],
                          [0.0  ,100.0]])
    
    poly_E = numpy.array([[65.0 ,35.0 ],
                          [10.0 ,35.0 ],
                          [0.0  ,45.0 ],
                          [0.0  ,30.0 ],
                          [25.0 ,17.5 ],
                          [82.5 ,17.5 ],
                          [65.0 ,35.0 ]])
    
    poly_A = numpy.array([[0.0  ,30.0 ],
                          [15.0 ,22.5 ],
                          [15.0 ,0.0  ],
                          [0.0  ,0.0  ],
                          [0.0  ,30.0 ]])
    
    poly_L = numpy.array([[15.0 ,22.5 ],
                          [25.0 ,17.5 ],
                          [67.5 ,17.5 ],
                          [67.5 ,12.5 ],
                          [50.0 ,12.5 ],
                          [50.0 ,0.0  ],
                          [15.0 ,0.0  ],
                          [15.0 ,22.5 ]])
    
    poly_P = numpy.array([[50.0 ,0.0  ],
                          [67.5 ,0.0  ],
                          [67.5 ,12.5 ],
                          [50.0 ,12.5 ],
                          [50.0 ,0.0  ]])
    
    poly_S = numpy.array([[67.5 ,0.0  ],
                          [67.5 ,17.5 ],
                          [82.5 ,17.5 ],
                          [92.5 ,7.5  ],
                          [82.5 ,7.5  ],
                          [82.5 ,0.0  ],
                          [67.5 ,0.0  ]])
    
    poly_Z = numpy.array([[82.5 ,0.0  ],
                          [82.5 ,7.5  ],
                          [92.5 ,7.5  ],
                          [100.0,0.0  ],
                          [82.5 ,0.0  ]])
    
    Poly_All_Labels=['U','E','A','L','P','S','Z']
    Poly_All=[poly_U,poly_E,poly_A,poly_L,poly_P,poly_S,poly_Z]
    Poly_All=map(lambda x : x/100.0,Poly_All)
    nP=len(Poly_All)
    
    Poly_Dict={}
    for iP in range(nP):
        Poly_Dict.update({Poly_All_Labels[iP] : Poly_All[iP]})
        
    return Poly_Dict

def GetTexture(Sand,Clay,DoPlot=False):
    
    import matplotlib.path as mplPath
    import matplotlib.patches as mplPatches

    Poly_Dict=DefineTexClassPolygons()
    Poly_All_Labels=Poly_Dict.keys()
    nP=len(Poly_Dict)
    
    #==============================================================================
    # Get Texture    
    #==============================================================================
    
    if numpy.size(Sand)<1:
        Sand=numpy.array([Sand,])
        Clay=numpy.array([Clay,])
    
    nTex=Sand.size
    TexClass=numpy.zeros(nTex,dtype='S3')
    TexClass[:]='nan'
    
    if numpy.size(Sand)>1:
        for iP in range(nP):
            TexLabel=Poly_All_Labels[iP]
            bbPath = mplPath.Path(Poly_Dict[TexLabel])
            TexClass[bbPath.contains_points(numpy.array([Sand,Clay]).transpose())]=TexLabel
    else:
        TexClass=['nan',]
    #==============================================================================
    # Draw
    #==============================================================================
    
    if DoPlot:
        
        fig=pylab.figure()
        ax=fig.add_subplot(111)
        
        for iP in range(nP):
            TexLabel=Poly_All_Labels[iP]
            bbPath = mplPath.Path(Poly_Dict[TexLabel])
            patch=mplPatches.PathPatch(bbPath,facecolor='w', lw=2)
        
            ax.add_patch(patch)
            ax.text(Poly_Dict[TexLabel].mean(axis=0)[0],Poly_Dict[TexLabel].mean(axis=0)[1],TexLabel)
        
        ax.plot(Sand,Clay,'r+')
        
        ax.set_aspect('equal')
        ax.set_xlabel('Sand (kg/kg)')
        ax.set_ylabel('Clay (kg/kg)')

        return TexClass,fig

    else:
        return TexClass

def DrawTexTriangle(Poly_Dict,fig=-1,ax=-1):

    import matplotlib.path as mplPath
    import matplotlib.patches as mplPatches
    
    if ax==-1:
        if fig==-1:
            fig=pylab.figure()
        ax=fig.add_subplot(111)
    
    Poly_All_Labels=Poly_Dict.keys()
    nP=len(Poly_Dict)
    
    for iP in range(nP):
        TexLabel=Poly_All_Labels[iP]
        bbPath = mplPath.Path(Poly_Dict[TexLabel])
        patch=mplPatches.PathPatch(bbPath,facecolor=None,fill=False, lw=2)
    
        ax.add_patch(patch)
        ax.text(Poly_Dict[TexLabel].mean(axis=0)[0],Poly_Dict[TexLabel].mean(axis=0)[1],TexLabel)
    
    ax.set_aspect('equal')
    ax.set_xlabel('Sand (kg/kg)')
    ax.set_ylabel('Clay (kg/kg)')
    
    return fig,ax


def TexCode2Name(TexCode):
    TexCode=numpy.array(TexCode)
    TexName=numpy.zeros_like(TexCode).astype('S16')
    TexCode_unique=numpy.array(['U','E','A','L','P','S','Z'])
    TexName_unique=numpy.array(['Zware Klei','Klei','Leem','Zandleem','Lichte Zandleem','Lemig Zand','Zand'])
    
    nTex=len(TexCode_unique)
    
    for iT in range(nTex):
        TexName[TexCode==[TexCode_unique[iT]]]=TexName_unique[iT]
    
    return TexName
    