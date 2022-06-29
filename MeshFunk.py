# -*- coding: utf-8 -*-
"""
Created on Tue Jan  7 16:32:30 2014

@author: jan
"""

#==============================================================================
# Fun with distmesh
# Offers perspectives to make Hydrus 2D and 3D
#==============================================================================

import distmesh
import numpy
import pylab
import scipy.spatial as spatial


#def fd10(p):
#    r, z = numpy.sqrt(p[:,0]**2 + p[:,1]**2), p[:,2]
#    d1, d2, d3 = r-1.0, z-1.0, -z-1.0
#    d4, d5 = numpy.sqrt(d1**2+d2**2), numpy.sqrt(d1**2+d3**2)
#    d = distmesh.dintersect(distmesh.dintersect(d1, d2), d3)
#    ix = (d1>0)*(d2>0); d[ix] = d4[ix]
#    ix = (d1>0)*(d3>0); d[ix] = d5[ix]
#    return distmesh.ddiff(d, distmesh.dsphere(p, 0,0,0, 0.5))
#    
#def fh10(p):
#    h1 = 4*numpy.sqrt((p**2).sum(1))-1.0
#    return numpy.minimum(h1, 2.0)
#    
#p, t = distmesh.distmeshnd(fd10, fh10, 0.1, (-1,-1,-1, 1,1,1))

fd = lambda p: distmesh.ddiff(distmesh.drectangle(p,-1,1,-1,1),distmesh.dcircle(p,0,0,0.5))
#fd = lambda p: distmesh.drectangle(p,-1,1,-1,1)
#fh = lambda p: -abs(p[:,1]-p[:,0]*0)**0.5-1
#fh = lambda p: (p[:,-1]-0.5)**2
fh = lambda p: abs(p[:,1])+abs(p[:,0])

p, t = distmesh.distmesh2d(fd, fh, 0.05, (-1,-1,1,1),[(-1,-1),(-1,1),(1,-1),(1,1)])

pylab.plot(p[:,0],p[:,1],'ro')

Hull=spatial.ConvexHull(p)
Hull_X=map(lambda x: p[x[0]],Hull.simplices)
Hull_Y=map(lambda x: p[x[1]],Hull.simplices)
pylab.plot(Hull_X,Hull_Y,'*')

Vor=spatial.Voronoi(p)
spatial.voronoi_plot_2d(Vor)
pylab.plot(Hull_X,Hull_Y,'r*')

pylab.show()

