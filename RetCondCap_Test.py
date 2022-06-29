# -*- coding: utf-8 -*-
"""
Created on Thu Oct 22 11:28:40 2015

@author: Administrator
"""

#==============================================================================
# Testing some hydrus functions
#==============================================================================

import numpy
import pylab

from RetentionConductivityCapacity_Funky import *

#==============================================================================
# Parameters
#==============================================================================

#h=numpy.linspace(-1000,100,10000)
h=-numpy.logspace(-1,7,10000)
h=-numpy.logspace(-4,1,10000)
#h=numpy.ones((100))*-100
hplot=-h

thetar=0.1
thetas=0.6
Ks=100.0
l=0.5
w1=0.2
omega1=0.1
alpha1=0.02
n1=1.2
m1=1+1.0/n1
alpha2=0.02
n2=1.2
m2=1+1.0/n2
pF0=6.8
a=-1.5

[thetar,thetas,Ks,l,w1,omega1,alpha1,n1,m1,alpha2,n2,m2,pF0,a]=map(lambda x : numpy.array([x,]),[thetar,thetas,Ks,l,w1,omega1,alpha1,n1,m1,alpha2,n2,m2,pF0,a])

#thetar=0.1*numpy.ones((100))
#thetas=0.6*numpy.ones((100))
#thetas[:50]=0.4
#Ks=100.0*numpy.ones((100))
#Ks[:50]=50
#l=0.5*numpy.ones((100))
#w1=0.2*numpy.ones((100))
#omega1=0.1*numpy.ones((100))
#alpha1=0.02*numpy.ones((100))
#n1=1.5*numpy.ones((100))
#m1=1+1.0/n1*numpy.ones((100))
#alpha2=0.2*numpy.ones((100))
#n2=2.5*numpy.ones((100))
#m2=1+1.0/n2*numpy.ones((100))
#pF0=100*numpy.ones((100))
#a=-1.5*numpy.ones((100))

#==============================================================================
# Test
#==============================================================================

### Retention
#
#fig=pylab.figure()
#ax=fig.add_subplot(111)
#
#Par=[thetar,thetas,alpha1,n1]
#theta=h2theta_VanGenuchten4(h,Par)
##print(theta)
#ax.plot(hplot,theta)
##ax.plot(theta)
#
##Par=[thetar,thetas,w1,alpha1,n1,alpha2,n2]
##theta=h2theta_Durner(h,Par)
###print(theta)
##ax.plot(hplot,theta)
###ax.plot(theta)
#
##Par=[thetar,thetas,alpha1,n1,pF0]
##theta=h2theta_PDI_AdsorptiveSaturation(h,Par)
###print(theta)
##ax.plot(hplot,theta)
###ax.plot(theta)
#
#Par=[thetar,thetas,alpha1,n1,pF0]
#theta=h2theta_PDIRet(h,Par)
##print(theta)
#ax.plot(hplot,theta)
##ax.plot(theta)
#
#Par=[thetar,thetas,w1,alpha1,n1,alpha2,n2,pF0]
#theta=h2theta_BimodalPDIRet(h,Par)
##print(theta)
#ax.plot(hplot,theta)
##ax.plot(theta)
#
#ax.set_xscale('log')

## Conductivity

fig=pylab.figure()
ax=fig.add_subplot(111)

Par=[Ks,l,alpha1,n1,m1]
K=h2K_Mualem5(h,Par)
#print(K)
ax.plot(hplot,K)
#ax.plot(K)

#Par=[Ks,l,w1,alpha1,n1,alpha2,n2]
#K=h2K_MualemBimodal(h,Par)
##print(K)
#ax.plot(hplot,K)
##ax.plot(K)


#Par=[thetar,thetas,Ks,l,omega1,alpha1,n1,pF0,a]
#K=h2K_PDI_Kfilm(h,Par)
##print(K)
#ax.plot(hplot,K)
##ax.plot(K)


Par=[thetar,thetas,Ks,l,omega1,alpha1,n1,pF0,a]
K=h2K_PDIK(h,Par,VaporConductivity=False)
#print(K)
ax.plot(hplot,K)
#ax.plot(K)


#Par=[thetar,thetas,Ks,l,w1,omega1,alpha1,n1,alpha2,n2,pF0,a]
#K=h2K_BimodalPDI_Kfilm(h,Par)
##print(K)
#ax.plot(hplot,K)
##ax.plot(K)


Par=[thetar,thetas,Ks,l,w1,omega1,alpha1,n1,alpha2,n2,pF0,a]
K=h2K_BimodalPDIK(h,Par,VaporConductivity=False)
#print(K)
ax.plot(hplot,K)
#ax.plot(K)


#Par=[thetar,thetas,alpha1,n1,pF0]
#K=h2K_PDI_Kvap(h,Par,SWR_Model=h2theta_PDIRet)
##print(K)
#ax.plot(hplot,K)
##ax.plot(K)


#Par=[thetar,thetas,w1,alpha1,n1,alpha2,n2,pF0]
#K=h2K_PDI_Kvap(h,Par,SWR_Model=h2theta_BimodalPDIRet)
##print(K)
#ax.plot(hplot,K)
##ax.plot(K)

ax.set_xscale('log')
ax.set_yscale('log')


## Capacity

#fig=pylab.figure()
#ax=fig.add_subplot(111)
#
#
#Par=[thetar,thetas,alpha1,n1]
#C=h2C_VanGenuchten4(h,Par)
##print(C)
#ax.plot(hplot,C)
##ax.plot(C)
#
#
#Par=[thetar,thetas,alpha1,n1,pF0]
#C=h2C_PDICap(h,Par)
##print(C)
#ax.plot(hplot,C)
##ax.plot(C)
#
#Par=[thetar,thetas,w1,alpha1,n1,alpha2,n2,pF0]
#C=h2C_BimodalPDICap(h,Par)
##theta=h2theta_BimodalPDIRet(h,Par)
##dthetadh=(theta[1:]-theta[:-1])/(h[1:]-h[:-1])
##hC=(h[1:]+h[:-1])/2
##print(C)
#ax.plot(hplot,C)
##ax.plot(C)
##ax.plot(-hC,dthetadh,'k.')
#
#ax.set_xscale('log')



#==============================================================================
# 
#==============================================================================

#pylab.show()
