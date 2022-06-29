# -*- coding: utf-8 -*-
"""
Created on Mon Jan 25 12:19:44 2016

@author: Administrator
"""

#==============================================================================
# Open output of a lifetime
#==============================================================================

import os
import pylab
import numpy
import copy 
import fnmatch

from Hydrus_Funky import *
from DIY_Hydrus_Funky import *
from RetentionConductivityCapacity_Funky import *


from matplotlib.backends.backend_pdf import PdfPages
figlist=[]
figsize1=[10./numpy.sqrt(2),10]
figsize2=[10,10./numpy.sqrt(2)]
#dpi=300
dpi=None

font = {'family' : 'monospace',
        'size'   : 8}
pylab.rc('font', **font)

pylab.ioff()

cm=pylab.get_cmap('Paired')

#==============================================================================
# Open that dataset
#==============================================================================

Root='C:\Users\Administrator\Documents\Arbeid\Nevenactiviteiten\Hydrus_MakeItRain'

filelist=[]
for root, dirnames, filenames in os.walk(Root):
  for filename in fnmatch.filter(filenames, '*.npz'):
      filelist.append(os.path.join(root, filename))
      
filename=filelist[0]
Savez=numpy.load(filename)

h=Savez['h']
theta=Savez['theta']
t=Savez['t']
dt=Savez['dt']
nit=Savez['nit']
H=Savez['H']
z=Savez['z']
q_BCtop=Savez['q_BCtop']
q_runoff=Savez['q_runoff']
h_SL=Savez['h_SL']
BC_top=Savez['BC_top']
q_top=Savez['q_top']
q_bottom=Savez['q_bottom']
h_BCbottom=Savez['h_BCbottom']
ModelSWRC=Savez['ModelSWRC']
ParSWRC=Savez['ParSWRC']
ModelKh=Savez['ModelKh']
ParKh=Savez['ParKh']
ModelCap=Savez['ModelCap']
ParCap=Savez['ParCap']

z_gwt=h[:,-1]
filter_GWT=z[None,:]<=z[-1]+z_gwt[:,None]
h_filt=copy.deepcopy(h)
h_filt[filter_GWT]=numpy.nan

allh=h_filt.flatten()
allh=allh[~numpy.isnan(allh)]
allh[allh>=-1e-1]=-1e-1

FluxAll=Flux_2D(z,h,dt,1,h2theta_BimodalPDIRet_FixedAlpha,ParSWRC,h2K_BimodalPDIK_FixedAlpha,ParKh)
t12=(t[1:]+t[:-1])/2
#==============================================================================
# Plot
#==============================================================================

fig=pylab.figure()
ax=fig.add_subplot(111)
#ax.imshow(h.transpose(),vmin=-100,vmax=5,interpolation='nearest',cmap='viridis')
im=ax.contourf(t,z,h.transpose(),numpy.linspace(-120,10,11),cmap='viridis_r')
pylab.rcParams['contour.negative_linestyle'] = 'solid'
ax.contour(t,z,h.transpose(),[0,],vmin=-100,vmax=10,colors='k')
ax.axis('tight')
ax.set_title('Tensiometric Head')
fig.colorbar(im)
figlist.append(fig)

fig=pylab.figure()
ax=fig.add_subplot(111)
#ax.imshow(h.transpose(),vmin=-100,vmax=5,interpolation='nearest',cmap='viridis')
im=ax.contourf(t,z,theta.transpose(),20,cmap='viridis_r')
pylab.rcParams['contour.negative_linestyle'] = 'solid'
ax.contour(t,z,theta.transpose(),[max(theta.flatten())-0.01,],colors='k')
ax.axis('tight')
ax.set_title('Water Content')
fig.colorbar(im)
figlist.append(fig)

fig=pylab.figure()
ax=fig.add_subplot(111)
#ax.imshow(h.transpose(),vmin=-100,vmax=5,interpolation='nearest',cmap='viridis')
im=ax.contourf(t12,z,FluxAll.transpose(),20,cmap='viridis_r',vmax=2e-5,vmin=-2e-5)
pylab.rcParams['contour.negative_linestyle'] = 'solid'
ax.contour(t12,z,FluxAll.transpose(),[0],colors='k')
ax.axis('tight')
ax.set_title('Flux')
fig.colorbar(im)
figlist.append(fig)


fig=pylab.figure()
ax=fig.add_subplot(111)
ax.hist(-allh,bins=numpy.logspace(-1,5,50),log=True)
ax.set_yscale('log')
ax.set_xscale('log')
ax.set_xlabel('h positivized (cm)')
figlist.append(fig)



h_filt[h_filt>=-1e-1]=-1e-1
fig=pylab.figure()
ax=fig.add_subplot(111)
#ax.plot((h_filt[1:,:]+h_filt[:-1,:])/2.0,dt,'k.',alpha=0.05)
im=ax.hexbin((-(h_filt[1:,:]+h_filt[:-1,:])/2.0).flatten(),(dt[:,None]*numpy.ones(z.size)[None,:]).flatten(), bins='log',xscale='log',yscale='log',cmap='binary',vmin=1.5,vmax=4)
#ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_title('h - dt')
ax.set_xlabel('h positivized (cm)')
ax.set_ylabel('dt (s)')
fig.colorbar(im)
figlist.append(fig)


h_filt[h_filt>=-1e-1]=-1e-1
fig=pylab.figure()
ax=fig.add_subplot(111)
#ax.plot((h_filt[1:,:]+h_filt[:-1,:])/2.0,dt,'k.',alpha=0.05)
im=ax.hexbin((-(h_filt[1:,:]+h_filt[:-1,:])/2.0).flatten(),(abs(FluxAll)).flatten(), bins='log',xscale='log',yscale='log',cmap='binary',vmin=1.5,vmax=4)
#ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_title('q - dt')
ax.set_xlabel('h positivized (cm)')
ax.set_ylabel('q (cm/s)')
fig.colorbar(im)
figlist.append(fig)


h_filt[h_filt>=-1e-1]=-1e-1
fig=pylab.figure()
ax=fig.add_subplot(111)
#ax.plot((h_filt[1:,:]+h_filt[:-1,:])/2.0,dt,'k.',alpha=0.05)
im=ax.hexbin((-(h_filt[1:,:]+h_filt[:-1,:])/2.0).flatten(),(abs(dt[:,None]*FluxAll)).flatten(), bins='log',xscale='log',yscale='log',cmap='binary',vmin=1.5,vmax=4)
#ax.set_xscale('log')
#ax.set_yscale('log')
ax.set_title('V - dt')
ax.set_xlabel('h positivized (cm)')
ax.set_ylabel('q * dt (cm)')
fig.colorbar(im)
figlist.append(fig)


##fig=pylab.figure()
##ax=fig.add_subplot(111)
###ax.imshow(h.transpose(),vmin=-100,vmax=5,interpolation='nearest',cmap='viridis')
##im=ax.contourf(t,z,h.transpose(),numpy.linspace(-120,10,11),cmap='viridis')
##pylab.rcParams['contour.negative_linestyle'] = 'solid'
##ax.contour(t,z,h.transpose(),[0,],vmin=-100,vmax=10,colors='k')
##ax.axis('tight')
##fig.colorbar(im)
##figlist.append(fig)
#
#
#fig=pylab.figure()
#ax=fig.add_subplot(211)
#ax.plot(nit,'k-',label='Replicate')
#ax.set_ylim([0,None])
#ax.set_title('number of iterations')
#ax=fig.add_subplot(212)
#ax.plot(dt,'k-',label='Replicate')
#ax.set_yscale('log')
#ax.set_title('dt')
#figlist.append(fig)
#
#fig=pylab.figure()
#ax=fig.add_subplot(111)
#ax.plot(t[1:],q_top,'k-',label='Replicate')
#ax.plot(t[1:],q_BCtop,'b-',label='BC')
##ax.plot(t_hydrus,q_top_hydrus,'r.-',label='Hydrus')
#ax.set_xlabel(TLevel_Header[0])
#ax.set_ylabel(TLevel_Header[3])
#ax.legend()
#ax.set_title(TLevel_Header[3])
#figlist.append(fig)
#
#fig=pylab.figure()
#ax=fig.add_subplot(111)
#ax.plot(t[1:],q_top,'k-',label='Replicate')
#ax.plot(t[1:],q_BCtop,'b-',label='BC')
##ax.plot(t_hydrus,q_top_hydrus,'r.-',label='Hydrus')
#ax.set_xlabel(TLevel_Header[0])
#ax.set_ylabel(TLevel_Header[3])
#ax.set_ylim([-2e-4,5e-5])
#ax.legend()
#ax.set_title(TLevel_Header[3])
#figlist.append(fig)
#
#fig=pylab.figure()
#ax=fig.add_subplot(111)
#ax.plot(t[1:],q_bottom,'k-',label='Replicate')
##ax.plot(t_hydrus,q_bottom_hydrus,'r.-',label='Hydrus')
#ax.set_xlabel(TLevel_Header[0])
#ax.set_ylabel(TLevel_Header[5])
#ax.legend()
#ax.set_title(TLevel_Header[5])
#figlist.append(fig)
#
#
#fig=pylab.figure()
#ax=fig.add_subplot(111)
#ax.plot(t[1:],q_runoff,'k-',label='Replicate')
##ax.plot(t_hydrus,q_runoff_hydrus,'r.-',label='Hydrus')
#ax.set_xlabel(TLevel_Header[0])
#ax.set_ylabel(TLevel_Header[14])
#ax.legend()
#ax.set_title(TLevel_Header[14])
#figlist.append(fig)
#
#
#
##nP_Hydrus=NodInf_T.size
##for iP in range(nP_Hydrus):
##    fig=pylab.figure()
##    ax=fig.add_subplot(111)
##    ax.plot(h[t==NodInf_T[iP],:][0,:],z,'k.-',label='Replicate')
##    ax.plot(NodInf_Data[iP][:,2],NodInf_Data[iP][:,1],'r.-',label='Hydrus')
##    ax.set_ylabel(NodInf_Header[0])
##    ax.set_xlabel(NodInf_Header[2])
##    ax.set_title('t = %s s'%NodInf_T[iP])
##    ax.legend()
##
##
##nP_Hydrus=NodInf_T.size
##for iP in range(nP_Hydrus):
##    fig=pylab.figure()
##    ax=fig.add_subplot(111)
##    ax.plot(theta[t==NodInf_T[iP],:][0,:],z,'k.-',label='Replicate')
##    ax.plot(NodInf_Data[iP][:,3],NodInf_Data[iP][:,1],'r.-',label='Hydrus')
##    ax.set_ylabel(NodInf_Header[1])
##    ax.set_xlabel(NodInf_Header[3])
##    ax.set_title('t = %s s'%NodInf_T[iP])
##    ax.legend()
#
#nO=ObsNode_Nodes.size
#for iO in range(nO):
#    fig=pylab.figure(figsize=figsize1)
#    ax=fig.add_subplot(111)
#    iZ=ObsNode_Nodes[iO]-1
#    ax.plot(t,h[:,iZ],'k-',label='Replicate')
#    #ax.plot(ObsNode_t,ObsNode_L[iO][:,0],'r-',label='Hydrus')
#    ax.set_xlabel('Time')
#    ax.set_ylabel('Head')
#    ax.legend(loc=4)
#    #ax.set_xlim([0,.1])
#    #ax.set_ylim([-1000,0])
#    ax.set_title('Observation node z = %s cm' %z[iZ])
#    figlist.append(fig)
#
#
#fig=pylab.figure()
#ax=fig.add_subplot(111)
#ax.axhline(0,color='k')
#ax.plot(t[1:][BC_top=='Dirichelet'],h_SL[:-1][BC_top=='Dirichelet'],'r.',label='Dirichelet')
#ax.plot(t[1:][BC_top=='Neumann'],h_SL[:-1][BC_top=='Neumann'],'b.',label='Neumann')
#ax.plot(t,h[:,0],'k--')
#ax.set_xlabel('h Surface Layer (cm)')
#ax.set_ylabel('Head')
#ax.set_ylim([-1e-2,1e-1])
#ax.legend()
#figlist.append(fig)
#

#==============================================================================
# Save
#==============================================================================

## Save data
import sys
basename=os.path.basename(sys.argv[0])[:-3]

#postfix='_%s_%s_%s_%s'%(FC,SFC,HC,K)
postfix=''


## save plots to pdf


pdfname = os.path.join(Root, basename + postfix + '.pdf')
pp = PdfPages(pdfname)
for fig in figlist:
    pp.savefig(fig)

pp.close()

#extension='.png'
#for iF in range(len(figlist)):
#    fig=figlist[iF]
#    figname = os.path.join(Root, basename+'_%s'%iF + postfix + extension)
#    fig.savefig(figname, bbox_inches='tight')

pylab.show()
#os.startfile(pdfname)
