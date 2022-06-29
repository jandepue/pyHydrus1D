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
import datetime

#from Hydrus_Funky import *
from DIY_Hydrus_Funky import *
from RetentionConductivityCapacity_Funky import *


from matplotlib.backends.backend_pdf import PdfPages
figlist=[]
figlist2=[]
figsize1=[10./numpy.sqrt(2),12]
figsize2=[10,10./numpy.sqrt(2)]
#dpi=300
dpi=None

font = {'family' : 'monospace',
        'size'   : 8}
pylab.rc('font', **font)

pylab.ioff()

cm=pylab.get_cmap('Paired')
linestyle=['-','--',':','-.']

#==============================================================================
# Open that dataset
#==============================================================================

#Root='C:\Users\Administrator\Documents\Arbeid\Nevenactiviteiten\Hydrus_MakeItRain'
Root='/home/supersoil/Documents/Jan/AlbonHydrus_DoTheWalk'

filelist=[]
labellist=[]
filt='*Z_REF_H'
#filt=''
for root, dirnames, filenames in os.walk(Root):
  for filename in fnmatch.filter(filenames, filt+'*.npz'):
      filelist.append(os.path.join(root, filename))
      labellist.append(filename.replace('.npz','').replace('AlbonHydrus_DoTheWalk_',''))
filelist=numpy.array(filelist)
labellist=numpy.array(labellist)
nF=len(filelist)

labellist_FC=numpy.array(map(lambda x : '_'.join(x.split('_')[:-1]),labellist))
labellist_K=numpy.array(map(lambda x : x.split('_')[-1],labellist))

FC_unique=numpy.unique(labellist_FC)
nFC=FC_unique.size

for iFC in range(nFC):
    
    K_unique=numpy.unique(labellist_K[labellist_FC==FC_unique[iFC]])
    nK=K_unique.size
    
    CompareK_Theta=[]
    CompareK_h=[]
    CompareK_Models=[]
    CompareK_Par=[]
    
    for iK in range(nK):
        filename=filelist[(labellist_FC==FC_unique[iFC]) & (labellist_K==K_unique[iK])][0]
        label=labellist[(labellist_FC==FC_unique[iFC]) & (labellist_K==K_unique[iK])][0]  
        
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
        
        if ModelSWRC==h2theta_VanGenuchten5:
            ModelSWRC=h2theta_VanGenuchten5
            ModelKh=h2K_Mualem5
            ModelCap=h2C_VanGenuchten5
        elif ModelSWRC==h2theta_BimodalPDIRet_FixedAlpha:
            ModelSWRC=h2theta_BimodalPDIRet_FixedAlpha
            ModelKh=h2K_BimodalPDIK_FixedAlpha
            ModelCap=h2C_BimodalPDICap_FixedAlpha
        
        iLayer=[]
        LayerLabel=['TOP','cSUB','SUB']
        for iZ in range(3):
            iLayer.append(numpy.argmax(ParSWRC[2,:]==numpy.unique(ParSWRC[2,:])[iZ]))
        iLayer.sort()
        nL=len(iLayer)
        
#        incl=0.0
#        FluxAll=Flux_2D(z,h,dt,incl,ModelSWRC,ParSWRC,ModelKh,ParKh)
#        t12=(t[1:]+t[:-1])/2
        
        #==============================================================================
        # Validation Data    
        #==============================================================================
        t0_sim='12/10/2014 22:00'
        t0_sim=datetime.datetime.strptime(t0_sim,'%d/%m/%Y %H:%M')
        validationname='/home/supersoil/Documents/Jan/git/Hydrus_DIY/DataLanduyt/Landuyt_SWC_Validation.txt'
        validation=numpy.genfromtxt(validationname,dtype='str',delimiter='\t',skip_header=1)
        validationt=numpy.array(map(lambda x: datetime.datetime.strptime(x,'%d/%m/%Y %H:%M'), validation[:,0]))
        validationdt=validationt-t0_sim
        validationdt_s=numpy.array(map(lambda x: x.seconds+x.days*24.0*3600, validationdt))
        validationswc=validation[:,1:].astype(float)
        filt=validationdt_s>0
        
        validationswc=validationswc[filt,:]
        validationdt_s=validationdt_s[filt]
        
        validationZ=[-10.0,-30.0,-80.0]
        validationiZ=[]
        validationiT=[]
        nV=len(validationZ)
        for iV in range(nV):
            validationiZ.append(numpy.argmin((z-validationZ[iV])**2))
            
        for T in validationdt_s:
            if t[-1]>=T:
                validationiT.append(numpy.argmin((t-T)**2))
            else:
                validationiT.append(-1)
                
        validationiT=numpy.array(validationiT)
        validation_ModelTheta=theta[validationiT,:]
        validation_ModelTheta[validationiT==-1]=numpy.nan
        validation_Modelh=h[validationiT,:]
        validation_Modelh[validationiT==-1]=numpy.nan
        
        
        #==============================================================================
        # Save
        #==============================================================================
        
        CompareK_Theta.append(validation_ModelTheta)
        CompareK_h.append(validation_Modelh)
        ComarpeK_T=validationdt_s
        CompareK_Models.append([ModelSWRC,ModelKh,ModelCap])
        CompareK_Par.append([ParSWRC,ParKh,ParCap])
        
        #==============================================================================
        # Plot
        #==============================================================================
        
#        fig=pylab.figure()
#        fig.suptitle(label)
#        ax=fig.add_subplot(111)
#        #ax.imshow(h.transpose(),vmin=-100,vmax=5,interpolation='nearest',cmap='viridis')
#        im=ax.contourf(t,z,h.transpose(),numpy.linspace(-120,10,11),cmap='viridis_r')
#        pylab.rcParams['contour.negative_linestyle'] = 'solid'
#        ax.contour(t,z,h.transpose(),[0,],vmin=-100,vmax=10,colors='k')
#        ax.axis('tight')
#        ax.set_title('Tensiometric Head')
#        fig.colorbar(im)
#        figlist.append(fig)
#        
#        fig=pylab.figure()
#        fig.suptitle(label)
#        ax=fig.add_subplot(111)
#        #ax.imshow(h.transpose(),vmin=-100,vmax=5,interpolation='nearest',cmap='viridis')
#        im=ax.contourf(t,z,theta.transpose(),20,cmap='viridis_r')
#        pylab.rcParams['contour.negative_linestyle'] = 'solid'
#        ax.contour(t,z,theta.transpose(),[max(theta.flatten())-0.01,],colors='k')
#        ax.axis('tight')
#        ax.set_title('Water Content')
#        fig.colorbar(im)
#        figlist.append(fig)
#        
#        fig=pylab.figure()
#        fig.suptitle(label)
#        ax=fig.add_subplot(111)
#        #ax.imshow(h.transpose(),vmin=-100,vmax=5,interpolation='nearest',cmap='viridis')
#        im=ax.contourf(t12,z,FluxAll.transpose(),[-1,-1e-5,1e-5,1],cmap='bwr')
#        im=ax.contourf(t12,z,FluxAll.transpose(),numpy.linspace(-1e-5,1e-5,20),cmap='bwr')
#    #    pylab.rcParams['contour.negative_linestyle'] = 'solid'
#    #    ax.contour(t12,z,FluxAll.transpose(),[0,],colors='k')
#        ax.axis('tight')
#        ax.set_title('Flux')
#        fig.colorbar(im)
#        figlist.append(fig)
#        
#        for iV in range(nV):
#            fig=pylab.figure()
#            fig.suptitle(label)
#            ax=fig.add_subplot(111)
#            ax.plot(t,theta[:,validationiZ[iV]],'k')
#            ax.plot(validationdt_s,validationswc[:,iV],'r')
#            ax.set_title(validationZ[iV])
#            ax.set_ylabel('SWC (m3/m3)')
#            ax.set_xlabel('Time (h)')
#            figlist.append(fig)
#    
#        cm=pylab.get_cmap('viridis')
#        fig=pylab.figure()
#        fig.suptitle(label)
#        ax=fig.add_subplot(111)
#        ax.plot([0,0.5,],[0,0.5,],'r-')
#        for iV in range(nV):
#            color=cm(iV/(nV-0.99999))
#            ax.plot(theta[:,validationiZ[iV]][validationiT],validationswc[:,iV],'o',color=color,mec=None,alpha=0.1)
#        ax.axis('equal')
#        ax.set_title(validationZ[iV])
#        ax.set_ylabel('SWC (m3/m3)')
#        ax.set_xlabel('Time (h)')
#        figlist.append(fig)
#        
#        fig=pylab.figure()
#        fig.suptitle(label)
#        ax=fig.add_subplot(111)
#        ax.hist(-allh,bins=numpy.logspace(-1,5,50),log=True)
#        ax.set_yscale('log')
#        ax.set_xscale('log')
#        ax.set_xlabel('h positivized (cm)')
#        figlist.append(fig)
#        
#        h_filt[h_filt>=-1e-1]=-1e-1
#        fig=pylab.figure()
#        fig.suptitle(label)
#        ax=fig.add_subplot(111)
#        #ax.plot((h_filt[1:,:]+h_filt[:-1,:])/2.0,dt,'k.',alpha=0.05)
#        im=ax.hexbin((-(h_filt[1:,:]+h_filt[:-1,:])/2.0).flatten(),(dt[:,None]*numpy.ones(z.size)[None,:]).flatten(), bins='log',xscale='log',yscale='log',cmap='binary',vmin=1.5,vmax=4)
#        #ax.set_xscale('log')
#        #ax.set_yscale('log')
#        ax.set_title('h - dt')
#        ax.set_xlabel('h positivized (cm)')
#        ax.set_ylabel('dt (s)')
#        ax.set_xlim([1e-2,1e6])
#        fig.colorbar(im)
#        figlist.append(fig)
#        
#        
#        h_filt[h_filt>=-1e-1]=-1e-1
#        fig=pylab.figure()
#        fig.suptitle(label)
#        ax=fig.add_subplot(111)
#        #ax.plot((h_filt[1:,:]+h_filt[:-1,:])/2.0,dt,'k.',alpha=0.05)
#        im=ax.hexbin((-(h_filt[1:,:]+h_filt[:-1,:])/2.0).flatten(),(abs(FluxAll)).flatten()+1e-10, bins='log',xscale='log',yscale='log',cmap='binary',vmin=1.5,vmax=4)
#    #    im=ax.hexbin((-(h_filt[1:,:]+h_filt[:-1,:])/2.0).flatten(),(FluxAll).flatten(), bins='log',xscale='log',cmap='binary',vmin=1.5,vmax=4)
#        #ax.set_xscale('log')
#        #ax.set_yscale('log')
#        ax.set_title('q - dt')
#        ax.set_xlabel('h positivized (cm)')
#        ax.set_ylabel('q (cm/s)')
#        ax.set_xlim([1e-2,1e6])
#        fig.colorbar(im)
#        figlist.append(fig)
#        
#        
#        h_filt[h_filt>=-1e-1]=-1e-1
#        fig=pylab.figure()
#        fig.suptitle(label)
#        ax=fig.add_subplot(111)
#        #ax.plot((h_filt[1:,:]+h_filt[:-1,:])/2.0,dt,'k.',alpha=0.05)
#        im=ax.hexbin((-(h_filt[1:,:]+h_filt[:-1,:])/2.0).flatten(),(abs(dt[:,None]*FluxAll)).flatten()+1e-10, bins='log',xscale='log',yscale='log',cmap='binary',vmin=1.5,vmax=4)
#        #ax.set_xscale('log')
#        #ax.set_yscale('log')
#        ax.set_title('V - dt')
#        ax.set_xlabel('h positivized (cm)')
#        ax.set_ylabel('q * dt (cm)')
#        ax.set_xlim([1e-2,1e6])
#        fig.colorbar(im)
#        figlist.append(fig)
#        
#        #
#        #fig=pylab.figure()
#        #ax=fig.add_subplot(211)
#        #ax.plot(nit,'k-',label='Replicate')
#        #ax.set_ylim([0,None])
#        #ax.set_title('number of iterations')
#        #ax=fig.add_subplot(212)
#        #ax.plot(dt,'k-',label='Replicate')
#        #ax.set_yscale('log')
#        #ax.set_title('dt')
#        #figlist.append(fig)
#        #
#        #fig=pylab.figure()
#        #ax=fig.add_subplot(111)
#        #ax.plot(t[1:],q_top,'k-',label='Replicate')
#        #ax.plot(t[1:],q_BCtop,'b-',label='BC')
#        ##ax.plot(t_hydrus,q_top_hydrus,'r.-',label='Hydrus')
#        #ax.set_xlabel(TLevel_Header[0])
#        #ax.set_ylabel(TLevel_Header[3])
#        #ax.legend()
#        #ax.set_title(TLevel_Header[3])
#        #figlist.append(fig)
#        #
#        #fig=pylab.figure()
#        #ax=fig.add_subplot(111)
#        #ax.plot(t[1:],q_top,'k-',label='Replicate')
#        #ax.plot(t[1:],q_BCtop,'b-',label='BC')
#        ##ax.plot(t_hydrus,q_top_hydrus,'r.-',label='Hydrus')
#        #ax.set_xlabel(TLevel_Header[0])
#        #ax.set_ylabel(TLevel_Header[3])
#        #ax.set_ylim([-2e-4,5e-5])
#        #ax.legend()
#        #ax.set_title(TLevel_Header[3])
#        #figlist.append(fig)
#        #
#        #fig=pylab.figure()
#        #ax=fig.add_subplot(111)
#        #ax.plot(t[1:],q_bottom,'k-',label='Replicate')
#        ##ax.plot(t_hydrus,q_bottom_hydrus,'r.-',label='Hydrus')
#        #ax.set_xlabel(TLevel_Header[0])
#        #ax.set_ylabel(TLevel_Header[5])
#        #ax.legend()
#        #ax.set_title(TLevel_Header[5])
#        #figlist.append(fig)
#        #
#        #
#        #fig=pylab.figure()
#        #ax=fig.add_subplot(111)
#        #ax.plot(t[1:],q_runoff,'k-',label='Replicate')
#        ##ax.plot(t_hydrus,q_runoff_hydrus,'r.-',label='Hydrus')
#        #ax.set_xlabel(TLevel_Header[0])
#        #ax.set_ylabel(TLevel_Header[14])
#        #ax.legend()
#        #ax.set_title(TLevel_Header[14])
#        #figlist.append(fig)
#        #
#        #
#        #
#        ##nP_Hydrus=NodInf_T.size
#        ##for iP in range(nP_Hydrus):
#        ##    fig=pylab.figure()
#        ##    ax=fig.add_subplot(111)
#        ##    ax.plot(h[t==NodInf_T[iP],:][0,:],z,'k.-',label='Replicate')
#        ##    ax.plot(NodInf_Data[iP][:,2],NodInf_Data[iP][:,1],'r.-',label='Hydrus')
#        ##    ax.set_ylabel(NodInf_Header[0])
#        ##    ax.set_xlabel(NodInf_Header[2])
#        ##    ax.set_title('t = %s s'%NodInf_T[iP])
#        ##    ax.legend()
#        ##
#        ##
#        ##nP_Hydrus=NodInf_T.size
#        ##for iP in range(nP_Hydrus):
#        ##    fig=pylab.figure()
#        ##    ax=fig.add_subplot(111)
#        ##    ax.plot(theta[t==NodInf_T[iP],:][0,:],z,'k.-',label='Replicate')
#        ##    ax.plot(NodInf_Data[iP][:,3],NodInf_Data[iP][:,1],'r.-',label='Hydrus')
#        ##    ax.set_ylabel(NodInf_Header[1])
#        ##    ax.set_xlabel(NodInf_Header[3])
#        ##    ax.set_title('t = %s s'%NodInf_T[iP])
#        ##    ax.legend()
#        #
#        #nO=ObsNode_Nodes.size
#        #for iO in range(nO):
#        #    fig=pylab.figure(figsize=figsize1)
#        #    ax=fig.add_subplot(111)
#        #    iZ=ObsNode_Nodes[iO]-1
#        #    ax.plot(t,h[:,iZ],'k-',label='Replicate')
#        #    #ax.plot(ObsNode_t,ObsNode_L[iO][:,0],'r-',label='Hydrus')
#        #    ax.set_xlabel('Time')
#        #    ax.set_ylabel('Head')
#        #    ax.legend(loc=4)
#        #    #ax.set_xlim([0,.1])
#        #    #ax.set_ylim([-1000,0])
#        #    ax.set_title('Observation node z = %s cm' %z[iZ])
#        #    figlist.append(fig)
#        #
#        #
#        #fig=pylab.figure()
#        #ax=fig.add_subplot(111)
#        #ax.axhline(0,color='k')
#        #ax.plot(t[1:][BC_top=='Dirichelet'],h_SL[:-1][BC_top=='Dirichelet'],'r.',label='Dirichelet')
#        #ax.plot(t[1:][BC_top=='Neumann'],h_SL[:-1][BC_top=='Neumann'],'b.',label='Neumann')
#        #ax.plot(t,h[:,0],'k--')
#        #ax.set_xlabel('h Surface Layer (cm)')
#        #ax.set_ylabel('Head')
#        #ax.set_ylim([-1e-2,1e-1])
#        #ax.legend()
#        #figlist.append(fig)
#        #
#        pylab.close('all')
    
    Tensionplot=-numpy.logspace(-2,7,1000)
    cm=pylab.get_cmap('Paired')
    for iL in range(nL):
        fig=pylab.figure(figsize=figsize1)
        ax21=fig.add_subplot(311)
        ax22=fig.add_subplot(312)
        ax23=fig.add_subplot(313)
        ax21.set_title('Soil Water Retention')
        ax22.set_title('Hydraulic Conductivity')
        ax23.set_title('Capacity')
        ax21.set_ylabel('Soil Water Content (m3/m3)')
        ax22.set_ylabel('Hydraulic Conductivity (cm/d)')
        ax23.set_xlabel('Matric Potential (cm)')
        ax23.set_ylabel('Capacity (m3/m3/cm)')
        ax23.set_xlabel('Matric Potential (cm)')
#        ax21.set_xlim([1e-2,1e5])
        ax21.set_ylim([0.0,0.5])
#        ax22.set_xlim([1e-2,1e5])
        ax22.set_ylim([1e-15,1e0])
#        ax23.set_xlim([1e-2,1e5])
        ax23.set_ylim([0.0,6e-3])
        ax21.set_xscale('log')
        ax22.set_xscale('log')
        ax23.set_xscale('log')
        ax22.set_yscale('log')
        for iK in range(nK):
#            color=cm(float(iK)/(nK-0.999))
            color='k'
            ls=linestyle[iK]
            ax21.plot(-Tensionplot ,CompareK_Models[iK][0](Tensionplot,map(lambda x : x[iLayer[iL]],CompareK_Par[iK][0])),ls=ls,color=color,label='%s'%K_unique[iK])
            ax22.plot(-Tensionplot ,CompareK_Models[iK][1](Tensionplot,map(lambda x : x[iLayer[iL]],CompareK_Par[iK][1])),ls=ls,color=color,label='%s'%K_unique[iK])
            ax23.plot(-Tensionplot ,CompareK_Models[iK][2](Tensionplot,map(lambda x : x[iLayer[iL]],CompareK_Par[iK][2])),ls=ls,color=color,label='%s'%K_unique[iK])
        ax21.legend(loc=1)
        ax22.legend(loc=1)
        ax23.legend(loc=1)
        fig.suptitle('%s_%s'%(FC_unique[iFC],LayerLabel[iL]))
        figlist2.append(fig)
    
    zplot=[0,-10,-30,-80]
    nZplot=len(zplot)
    fig=pylab.figure(figsize=figsize1)
    for iZ in range(nZplot):
        ax=fig.add_subplot(nZplot,1,iZ+1)
        ax.set_title('z = %s cm'%zplot[iZ])
        ax.set_xlabel('time (h)')
        ax.set_ylabel('h (cm)')
        for iK in range(nK):
#            color=cm(float(iK)/(nK-0.999))
            color='k'
            ls=linestyle[iK]
            ax.plot(ComarpeK_T ,CompareK_h[iK][:,numpy.argmin((z-zplot[iZ])**2)],ls=ls,color=color,label='%s'%K_unique[iK])
        ax.legend(loc=1)
    fig.suptitle('%s_%s'%(FC_unique[iFC],LayerLabel[iL]))
    figlist2.append(fig)
    
    zplot=[0,-10,-30,-80]
    nZplot=len(zplot)
    fig=pylab.figure(figsize=figsize1)
    for iZ in range(nZplot):
        ax=fig.add_subplot(nZplot,1,iZ+1)
        ax.set_title('z = %s cm'%zplot[iZ])
        ax.set_xlabel('time (h)')
        ax.set_ylabel('theta (m3/m3)')
        for iK in range(nK):
            color=cm(float(iK)/(nK-0.999))
            ls='-'
#            color='k'
#            ls=linestyle[iK]
            ax.plot(ComarpeK_T ,CompareK_Theta[iK][:,numpy.argmin((z-zplot[iZ])**2)],ls=ls,color=color,label='%s'%K_unique[iK])
        ax.legend(loc=1)
    fig.suptitle('%s_%s'%(FC_unique[iFC],LayerLabel[iL]))
    figlist2.append(fig)
    
    zplot=[-10,-30,-80]
    nZplot=len(zplot)
    fig=pylab.figure(figsize=figsize1)
    for iZ in range(nZplot):
        ax=fig.add_subplot(nZplot,1,iZ+1)
        ax.set_title('z = %s cm'%zplot[iZ])
        ax.set_xlabel('time (h)')
        ax.set_ylabel('theta (m3/m3)')
        ax.plot([0,0.6],[0.0,0.6],'r-')
        for iK in range(nK):
            color=cm(float(iK)/(nK-0.999))
            ax.plot(validationswc ,CompareK_Theta[iK][:,numpy.argmin((z-zplot[iZ])**2)],'.',color=color,label='%s'%K_unique[iK],alpha=0.1)
        ax.legend(loc=2)
        ax.set_aspect('equal')
    fig.suptitle('%s_%s'%(FC_unique[iFC],LayerLabel[iL]))
    figlist2.append(fig)
    
#    pylab.close('all')
    
        
    
#==============================================================================
# Save
#==============================================================================

## Save data
import sys
basename=os.path.basename(sys.argv[0])[:-3]

#postfix='_%s_%s_%s_%s'%(FC,SFC,HC,K)
postfix=''


## save plots to pdf


#pdfname = os.path.join(Root, basename + postfix + '.pdf')
#pp = PdfPages(pdfname)
#for fig in figlist:
#    pp.savefig(fig)
#
#pp.close()

pdfname = os.path.join(Root, basename + postfix + '2.pdf')
pp = PdfPages(pdfname)
for fig in figlist2:
    pp.savefig(fig)

pp.close()

#extension='.png'
#for iF in range(len(figlist)):
#    fig=figlist[iF]
#    figname = os.path.join(Root, basename+'_%s'%iF + postfix + extension)
#    fig.savefig(figname, bbox_inches='tight')

#pylab.show()
os.startfile(pdfname)
