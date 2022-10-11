# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 13:37:36 2015

@author: Administrator
"""

#==============================================================================
# Doublecheck outcome with hydrus + massbalance
#==============================================================================

import os
import pylab
import numpy
import copy 

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

#==============================================================================
# Open Control
#==============================================================================

ProjectFolder='C:\\Users\\Administrator\\Documents\\Arbeid\\Nevenactiviteiten\\Hydrus_Matlab\\test_MakeItRain'

TLevelOutName=os.path.join(ProjectFolder,'T_Level.out')
TlevelData,TLevel_Header=ReadTLevelOut(TLevelOutName)

NodInfOutName=os.path.join(ProjectFolder,'Nod_Inf.out')
NodInf_T,NodInf_Data,NodInf_Header=ReadNodInfOut(NodInfOutName)

ObsOutName=os.path.join(ProjectFolder,'Obs_Node.out')
ObsNode_Nodes,ObsNode_t,ObsNode_L=ReadObsNodeOut(ObsOutName)

t_hydrus=TlevelData[:,0]
dT_Hydrus=numpy.append(t_hydrus[0],t_hydrus[1:]-t_hydrus[:-1])

#==============================================================================
# Initial
#==============================================================================

nz=201
zmin=0
zmax=-100
z=numpy.linspace(zmin,zmax,nz)    # cm // choose z-axis upwards or downwards!
#dz=(z[2:]-z[:-2])/2

t0=0.
t_end=1*24*60*60.                 # s
#t_end=2e4                # s
dt_init=0.1                 # s
dt_min=0.001                # s
dt_max=1e3                 # s

h0=numpy.ones(nz)*-100      # cm
#h0[0]=5                     # cm
#h0[0]=0                     # cm
#h0[-1]=-100                  # cm

#t_rain_L=[5e4,2*24*60*60.]
#q_rain_L=[-4e-4,3e-4] # cm/s
t_rain_L=[5e4,8e4,2*24*60*60.]
q_rain_L=[-4e-4,3e-4,-8e-4] # cm/s

#q_bottom=0  # cm/s

#==============================================================================
# Iteration parameters
#==============================================================================

#theta_tol=0.001     # m3/m3
#h_tol=0.1             # cm
#q_tol=1e-5
#pF_tol=0.1

theta_tol=0.001     # m3/m3
h_tol=0.01             # cm
q_tol=1e5
pF_tol=1e5

#theta_tol_rel=1e-2     # m3/m3
#h_tol_rel=1e-2             # cm
#q_tol_rel=1e-2

nit_decr_dt=5
decr_dt=0.8
nit_incr_dt=2
incr_dt=1.1

nit_term=20
dt_term=0.8

#==============================================================================
# Soil properties
#==============================================================================

#incl=numpy.pi/2      # inclinatie v z-as tov verticale
incl=0
n=numpy.ones(nz)*1.56
#m=numpy.ones(nz)*0.6
m=1-1/n
alpha=numpy.ones(nz)*0.036  # 1/cm
thetar=numpy.ones(nz)*0.078
thetas=numpy.ones(nz)*0.43
l=numpy.ones(nz)*0.5
Ks=numpy.ones(nz)*2.88889e-4   # cm/s

ModelSWRC=h2theta_VanGenuchten5
ParSWRC=[thetar,thetas,alpha,n,m]
ModelKh=h2K_Mualem5
ParKh=[Ks,l,alpha,n,m]
ModelCap=h2C_VanGenuchten5
ParCap=[thetar,thetas,alpha,n,m]


ParSWRC_2D=map(lambda x : x[None,:],ParSWRC)
ParKh_2D=map(lambda x : x[None,:],ParKh)


#==============================================================================
# Iteration
#==============================================================================

t_request=[4000,80000] # list of forced t
t_request=t_request+t_rain_L
t_request=numpy.sort(numpy.unique(t_request))

dt=[]       # List with dt (s)
t=[t0]      # List with t (s)
h=[h0]      # List with h (m)
nit=[]      # List with n it
q_runoff=[]
BC_top=[]
h_SL=[0,]

dt_now=dt_init
end_T=0
iT=0
end_IT=0

BC_switch_prev=0
h_SurfaceLayer=numpy.array([0.,])
h_SurfaceLayer_MAX=10.0
#h_SurfaceLayer_MAX=0.01
h_SurfaceLayer_MIN=-1e5

while end_T==0:
#    if dt_now > dt_max:
#        dt_now=dt_max
    
    h_prevT=copy.deepcopy(h[iT])
    t_now=t[iT]+dt_now
    
    
    ## Adapt T_request
            
    # check t_request
    if numpy.size(t_request)>0:
        if t_now > t_request[0]:
            dt_now=t_request[0]-t[iT]
            t_now=t_request[0]
            t_request=t_request[1:]
    
    # finish when t_end is reached
    if t_now>=t_end:
        dt_now=t_end-t[iT]
        t_now=t_end
        end_T=1
        
    ## q_rain + BC_switch
    
#    q_rain=numpy.array(q_rain_L)[numpy.array(t_rain_L)>=t_now][0]
    q_rain=numpy.array(q_rain_L)[numpy.array(t_rain_L)>t[iT]][0]
    
    # avoid overshoot at top BC
    if (h_SurfaceLayer>0) & (q_rain > 0):
        dt_max_qsurface=h_SurfaceLayer[0]/q_rain
        if dt_max_qsurface<dt_now:
            dt_now=dt_max_qsurface
            t_now=t[iT]+dt_now    
    
    
    if h_SurfaceLayer<=0:
        if h_prevT[0]>=0:
            dz=z[1]-z[0]
            dh=h_prevT[1]-h_prevT[0]
    #        dh=h_prevT[1]-h_SurfaceLayer
            K_PrevT=ModelKh(h_prevT,ParKh)
#            K12=(ParKh[0][0]+K_PrevT[0])/2
            K12=(K_PrevT[1]+K_PrevT[0])/2
    #        K12=K_PrevT[0]
#            K12=ParKh[0][0]
            
            qlim=-K12*(dh/dz+numpy.cos(incl))
            
    #        BC_switch = abs(qlim) >= abs(q_rain) # 0 : Dirichelet // 1 : Neumann
            BC_switch = numpy.array([(qlim <= q_rain),]) # 0 : Dirichelet // 1 : Neumann
    #        print(qlim)
        else:
            if q_rain>0: # upward flux
                BC_switch = numpy.array([(h_prevT[0]>h_SurfaceLayer_MIN),]) # 0 : Dirichelet // 1 : Neumann
            else: # downward or no flux
                BC_switch = numpy.array([True,]) # 0 : Dirichelet // 1 : Neumann
            
    else:
        BC_switch=numpy.array([ False,])
    
    ResetdT=numpy.array([ False,])
    if any(t[iT] == numpy.array(t_rain_L)):
        ResetdT=numpy.array([ True,])
    
    if ((iT>0) & (BC_switch_prev!=BC_switch) & end_IT==1)|(ResetdT):
        dt_now=dt_init
        t_now=t[iT]+dt_now

    # Override everything and use Hydrus time
#    dt_now=dT_Hydrus[iT]    # use Hydrus dt
#    t_now=t[iT]+dt_now     # use Hydrus dt
    

    print('Time: %1.1f / %1.1f - BC: %s - q_rain = %1.2e - h_SL = %1.2e' %(t_now,t_end,['Dirichelet','Neumann'][BC_switch[0]],q_rain,h_SurfaceLayer))
    
    
    ## initialize h_prevIT
    
    if iT > 0:
        dt_prev=dt[iT-1]
        h_prevIT=h[iT]+dt_now/dt_prev*(h[iT]-h[iT-1]) ## Linear extrapolation of previous solution
    else:
        h_prevIT=copy.deepcopy(h[iT])
#    h_prevIT=copy.deepcopy(h[iT])
    
    
    ## Start Iteration    
    
    nit_now=0
    end_IT=0
    
    while (end_IT==0) & (nit_now<=nit_term):
#        print(nit_now)
        
        
        # GET P & F WITHOUT BOUNDARY CONDITIONS
        P=P_j1k(h_prevIT,z,dt_now,incl,ModelSWRC,ParSWRC,ModelKh,ParKh,ModelCap,ParCap)
        F=F_j1k(h_prevIT,h_prevT,z,dt_now,incl,ModelSWRC,ParSWRC,ModelKh,ParKh,ModelCap,ParCap)

        # TOP BOUNDARY CONDITION

        if BC_switch :
#            print('Neumann')
            BC='Neumann'
            P,F=BC_NeumannTop(P,F,h_prevIT,h_prevT,dt_now,q_rain,z,incl,ModelSWRC,ParSWRC,ModelKh,ParKh,ModelCap,ParCap)
        else:
#            print('Dirichelet')
            BC='Dirichelet'
            if (q_rain>0) & (h_SurfaceLayer<=0):
                h_BC=h_SurfaceLayer_MIN
            else:
                h_BC=h_SurfaceLayer
            P,F=BC_DiricheletTop(P,F,h_prevIT,h_prevT,dt_now,h_BC,z,ModelKh,ParKh)
            
        
        # BOTTOM BOUNDARY CONDITION
#        h_BC=h0[-1]
#        P,F=BC_DiricheletBottom(P,F,h_prevIT,h_prevT,dt_now,h_BC,z,ModelKh,ParKh)
    
#        P,F=BC_NeumannBottom(P,F,h_prevIT,h_prevT,dt_now,q_bottom,z,incl,ModelSWRC,ParSWRC,ModelKh,ParKh,ModelCap,ParCap)
        
        P,F=BC_FreeDrainageBottom(P,F,h_prevIT,h_prevT,dt_now,z,incl,ModelSWRC,ParSWRC,ModelKh,ParKh,ModelCap,ParCap)
        
        # CALCULATE H_NOW
        h_now=numpy.dot(numpy.linalg.inv(P),F).squeeze()
#        h_now=gausselim(P,F).squeeze()
        
        
        # CHECK THETA_DIFF & H_DIFF
        nit_now=nit_now+1
        if nit_now > 1:
            theta_diff=abs(ModelSWRC(h_now,ParSWRC)-ModelSWRC(h_prevIT,ParSWRC))
            h_diff=abs(h_now-h_prevIT)
            q_diff=abs(Flux_1D(z,h_now,incl,ModelKh,ParKh)-Flux_1D(z,h_prevIT,incl,ModelKh,ParKh))
            pF_diff=abs(numpy.log10(abs(h_now))-numpy.log10(abs(h_prevIT)))
            pF_diff=pF_diff[abs(h_diff)>0.5]
            pF_diff=pF_diff[~numpy.isnan(pF_diff)]
            
#            theta_diff_rel=abs(ModelSWRC(h_now,ParSWRC)/(1e-20+ModelSWRC(h_prevIT,ParSWRC))-1)
#            h_diff_rel=abs((h_now)/(1e-20+h_prevIT)-1)
#            q_diff_rel=abs(Flux_1D(z,h_now,incl,ModelKh,ParKh)/(1e-20+Flux_1D(z,h_prevIT,incl,ModelKh,ParKh))-1)
#            h_diff_rel=h_diff_rel[abs(h_prevIT)>0.5]
#            [theta_diff_rel,h_diff_rel,q_diff_rel]=map(lambda x : x[~numpy.isnan(x)],[theta_diff_rel,h_diff_rel,q_diff_rel])
            
#            print('theta : %s | h : %s | q : %s | pF : %s '%(all(theta_diff < theta_tol) , all(h_diff < h_tol), all(q_diff < q_tol) , all(pF_diff < pF_tol)))
            
            if all(theta_diff < theta_tol) & all(h_diff < h_tol) & all(q_diff < q_tol) & all(pF_diff < pF_tol):
#            if all(theta_diff < theta_tol) & all(h_diff < h_tol):
#            if all(theta_diff_rel < theta_tol_rel) & all(h_diff_rel < h_tol_rel) & all(q_diff_rel < q_tol_rel):
                end_IT=1
        h_prevIT=copy.deepcopy(h_now)
            
    if end_IT==1:
        # Save h
        h.append(h_now)
        nit.append(nit_now)
        t.append(t_now)
        dt.append(dt_now)
        BC_top.append(BC)
        BC_switch_prev=BC_switch
        
        # update h_surface
        if BC=='Neumann':
            if h_now[0]<h_SurfaceLayer_MIN:
                h_now[0]=h_SurfaceLayer_MIN
            
        if BC=='Dirichelet':
            qtop=Flux_Top(z,numpy.array([h_prevT,h_now]),dt_now,incl,ModelSWRC,ParSWRC,ModelKh,ParKh)
            h_SurfaceLayer=h_SurfaceLayer+(qtop-q_rain)*dt_now
        
        # Runoff
        h_runoff=max([h_SurfaceLayer-h_SurfaceLayer_MAX,0])
        h_SurfaceLayer=h_SurfaceLayer-h_runoff
        h_SurfaceLayer=max([h_SurfaceLayer,numpy.array([0.,])])
        
        h_SL.append(h_SurfaceLayer)
        q_runoff.append((h_runoff)/dt_now)
        
        # New time
        iT+=1
        if (nit_now <= nit_incr_dt) & (dt_now<dt_max):
            print('Increasing time step')
            dt_now=dt_now*incr_dt
        elif (nit_now >= nit_decr_dt) & (dt_now>dt_min):
            print('Decreasing time step')
            dt_now=dt_now*decr_dt
            
    elif(dt_now>dt_min):
        print('Decreasing time step drastically')
        dt_now=dt_now*dt_term
        print dt_now
    else:
        print('Not converging, minimal time step reached, ending calculation')
        end_T=1
    
h=numpy.array(h)
theta=ModelSWRC(h,ParSWRC_2D)
t=numpy.array(t)
dt=numpy.array(dt)
nit=numpy.array(nit)
H=h+(-z*numpy.cos(incl))
q_runoff=numpy.array(q_runoff)
h_SL=numpy.array(h_SL)
BC_top=numpy.array(BC_top)
#dz=z[1]-z[0]
#dh=h[:,1]-h[:,0]
#K12=(ModelKh(h[:,1],map(lambda x : x[None,1],ParKh))+ModelKh(h[:,0],map(lambda x : x[None,0],ParKh)))/2
#K12_T12=(K12[1:]+K12[:-1])/2
#dh_T12=(dh[1:]+dh[:-1])/2
#dV=theta[:,0]*dz/2
#dV_dt12=(dV[:-1]-dV[1:])/dt
#q_top=-dV_dt12-K12_T12*(dh_T12/dz+numpy.cos(incl))
##q_top=-dV_dt12


t_T12=(t[1:]+t[:-1])/2

#fig=pylab.figure()
#ax=fig.add_subplot(111)
#ax.plot(t_T12,q_top,'b.-',label='q_top')
#ax.plot(t_T12,q_runoff,'r.-',label='q_runoff')
#ax.plot(t_T12,q_top-q_runoff,'g.-',label='q_tot')
#ax.legend()

#==============================================================================
# Flux
#==============================================================================

q_top_hydrus=TlevelData[:,3]

q_top=Flux_Top(z,h,dt,incl,ModelSWRC,ParSWRC,ModelKh,ParKh)


q_bottom_hydrus=TlevelData[:,5]

q_bottom=Flux_Bottom(z,h,dt,incl,ModelSWRC,ParSWRC,ModelKh,ParKh)


q_runoff_hydrus=TlevelData[:,14]

#==============================================================================
# Check
#==============================================================================


fig=pylab.figure()
ax=fig.add_subplot(111)
#ax.imshow(h.transpose(),vmin=-100,vmax=5,interpolation='nearest',cmap='viridis')
im=ax.contourf(t,z,h.transpose(),numpy.linspace(-120,10,11),cmap='viridis')
pylab.rcParams['contour.negative_linestyle'] = 'solid'
ax.contour(t,z,h.transpose(),[0,],vmin=-100,vmax=10,colors='k')
ax.axis('tight')
ax.set_title('Tensiometric Head')
fig.colorbar(im)
figlist.append(fig)

fig=pylab.figure()
ax=fig.add_subplot(111)
#ax.imshow(h.transpose(),vmin=-100,vmax=5,interpolation='nearest',cmap='viridis')
im=ax.contourf(t,z,theta.transpose(),20,cmap='viridis')
pylab.rcParams['contour.negative_linestyle'] = 'solid'
ax.contour(t,z,theta.transpose(),[max(thetas)-0.01,],colors='k')
ax.axis('tight')
ax.set_title('Water Content')
fig.colorbar(im)
figlist.append(fig)

#fig=pylab.figure()
#ax=fig.add_subplot(111)
##ax.imshow(h.transpose(),vmin=-100,vmax=5,interpolation='nearest',cmap='viridis')
#im=ax.contourf(t,z,h.transpose(),numpy.linspace(-120,10,11),cmap='viridis')
#pylab.rcParams['contour.negative_linestyle'] = 'solid'
#ax.contour(t,z,h.transpose(),[0,],vmin=-100,vmax=10,colors='k')
#ax.axis('tight')
#fig.colorbar(im)
#figlist.append(fig)


fig=pylab.figure()
ax=fig.add_subplot(111)
ax.plot(nit,'k.-',label='Replicate')
ax.set_ylim([0,None])
ax.set_title('number of iterations')
figlist.append(fig)

fig=pylab.figure()
ax=fig.add_subplot(111)
ax.plot(dt,'k.-',label='Replicate')
ax.set_yscale('log')
ax.set_title('dt')
figlist.append(fig)

fig=pylab.figure()
ax=fig.add_subplot(111)
ax.plot(t[1:],q_top,'k.-',label='Replicate')
ax.plot(t_hydrus,q_top_hydrus,'r.-',label='Hydrus')
ax.set_xlabel(TLevel_Header[0])
ax.set_ylabel(TLevel_Header[3])
ax.legend()
ax.set_title(TLevel_Header[3])
figlist.append(fig)

fig=pylab.figure()
ax=fig.add_subplot(111)
ax.plot(t[1:],q_bottom,'k.-',label='Replicate')
ax.plot(t_hydrus,q_bottom_hydrus,'r.-',label='Hydrus')
ax.set_xlabel(TLevel_Header[0])
ax.set_ylabel(TLevel_Header[5])
ax.legend()
ax.set_title(TLevel_Header[5])
figlist.append(fig)


fig=pylab.figure()
ax=fig.add_subplot(111)
ax.plot(t[1:],q_runoff,'k.-',label='Replicate')
ax.plot(t_hydrus,q_runoff_hydrus,'r.-',label='Hydrus')
ax.set_xlabel(TLevel_Header[0])
ax.set_ylabel(TLevel_Header[14])
ax.legend()
ax.set_title(TLevel_Header[14])
figlist.append(fig)



#nP_Hydrus=NodInf_T.size
#for iP in range(nP_Hydrus):
#    fig=pylab.figure()
#    ax=fig.add_subplot(111)
#    ax.plot(h[t==NodInf_T[iP],:][0,:],z,'k.-',label='Replicate')
#    ax.plot(NodInf_Data[iP][:,2],NodInf_Data[iP][:,1],'r.-',label='Hydrus')
#    ax.set_ylabel(NodInf_Header[0])
#    ax.set_xlabel(NodInf_Header[2])
#    ax.set_title('t = %s s'%NodInf_T[iP])
#    ax.legend()
#
#
#nP_Hydrus=NodInf_T.size
#for iP in range(nP_Hydrus):
#    fig=pylab.figure()
#    ax=fig.add_subplot(111)
#    ax.plot(theta[t==NodInf_T[iP],:][0,:],z,'k.-',label='Replicate')
#    ax.plot(NodInf_Data[iP][:,3],NodInf_Data[iP][:,1],'r.-',label='Hydrus')
#    ax.set_ylabel(NodInf_Header[1])
#    ax.set_xlabel(NodInf_Header[3])
#    ax.set_title('t = %s s'%NodInf_T[iP])
#    ax.legend()

nO=ObsNode_Nodes.size
for iO in range(nO):
    fig=pylab.figure(figsize=figsize1)
    ax=fig.add_subplot(111)
    iZ=ObsNode_Nodes[iO]-1
    ax.plot(t,h[:,iZ],'k-',label='Replicate')
    ax.plot(ObsNode_t,ObsNode_L[iO][:,0],'r-',label='Hydrus')
    ax.set_xlabel('Time')
    ax.set_ylabel('Head')
    ax.legend(loc=4)
    #ax.set_xlim([0,.1])
    #ax.set_ylim([-1000,0])
    ax.set_title('Observation node z = %s cm' %z[iZ])
    figlist.append(fig)


fig=pylab.figure()
ax=fig.add_subplot(111)
ax.axhline(0,color='k')
ax.plot(t[1:][BC_top=='Dirichelet'],h_SL[:-1][BC_top=='Dirichelet'],'r.',label='Dirichelet')
ax.plot(t[1:][BC_top=='Neumann'],h_SL[:-1][BC_top=='Neumann'],'b.',label='Neumann')
ax.plot(t,h[:,0],'k--')
ax.set_xlabel('h Surface Layer (cm)')
ax.set_ylabel('Head')
ax.set_ylim([-0.1,10])
ax.legend()
figlist.append(fig)


## save plots to pdf

import sys
basename=os.path.basename(sys.argv[0])[:-3]

Root='C:\Users\Administrator\Documents\Arbeid\Figuurtjes'
postfix=''

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

#pylab.show()
os.startfile(pdfname)
