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

ProjectFolder='C:\\Users\\Administrator\\Documents\\Arbeid\\Nevenactiviteiten\\Hydrus_Matlab\\test'

TLevelOutName=os.path.join(ProjectFolder,'T_Level.out')
TlevelData,TLevel_Header=ReadTLevelOut(TLevelOutName)

NodInfOutName=os.path.join(ProjectFolder,'Nod_Inf.out')
NodInf_T,NodInf_Data,NodInf_Header=ReadNodInfOut(NodInfOutName)

dT_Hydrus=numpy.append(TlevelData[0,0],TlevelData[1:,0]-TlevelData[:-1,0])
dV_Hydrus=numpy.append(0,TlevelData[1:,16]-TlevelData[:-1,16])

Balance_Hydrus=dT_Hydrus*TlevelData[:,3]+dT_Hydrus*TlevelData[:,5]-dV_Hydrus

#==============================================================================
# Initial
#==============================================================================

nz=101
zmin=0
zmax=-100
z=numpy.linspace(zmin,zmax,nz)    # cm // choose z-axis upwards or downwards!
dz=(z[2:]-z[:-2])/2

t0=0.
t_end=1*24*60*60.                 # s
dt_init=0.1                 # s
dt_min=0.001                # s
dt_max=1e5                   # s

h0=numpy.ones(nz)*-100.0     # cm
#h0[0]=-5                     # cm
#h0[-1]=-100                  # cm

q_bottom=0  # cm/s
q_top=-2e-4 # cm/s

#==============================================================================
# Iteration parameters
#==============================================================================

theta_tol=0.001     # m3/m3
h_tol=0.1         # cm

nit_decr_dt=7
decr_dt=0.8
nit_incr_dt=3
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

#==============================================================================
# Iteration
#==============================================================================

t_request=[4000,80000] # list of forced t

t_request=numpy.sort(t_request)

dt=[]       # List with dt (s)
t=[t0]      # List with t (s)
h=[h0]      # List with h (m)
nit=[]      # List with n it

dt_now=dt_init
end_T=0
iT=0

while end_T==0:
    t_now=t[iT]+dt_now
    
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
    

#    dt_now=dT_Hydrus[iT]    
#    t_now=t[iT]+dt_now

    print('Time: %s / %s ' %(t_now,t_end))
        
    nit_now=0
    h_prevT=copy.deepcopy(h[iT])
    if iT > 0:
        dt_prev=dt[iT-1]
        h_prevIT=h[iT]+dt_now/dt_prev*(h[iT]-h[iT-1])
    else:
        h_prevIT=copy.deepcopy(h[iT])

    end_IT=0
    
    while (end_IT==0) & (nit_now<=nit_term):
#        print(nit_now)
        
        # GET P & F WITHOUT BOUNDARY CONDITIONS
        P=P_j1k(h_prevIT,z,dt_now,incl,ModelSWRC,ParSWRC,ModelKh,ParKh,ModelCap,ParCap)
        F=F_j1k(h_prevIT,h_prevT,z,dt_now,incl,ModelSWRC,ParSWRC,ModelKh,ParKh,ModelCap,ParCap)

        # TOP BOUNDARY CONDITION
#        h_BC=h0[0]
#        P,F=BC_DiricheletTop(P,F,h_prevIT,h_prevT,dt_now,h_BC,z,ModelKh,ParKh)
        
        q_BC=q_top
        P,F=BC_NeumannTop(P,F,h_prevIT,h_prevT,dt_now,q_BC,z,incl,ModelSWRC,ParSWRC,ModelKh,ParKh,ModelCap,ParCap)
        
        # BOTTOM BOUNDARY CONDITION
#        h_BC=h0[-1]
#        P,F=BC_DiricheletBottom(P,F,h_prevIT,h_prevT,dt_now,h_BC,z,ModelKh,ParKh)
    
#        P,F=BC_NeumannBottom(P,F,h_prevIT,h_prevT,dt_now,q_BC,z,incl,ModelSWRC,ParSWRC,ModelKh,ParKh,ModelCap,ParCap)
        P,F=BC_FreeDrainageBottom(P,F,h_prevIT,h_prevT,dt_now,z,incl,ModelSWRC,ParSWRC,ModelKh,ParKh,ModelCap,ParCap)
        
        # CALCULATE H_NOW
#        h_now=numpy.dot(numpy.linalg.inv(P),F).squeeze()
#        h_now=numpy.linalg.solve(P,F).squeeze()
        h_now=gausselim(P,F).squeeze()
        
        # CHECK THETA_DIFF & H_DIFF
        nit_now=nit_now+1
        if nit_now > 1:
            theta_diff=abs(ModelSWRC(h_now,ParSWRC)-ModelSWRC(h_prevIT,ParSWRC));
            h_diff=abs(h_now-h_prevIT);
            if all(theta_diff < theta_tol) & all(h_diff < h_tol):
                end_IT=1
        h_prevIT=copy.deepcopy(h_now)
            
    if end_IT==1:
        # Save h
        h.append(h_now)
        nit.append(nit_now)
        t.append(t_now)
        dt.append(dt_now)
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
        
    else:
        print('Not converging, minimal time step reached, ending calculation')
        end_T=1

h=numpy.array(h)
t=numpy.array(t)
dt=numpy.array(dt)



#TopBC_Funk=BC_NeumannBottom
#q_BC=-1e-5
#TopBC_Par=(q_BC,z,incl,ModelSWRC,ParSWRC,ModelKh,ParKh,ModelCap,ParCap)
#
##TopBC_Funk=BC_DiricheletTop
##h_BC=h0[0]
##TopBC_Par=(h_BC,z,ModelKh,ParKh)
#
#BottomBC_Funk=BC_DiricheletBottom
#h_BC=h0[-1]
#BottomBC_Par=(h_BC,z,ModelKh,ParKh)
#
#t,dt,h=DoHydrus(h0,z,t0,t_end,dt_init,incl,ModelSWRC,ParSWRC,ModelKh,ParKh,ModelCap,ParCap,TopBC_Funk,TopBC_Par,BottomBC_Funk,BottomBC_Par,h_tol,theta_tol,nit_term,dt_min,dt_max,nit_decr_dt,nit_incr_dt,decr_dt,incr_dt,dt_term,t_request=t_request)

ParSWRC_2D=map(lambda x : x[None,:],ParSWRC)
ParKh_2D=map(lambda x : x[None,:],ParKh)

theta=ModelSWRC(h,ParSWRC_2D)
H=h+(-z*numpy.cos(incl))

V=Storage(theta,z)
V_tot=V.sum(axis=1)
dV_tot=V_tot[1:]-V_tot[:-1]

dz_top=(z[1]-z[0])
K_prof=ModelKh(h,ParKh_2D)

K_top=(K_prof[1:,0]+K_prof[1:,1]+K_prof[:-1,0]+K_prof[:-1,1])/4.0
flux_top=-K_top*(((h[:-1,1]-h[:-1,0])+(h[1:,1]-h[1:,0]))/2+dz_top*numpy.cos(incl))/dz_top
Balance_Replicate=flux_top*dt-dV_tot

#K_top=(K_prof[:,0]+K_prof[:,1])/2.0
#flux_top=-K_top*((h[:,1]-h[:,0])+dz_top*numpy.cos(incl))/dz_top
#flux_between=(flux_top[1:]+flux_top[:-1])/2.0
#Balance_Replicate=flux_between*dt-dV_tot



#==============================================================================
# vTop
#==============================================================================

vtop_hydrus=TlevelData[:,3]

dz12=z[1]-z[0]
dh12=h[:,1]-h[:,0]
K12=(ModelKh(h[:,1],map(lambda x : x[None,1],ParKh))+ModelKh(h[:,0],map(lambda x : x[None,0],ParKh)))/2
q12=(-K12*(dh12/dz12+numpy.cos(incl)))[1:]
dV=(theta[1:,0] - theta[:-1,0])

q_top=(dV*dz12)/(dt*2.0)+q12


vbottom_hydrus=TlevelData[:,5]

dz12=z[-1]-z[-2]
dh12=h[:,-1]-h[:,-2]
K12=(ModelKh(h[:,-2],map(lambda x : x[None,-2],ParKh))+ModelKh(h[:,-1],map(lambda x : x[None,-1],ParKh)))/2
q12=(-K12*(dh12/dz12+numpy.cos(incl)))[1:]
dV=(theta[1:,-1] - theta[:-1,-1])

q_bottom=-(dV*dz12)/(dt*2.0)+q12


#==============================================================================
# Check
#==============================================================================

fig=pylab.figure()
ax=fig.add_subplot(111)
ax.plot(t[1:],q_top,'k.-',label='Replicate')
ax.plot(TlevelData[:,0],vtop_hydrus,'r.-',label='Hydrus')
ax.set_xlabel(TLevel_Header[0])
ax.set_ylabel(TLevel_Header[3])

fig=pylab.figure()
ax=fig.add_subplot(111)
ax.plot(t[1:],q_bottom,'k.-',label='Replicate')
ax.plot(TlevelData[:,0],vbottom_hydrus,'r.-',label='Hydrus')
ax.set_xlabel(TLevel_Header[0])
ax.set_ylabel(TLevel_Header[5])


fig=pylab.figure()
ax=fig.add_subplot(111)
ax.plot(t,V_tot,'k.-',label='Replicate')
ax.plot(TlevelData[:,0],TlevelData[:,16],'r.-',label='Hydrus')
ax.set_xlabel(TLevel_Header[0])
ax.set_ylabel(TLevel_Header[16])


nP_Hydrus=NodInf_T.size
for iP in range(nP_Hydrus):
    fig=pylab.figure()
    ax=fig.add_subplot(111)
    ax.plot(h[t==NodInf_T[iP],:][0,:],z,'k.-',label='Replicate')
    ax.plot(NodInf_Data[iP][:,2],NodInf_Data[iP][:,1],'r.-',label='Hydrus')
    ax.set_ylabel(NodInf_Header[0])
    ax.set_xlabel(NodInf_Header[2])
    ax.set_title('t = %s s'%NodInf_T[iP])
    ax.legend()


nP_Hydrus=NodInf_T.size
for iP in range(nP_Hydrus):
    fig=pylab.figure()
    ax=fig.add_subplot(111)
    ax.plot(theta[t==NodInf_T[iP],:][0,:],z,'k.-',label='Replicate')
    ax.plot(NodInf_Data[iP][:,3],NodInf_Data[iP][:,1],'r.-',label='Hydrus')
    ax.set_ylabel(NodInf_Header[1])
    ax.set_xlabel(NodInf_Header[3])
    ax.set_title('t = %s s'%NodInf_T[iP])
    ax.legend()


fig=pylab.figure()
ax=fig.add_subplot(111)
ax.plot(t[1:],Balance_Replicate,'k.-',label='Replicate')
ax.plot(TlevelData[:,0],Balance_Hydrus,'r.-',label='Hydrus')
ax.set_xlabel(TLevel_Header[0])
#ax.set_ylabel(TLevel_Header[16])
ax.legend()


pylab.show()