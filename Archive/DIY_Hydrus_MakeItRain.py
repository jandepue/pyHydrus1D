# -*- coding: utf-8 -*-
"""
Created on Tue Oct 27 15:50:37 2015

@author: Administrator
"""

#==============================================================================
# Hydrus Template
#==============================================================================

import os
import pylab
import numpy
import copy

# from Hydrus_Funky import *
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
# Initial
#==============================================================================

nz=201
zmin=0
zmax=-100
z=numpy.linspace(zmin,zmax,nz)    # cm // choose z-axis upwards or downwards!
#dz=(z[2:]-z[:-2])/2

t0=0.
t_end=1*24*60*60.                 # s
t_end=2e4                # s
dt_init=0.1                 # s
dt_min=0.001                # s
dt_max=1e2                 # s

h0=numpy.ones(nz)*-100      # cm
#h0[0]=5                     # cm
#h0[0]=0                     # cm
#h0[-1]=-100                  # cm

q_rain=-4e-4 # cm/s

#q_bottom=0  # cm/s

#==============================================================================
# Iteration parameters
#==============================================================================

theta_tol=0.0001     # m3/m3
h_tol=0.1             # cm
q_tol=1e-5
pF_tol=0.1

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
t_request=numpy.sort(t_request)

dt=[]       # List with dt (s)
t=[t0]      # List with t (s)
h=[h0]      # List with h (m)
nit=[]      # List with n it
q_runoff=[]
BC_top=[]

dt_now=dt_init
end_T=0
iT=0

BC_switch_prev=0
h_SurfaceLayer=0
h_SurfaceLayer_MAX=10

while end_T==0:
#    if dt_now > dt_max:
#        dt_now=dt_max

    h_prevT=copy.deepcopy(h[iT])

    t_now=t[iT]+dt_now

    dz=z[1]-z[0]
#    dh=h_prevT[1]-0
    dh=h_prevT[1]-h_SurfaceLayer
    K_PrevT=ModelKh(h_prevT,ParKh)
    K12=(ParKh[0][0]+K_PrevT[1])/2
    qlim=-K12*(dh/dz+numpy.cos(incl))

    BC_switch = abs(qlim) >= abs(q_rain) # 0 : Dirichelet // 1 : Neumann

    if (iT>0) & (BC_switch_prev!=BC_switch):
        dt_now=dt_init

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

    print('Time: %1.2f / %1.2f ' %(t_now,t_end))

    nit_now=0
    if iT > 0:
        dt_prev=dt[iT-1]
        h_prevIT=h[iT]+dt_now/dt_prev*(h[iT]-h[iT-1])
    else:
        h_prevIT=copy.deepcopy(h[iT])
#    h_prevIT=copy.deepcopy(h[iT])

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
            # P,F=BC_NeumannTop2(P,F,h_prevIT,h_prevT,dt_now,q_rain,z,incl,ModelSWRC,ParSWRC,ModelKh,ParKh,ModelCap,ParCap)
        else:
#            print('Dirichelet')
            BC='Dirichelet'
            h_BC=h_SurfaceLayer
            P,F=BC_DiricheletTop(P,F,h_prevIT,h_prevT,dt_now,h_BC,z,ModelKh,ParKh)


        # BOTTOM BOUNDARY CONDITION
#        h_BC=h0[-1]
#        P,F=BC_DiricheletBottom(P,F,h_prevIT,h_prevT,dt_now,h_BC,z,ModelKh,ParKh)

#        P,F=BC_NeumannBottom(P,F,h_prevIT,h_prevT,dt_now,q_bottom,z,incl,ModelSWRC,ParSWRC,ModelKh,ParKh,ModelCap,ParCap)

        P,F=BC_FreeDrainageBottom(P,F,h_prevIT,h_prevT,dt_now,z,incl,ModelSWRC,ParSWRC,ModelKh,ParKh,ModelCap,ParCap)

        # CALCULATE H_NOW
        h_now=numpy.dot(numpy.linalg.inv(P),F).squeeze()


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
        if BC=='Dirichelet':
            qtop=Flux_Top(z,numpy.array([h_prevT,h_now]),dt_now,incl,ModelSWRC,ParSWRC,ModelKh,ParKh)
            h_SurfaceLayer=h_SurfaceLayer+(qtop-q_rain)*dt_now

        # Runoff
        h_runoff=max([h_SurfaceLayer-h_SurfaceLayer_MAX,0])
        h_SurfaceLayer=h_SurfaceLayer-h_runoff
        q_runoff.append((h_runoff)/dt_now)

        # New time
        iT+=1
        if iT<5:
            ddt=1
        else:
            ddt=all(numpy.array(dt)[-4:]-numpy.array(dt)[-5:-1] >=0)

        if (nit_now <= nit_incr_dt) & (dt_now<dt_max) & (ddt==1):
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
theta=ModelSWRC(h,ParSWRC_2D)
t=numpy.array(t)
dt=numpy.array(dt)
H=h+(-z*numpy.cos(incl))
q_runoff=numpy.array(q_runoff)


#dz=z[1]-z[0]
#dh=h[:,1]-h[:,0]
#K12=(ModelKh(h[:,1],map(lambda x : x[None,1],ParKh))+ModelKh(h[:,0],map(lambda x : x[None,0],ParKh)))/2
#K12_T12=(K12[1:]+K12[:-1])/2
#dh_T12=(dh[1:]+dh[:-1])/2
#dV=theta[:,0]*dz/2
#dV_dt12=(dV[:-1]-dV[1:])/dt
#q_top=-dV_dt12-K12_T12*(dh_T12/dz+numpy.cos(incl))
##q_top=-dV_dt12

q_top=Flux_Top(z,h,dt,incl,ModelSWRC,ParSWRC,ModelKh,ParKh)

t_T12=(t[1:]+t[:-1])/2

fig=pylab.figure()
ax=fig.add_subplot(111)
ax.plot(t_T12,q_top,'b.-',label='q_top')
ax.plot(t_T12,q_runoff,'r.-',label='q_runoff')
ax.plot(t_T12,q_top-q_runoff,'g.-',label='q_tot')
ax.legend()

fig=pylab.figure()
ax=fig.add_subplot(111)
ax.plot(t,h[:,0],'.-')
ax.plot(t,h[:,1],'.-')
ax.plot(t,h[:,2],'.-')


#fig=pylab.figure()
#ax=fig.add_subplot(111)
#ax.plot(h[-1,:],z,'.-')

pylab.show()
