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

nz=101
zmin=0
zmax=-100
z=numpy.linspace(zmin,zmax,nz)    # cm // choose z-axis upwards or downwards!
dz=(z[2:]-z[:-2])/2

t0=0.
t_end=1*24*60*60.                 # s
dt_init=0.1                 # s
dt_min=0.001                # s
dt_max=1e4                   # s

h0=numpy.ones(nz)*-100      # cm
h0[0]=5                     # cm
#h0[-1]=-100                  # cm

q_bottom=0  # cm/s

#==============================================================================
# Iteration parameters
#==============================================================================

theta_tol=0.001     # m3/m3
h_tol=1             # cm

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
        h_BC=h0[0]
        P,F=BC_DiricheletTop(P,F,h_prevIT,h_prevT,dt_now,h_BC,z,ModelKh,ParKh)

        # BOTTOM BOUNDARY CONDITION
#        h_BC=h0[-1]
#        P,F=BC_DiricheletBottom(P,F,h_prevIT,h_prevT,dt_now,h_BC,z,ModelKh,ParKh)

        P,F=BC_NeumannBottom(P,F,h_prevIT,h_prevT,dt_now,q_bottom,z,incl,ModelSWRC,ParSWRC,ModelKh,ParKh,ModelCap,ParCap)

        # CALCULATE H_NOW
        h_now=numpy.dot(numpy.linalg.inv(P),F).squeeze()

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

# ParSWRC_2D=map(lambda x : x[None,:],ParSWRC)
# ParKh_2D=map(lambda x : x[None,:],ParKh)
ParKh_2D=[x[None,:] for x in ParKh]
ParSWRC_2D=[x[None,:] for x in ParSWRC]

h=numpy.array(h)
theta=ModelSWRC(h,ParSWRC_2D)
t=numpy.array(t)
dt=numpy.array(dt)
H=h+(-z*numpy.cos(incl))
