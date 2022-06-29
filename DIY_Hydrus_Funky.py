# -*- coding: utf-8 -*-
"""
Created on Tue Aug 27 11:07:42 2013

@author: jan
"""

import numpy
import pylab
import scipy.linalg as linalg
import copy

from RetentionConductivityCapacity_Funky import *

#==============================================================================
# Full Hydrus
#==============================================================================

def DoHydrus(h0,z,t0,t_end,dt_init,incl,ModelSWRC,ParSWRC,ModelKh,ParKh,ModelCap,ParCap,TopBC_Funk,TopBC_Par,BottomBC_Funk,BottomBC_Par,h_tol,theta_tol,nit_term,dt_min,dt_max,nit_decr_dt,nit_incr_dt,decr_dt,incr_dt,dt_term,t_request=[]):
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
#            print(nit_now)
            
            # GET P & F WITHOUT BOUNDARY CONDITIONS
            P=P_j1k(h_prevIT,z,dt_now,incl,ModelSWRC,ParSWRC,ModelKh,ParKh,ModelCap,ParCap)
            F=F_j1k(h_prevIT,h_prevT,z,dt_now,incl,ModelSWRC,ParSWRC,ModelKh,ParKh,ModelCap,ParCap)
    
            # TOP BOUNDARY CONDITION
            P,F=TopBC_Funk(P,F,h_prevIT,h_prevT,dt_now,*TopBC_Par)
#            h_BC=h0[0]
#            P,F=BC_DiricheletTop(P,F,h_prevIT,h_BC,z,ModelKh,ParKh)
#            P,F=BC_NeumannTop(P,F,h_prevIT,h_prevT,dt_now,q_top,z,incl,ModelSWRC,ParSWRC,ModelKh,ParKh,ModelCap,ParCap)
            
            # BOTTOM BOUNDARY CONDITION
            P,F=BottomBC_Funk(P,F,h_prevIT,h_prevT,dt_now,*BottomBC_Par)
#            h_BC=h0[-1]
#            P,F=BC_DiricheletBottom(P,F,h_prevIT,h_BC,z,ModelKh,ParKh)
#            P,F=BC_NeumannBottom(P,F,h_prevIT,h_prevT,dt_now,q_bottom,z,incl,ModelSWRC,ParSWRC,ModelKh,ParKh,ModelCap,ParCap)
            
            # CALCULATE H_NOW
#            h_now=numpy.dot(numpy.linalg.inv(P),F).squeeze()
#            h_now=numpy.linalg.solve(P,F).squeeze()
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
            if (nit_now < nit_incr_dt) & (dt_now<dt_max):
                print('Increasing time step')
                dt_now=dt_now*incr_dt
            elif (nit_now > nit_decr_dt) & (dt_now>dt_min):
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

    return t,dt,h
    
#==============================================================================
# Linear algebra
#==============================================================================

#import numpy as np
#from scipy.linalg import lu, inv

def gausselim(A,B):
    """
    https://gist.github.com/mikofski/11192605
    Solve Ax = B using Gaussian elimination and LU decomposition.
    A = LU   decompose A into lower and upper triangular matrices
    LUx = B  substitute into original equation for A
    Let y = Ux and solve:
    Ly = B --> y = (L^-1)B  solve for y using "forward" substitution
    Ux = y --> x = (U^-1)y  solve for x using "backward" substitution
    :param A: coefficients in Ax = B
    :type A: numpy.ndarray of size (m, n)
    :param B: dependent variable in Ax = B
    :type B: numpy.ndarray of size (m, 1)
    """
    # LU decomposition with pivot
    pl, u = linalg.lu(A, permute_l=True)
    # forward substitution to solve for Ly = B
    y = numpy.zeros(B.size)
    for m, b in enumerate(B.flatten()):
        y[m] = b
        # skip for loop if m == 0
        if m:
            for n in xrange(m):
                y[m] -= y[n] * pl[m,n]
        y[m] /= pl[m, m]

    # backward substitution to solve for y = Ux
    x = numpy.zeros(B.size)
    lastidx = B.size - 1  # last index
    for midx in xrange(B.size):
        m = B.size - 1 - midx  # backwards index
        x[m] = y[m]
        if midx:
            for nidx in xrange(midx):
                n = B.size - 1  - nidx
                x[m] -= x[n] * u[m,n]
        x[m] /= u[m, m]
    return x

#==============================================================================
# Funky
#==============================================================================

def F_j1k(h_prevIT,h_prevT,z,dt,incl,ModelSWRC,ParSWRC,ModelKh,ParKh,ModelCap,ParCap):
    import numpy
    
    nz=numpy.shape(z)[0]
    dz=(z[2:]-z[:-2])/2
    
    Cap_prevIT=ModelCap(h_prevIT,ParCap)
    theta_prevIT=ModelSWRC(h_prevIT,ParSWRC)
    theta_prevT=ModelSWRC(h_prevT,ParSWRC)
    K_PrevIT=ModelKh(h_prevIT,ParKh)
    K_Centr=(K_PrevIT[2:]-K_PrevIT[:-2])/2
    
    Cap_prevIT=Cap_prevIT[1:-1]
    h_prevIT=h_prevIT[1:-1]
    theta_prevIT=theta_prevIT[1:-1]
    theta_prevT=theta_prevT[1:-1]
    
    F=numpy.zeros((nz,1))
    F[1:-1,0]=dz/dt*Cap_prevIT*h_prevIT-dz/dt*(theta_prevIT-theta_prevT)+K_Centr*numpy.cos(incl) # WITHOUT BOUNDARY CONDITIONS
    
    return F
    
    
def P_j1k(h_prevIT,z,dt,incl,ModelSWRC,ParSWRC,ModelKh,ParKh,ModelCap,ParCap):
    nz=numpy.shape(z)[0]
    dz=(z[2:]-z[:-2])/2
    dzbw=z[1:-1]-z[:-2]
    dzfw=z[2:]-z[1:-1]
    
    Cap_prevIT=ModelCap(h_prevIT,ParCap)
    K_PrevIT=ModelKh(h_prevIT,ParKh)
    
    Cap_prevIT=Cap_prevIT[1:-1]
    K_bw=(K_PrevIT[1:-1]+K_PrevIT[:-2])/(2*dzbw)
    K_fw=(K_PrevIT[2:]+K_PrevIT[1:-1])/(2*dzfw)
    
    D=dz/dt*Cap_prevIT+K_bw+K_fw
    E=-K_fw[:-1]
#    E=-K_bw[1:]
    
    D=numpy.diag(D,0)
    E1=numpy.diag(E,1)
    E2=numpy.diag(E,-1)
    P= numpy.zeros((nz,nz))
    P[1:-1,1:-1]=D+E1+E2
    
    return P

#==============================================================================
# Boundary Conditions
#==============================================================================

## Dirichelet

def BC_DiricheletTop(P,F,h_prevIT,h_prevT,dt_now,h_BC,z,ModelKh,ParKh):
    P[0,0]=1.0
    K_PrevIT=ModelKh(h_prevIT,ParKh)
    F[1]=F[1]+h_BC*(K_PrevIT[0]+K_PrevIT[1])/(2*(z[1]-z[0]))
    F[0]=h_BC
    
    return P,F


def BC_DiricheletBottom(P,F,h_prevIT,h_prevT,dt_now,h_BC,z,ModelKh,ParKh):
    P[-1,-1]=1.0
    K_PrevIT=ModelKh(h_prevIT,ParKh)
    F[-2]=F[-2]+h_BC*(K_PrevIT[-2]+K_PrevIT[-1])/(2*(z[-1]-z[-2]))
    F[-1]=h_BC

    return P,F


## Neumann

def BC_NeumannTop(P,F,h_prevIT,h_prevT,dt_now,q_BC,z,incl,ModelSWRC,ParSWRC,ModelKh,ParKh,ModelCap,ParCap):
#    q_BC=-0.01; # cm/s !DIRECTION Z-AXIS!
#    q_BC=0.0002 # cm/s !DIRECTION Z-AXIS!
    dz=z[1]-z[0]
    K_PrevIT=ModelKh(h_prevIT,ParKh)
    K12=(K_PrevIT[0]+K_PrevIT[1])/2
    theta_prevIT=ModelSWRC(h_prevIT,ParSWRC)[0]
    theta_prevT=ModelSWRC(h_prevT,ParSWRC)[0]
    Cap_prevIT=ModelCap(h_prevIT,ParCap)[0]

    P[0,0] = dz/(2*dt_now) * Cap_prevIT + K12/dz
    P[0,1]= -K12/dz
    P[1,0]= -K12/dz
    F[0] = dz/(2*dt_now) * Cap_prevIT * h_prevIT[0] - dz/(2*dt_now)*(theta_prevIT-theta_prevT) + K12*numpy.cos(incl) + q_BC

    return P,F
    

def BC_NeumannBottom(P,F,h_prevIT,h_prevT,dt_now,q_BC,z,incl,ModelSWRC,ParSWRC,ModelKh,ParKh,ModelCap,ParCap):
#    q_BC=-0.01; # cm/s !DIRECTION Z-AXIS!
#    q_BC=0.0002 # cm/s !DIRECTION Z-AXIS!
    dz=z[-1]-z[-2]
    K_PrevIT=ModelKh(h_prevIT,ParKh)
    K12=(K_PrevIT[-2]+K_PrevIT[-1])/2;
    theta_prevIT=ModelSWRC(h_prevIT,ParSWRC)[-1]
    theta_prevT=ModelSWRC(h_prevT,ParSWRC)[-1]
    Cap_prevIT=ModelCap(h_prevIT,ParCap)[-1]

    P[-1,-1] = dz/(2*dt_now) * Cap_prevIT + K12/dz
    P[-1,-2]= -K12/dz
    P[-2,-1]= -K12/dz
    F[-1] = dz/(2*dt_now) * Cap_prevIT * h_prevIT[-1] - dz/(2*dt_now)*(theta_prevIT-theta_prevT) - K12*numpy.cos(incl) - q_BC
    
    return P,F


#def BC_NeumannTop2(P,F,h_prevIT,h_prevT,dt_now,q_BC,z,incl,ModelSWRC,ParSWRC,ModelKh,ParKh,ModelCap,ParCap):
##    q_BC=-0.01; # cm/s !DIRECTION Z-AXIS!
##    q_BC=0.0002 # cm/s !DIRECTION Z-AXIS!
#    dz=z[1]-z[0]
#    K_PrevIT=ModelKh(h_prevIT,ParKh)
#    K12=(K_PrevIT[0]+K_PrevIT[1])/2
#    theta_prevIT=ModelSWRC(h_prevIT,ParSWRC)[0]
#    theta_prevT=ModelSWRC(h_prevT,ParSWRC)[0]
#    Cap_prevIT=ModelCap(h_prevIT,ParCap)[0]
#
#    P[0,0] = K12/dz
#    P[0,1]= -K12/dz
#    P[1,0]= -K12/dz
#    F[0] = -K12*numpy.cos(incl) + q_BC
#
#    return P,F
#
#def BC_NeumannBottom2(P,F,h_prevIT,h_prevT,dt_now,q_BC,z,incl,ModelSWRC,ParSWRC,ModelKh,ParKh,ModelCap,ParCap):
##    q_BC=-0.01; # cm/s !DIRECTION Z-AXIS!
##    q_BC=0.0002 # cm/s !DIRECTION Z-AXIS!
#    dz=z[-1]-z[-2]
#    K_PrevIT=ModelKh(h_prevIT,ParKh)
#    K12=(K_PrevIT[-2]+K_PrevIT[-1])/2;
#    theta_prevIT=ModelSWRC(h_prevIT,ParSWRC)[-1]
#    theta_prevT=ModelSWRC(h_prevT,ParSWRC)[-1]
#    Cap_prevIT=ModelCap(h_prevIT,ParCap)[-1]
#
#    P[-1,-1] = K12/dz
#    P[-1,-2]= -K12/dz
#    P[-2,-1]= -K12/dz
#    F[-1] = -K12*numpy.cos(incl) + q_BC
#    
#    return P,F

## FREE DRAINAGE

def BC_FreeDrainageBottom(P,F,h_prevIT,h_prevT,dt_now,z,incl,ModelSWRC,ParSWRC,ModelKh,ParKh,ModelCap,ParCap):
    K_PrevIT=ModelKh(h_prevIT,ParKh)
    q_BC=-K_PrevIT[-1]*numpy.cos(incl) # cm/s !DIRECTION Z-AXIS!
    dz=z[-1]-z[-2]
    K12=(K_PrevIT[-2]+K_PrevIT[-1])/2
    theta_prevIT=ModelSWRC(h_prevIT,ParSWRC)[-1]
    theta_prevT=ModelSWRC(h_prevT,ParSWRC)[-1]
    Cap_prevIT=ModelCap(h_prevIT,ParCap)[-1]

    P[-1,-1] = dz/(2*dt_now) * Cap_prevIT + K12/dz
    P[-1,-2]=-K12/dz
    P[-2,-1]=-K12/dz
    F[-1] = dz/(2*dt_now) * Cap_prevIT * h_prevIT[-1] - dz/(2*dt_now)*(theta_prevIT-theta_prevT) - K12*numpy.cos(incl) - q_BC
    
    return P,F
    
#==============================================================================
# Atmospheric issues
#==============================================================================

def ha2Hr_Feddes(ha,T=300.0):
    ha=ha/100.0 # Convert cm to m!
    M=0.018015 # kg/mokl
    g=9.81 # m/s**2
    R=8.314 # J/(mol Kg)
    
    Hr=numpy.exp(-ha*M*g/(R*T))
    return Hr

def Hr2ha_Feddes(Hr,T=300.0):
    M=0.018015 # kg/mokl
    g=9.81 # m/s**2
    R=8.314 # J/(mol Kg)
    
    ha=-R*T/(M*g)*numpy.log(Hr)
    
    ha=ha*100.0 # Convert m to cm!
    return ha # cm

def PennmanMonteith(T,Hr,Rn,Uz,zm=1.5,zh=1.5):
#    T = Temperature (degreeC)
#    Rn = Net Radiation (MJ/m2/d)
#    Hr =Relative Humidity (-)
#    Uz = measured wind speed (m/s)
#    ed = actual vapor pressure (kPa)
#    zm = 1.5 # windspeed measurement height (m)
#    zh = 1.5 # Temperature measurement height (m)

#    lambd = Latent heat (MJ/kg)
#    delta = slope of vapor pressure curve (kPa/degreeC)
#    gamma = psychrometric constant (kPa/degreeC)
#    rho = # Atmospheric density (kg/m3)
#    ea = saturation vapor pressure (kPa)
#    ra = aerodynamic resistance (s/m)

    G = 0.0 # Soil heat flux (MJ/m2/d)
    cp = 1.013 # specific heat of moist air (kJ/kg/degreeC)
    rc = 0.0 # crop canopy resistance (s/m)
    P = 101.3 # Atmospheric Pressure (kPa)
    
    lambd = 2.501 - (2.361e-3)*T
    ea = 0.611*numpy.exp(17.27*T/(T+237.3))
    delta = 4098.0*ea/(T+237.3)**2
    gamma = 0.00163*P/lambd
    rho = 3.486 * P/(1.01*(T+273))
    ed = ea*Hr
    
    hc = 0.12 # Crop height (m)
    d = 2.0/3.0 * hc # zero plane displacement (m)
    zom = 0.123 * hc # roughness parameter
    zoh = 0.1*zom
    k = 0.41 # von karmann constant
    
    ra = numpy.log((zm-d)/zom)*numpy.log((zh-d)/zoh)/(k**2*Uz)
    
    ET0=1.0/lambd*(delta*(Rn - G)/(delta+gamma*(1.0+rc/ra))+rho*cp*(ea-ed/(ra*(delta+gamma*(1.0+rc/ra)))))/10.0 # cm/d
    
    # Condensation
    ET0[ET0<0]=0.0 # No Condensation
#    ET0[ET0<0]=ET0[ET0<0]*0.2 # Correction factor Condensation (Penman Monteith overestimates)
    
    return ET0
    
#==============================================================================
# Mass Balance
#==============================================================================

def Flux_1D(z,h,incl,ModelKh,ParKh):
    h=numpy.array(h)
    dz12=z[1:]-z[:-1]
    dh12=h[1:]-h[:-1]
    K12=(ModelKh(h[1:],map(lambda x : x[1:],ParKh))+ModelKh(h[:-1],map(lambda x : x[:-1],ParKh)))/2
    q12=(-K12*(dh12/dz12+numpy.cos(incl)))
    return q12
    
    
    
def Flux_Top(z,h,dt,incl,ModelSWRC,ParSWRC,ModelKh,ParKh):
    theta=ModelSWRC(h,ParSWRC)
    dz12=z[1]-z[0]
    dh12=h[:,1]-h[:,0]
    K12=(ModelKh(h[:,1],map(lambda x : x[None,1],ParKh))+ModelKh(h[:,0],map(lambda x : x[None,0],ParKh)))/2
    q12=(-K12*(dh12/dz12+numpy.cos(incl)))[1:]
    dV=(theta[1:,0] - theta[:-1,0])
    
    q_top=(dV*dz12)/(dt*2.0)+q12
    return q_top

def Flux_Bottom(z,h,dt,incl,ModelSWRC,ParSWRC,ModelKh,ParKh):
    theta=ModelSWRC(h,ParSWRC)
    dz12=z[-2]-z[-1]
    dh12=h[:,-2]-h[:,-1]
    K12=(ModelKh(h[:,-2],map(lambda x : x[None,-2],ParKh))+ModelKh(h[:,-1],map(lambda x : x[None,-1],ParKh)))/2
    q12=(-K12*(dh12/dz12+numpy.cos(incl)))[1:]
    dV=(theta[1:,-1] - theta[:-1,-1])
    
    q_bottom=-(dV*dz12)/(dt*2.0)+q12
    return q_bottom

def Flux_2D(z,h,dt,incl,ModelSWRC,ParSWRC,ModelKh,ParKh):
    dz=(z[1:]-z[:-1])
    h_betweenT=(h[:-1,:]+h[1:,:])/2     # average head between time steps: nT-1 x nZ
    q=numpy.zeros_like(h_betweenT)      # flux: nT-1 x nZ
    Kh=ModelKh(h_betweenT,ParKh)
    
    K_fw=(Kh[...,:-1]+Kh[...,1:])/2 # + 1/2
    h_diff=h_betweenT[...,1:]-h_betweenT[...,:-1]+dz*numpy.cos(incl)
    
    q[...,0]=Flux_Top(z,h,dt,incl,ModelSWRC,ParSWRC,ModelKh,ParKh)
    q[...,1:-1]=(-K_fw[...,1:]*(h_diff[...,1:]/dz[1:])*dz[:-1]-K_fw[...,:-1]*(h_diff[...,:-1]/dz[:-1])*dz[1:])/(dz[1:]+dz[:-1])
    q[...,-1]=Flux_Bottom(z,h,dt,incl,ModelSWRC,ParSWRC,ModelKh,ParKh)
    return q
    
def Flux(z,h,theta,dt,incl,ModelKh,ParKh):
    dz=(z[1:]-z[:-1])
    h_betweenT=(h[:-1,:]+h[1:,:])/2     # average head between time steps: nT-1 x nZ
    q=numpy.zeros_like(h_betweenT)      # flux: nT-1 x nZ
    Kh=ModelKh(h_betweenT,ParKh)
    
    K_fw=(Kh[...,:-1]+Kh[...,1:])/2 # + 1/2
    h_diff=h_betweenT[...,1:]-h_betweenT[...,:-1]+dz*numpy.cos(incl)
    
    q[...,0]=-K_fw[...,0]*(h_diff[...,0]/dz[0])-dz[0]/2*(theta[1:,0]-theta[:-1,0])/dt
    q[...,1:-1]=(-K_fw[...,1:]*(h_diff[...,1:]/dz[1:])*dz[:-1]-K_fw[...,:-1]*(h_diff[...,:-1]/dz[:-1])*dz[1:])/(dz[1:]+dz[:-1])
    q[...,-1]=-K_fw[...,-1]*(h_diff[...,-1]/dz[-1])-dz[-1]/2*(theta[1:,-1]-theta[:-1,-1])/dt
    return q
    
def Storage(theta,z):
    dz=abs((z[1:]-z[:-1])/2)
    dz_nodes=numpy.zeros_like(z)
    dz_nodes[:-1]=dz
    dz_nodes[1:]=dz_nodes[1:]+dz
#    dz_nodes=dz_nodes
            
    S=theta*dz_nodes[None,:]
    return S
    
def dStorage_dt(theta,z,dt):
    S=Storage(theta,z)
    dS=S[1:,:]-S[:-1,:]
    dSdt=dS/dt[:,None]  # dSdt : nT-1 x nZ
    return dSdt

def MassBalance(z,h,theta,dt,incl,ModelKh,ParKh):
    q=Flux(z,h,theta,dt,incl,ModelKh,ParKh)
    q_int=((q[:,0]-q[:,-1])*dt).cumsum()
    S=Storage(theta,z).sum(axis=1)
    
    Balance=S[1:]-S[0]+q_int
    return Balance
    
def MassBalance_RelErr(z,h,theta,dt,incl,ModelKh,ParKh):
    q=Flux(z,h,theta,dt,incl,ModelKh,ParKh)
    q_int=((q[:,0]-q[:,-1])*dt).cumsum()
    S=Storage(theta,z).sum(axis=1)
    
    Balance=S[1:]-S[0]+q_int
    RelErr=abs(Balance)/(numpy.maximum(abs(S[1:]-S[0]),abs(q_int)))
    return Balance,RelErr
