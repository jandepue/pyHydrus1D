# -*- coding: utf-8 -*-
"""
Created on Tue Oct 20 18:01:16 2015

@author: Administrator
"""

#==============================================================================
# All Retention/Conductivity/Capacity Functions you'll ever need
#==============================================================================
import numpy

def PositiveH_Correction(h,X,replaceX):
    if numpy.size(h)>1:
        if numpy.size(replaceX)>1:
            if numpy.ndim(X)==1:
                X[numpy.where(h>=0)]=replaceX[numpy.where(h>=0)]
            else:
                replaceX=numpy.tile(replaceX,(numpy.shape(X)[0],1))
                X[numpy.where(h>=0)]=replaceX[numpy.where(h>=0)]    
        else:
            X[numpy.where(h>=0)]=replaceX
    else:
        if h >= 0:
            X=replaceX
    return X

def GenericReplacement(BooleanReplace,X,replaceX):
    
    if numpy.size(BooleanReplace)>1:
        if numpy.size(replaceX)>1:
            if numpy.ndim(X)==1:
                X[numpy.where(BooleanReplace)]=replaceX[numpy.where(BooleanReplace)]
            else:
                replaceX=numpy.tile(replaceX,(numpy.shape(X)[0],1))
                X[numpy.where(BooleanReplace)]=replaceX[numpy.where(BooleanReplace)]    
        else:
            X[numpy.where(BooleanReplace)]=replaceX
    else:
        if BooleanReplace:
            X=replaceX
    return X

#==============================================================================
# Retention
#==============================================================================

def h2theta_VanGenuchten4(h,Par):
    thetar,thetas,alpha,n=Par
    m=1.0-1.0/n
    
    theta=thetar+(thetas-thetar)/(1+(alpha*abs(h))**n)**m    
    
    theta=PositiveH_Correction(h,theta,thetas)
    return theta

def h2theta_VanGenuchten5(h,Par):
    thetar,thetas,alpha,n,m=Par
    
    theta=thetar+(thetas-thetar)/(1+(alpha*abs(h))**n)**m
    
    theta=PositiveH_Correction(h,theta,thetas)
    return theta

def h2theta_Durner(h,Par):
    # Double van Genuchten
    thetar,thetas,w1,alpha1,n1,alpha2,n2=Par
    m1=1.0-1.0/n1
    m2=1.0-1.0/n2
    w2=1.0-w1
    
    theta=thetar+(thetas-thetar)*(w1*h2theta_VanGenuchten4(h,[0.0,1.0,alpha1,n1])+w2*h2theta_VanGenuchten4(h,[0.0,1.0,alpha2,n2]))

    theta=PositiveH_Correction(h,theta,thetas)
    return theta
    
    
def h2theta_PDI_AdsorptiveSaturation(h,Par):
    import numpy
    if len(Par)==4:
        thetar,thetas,alpha,n=Par
        pF0=6.8
    elif len(Par)==5:
        thetar,thetas,alpha,n,pF0=Par
#    pFa=numpy.log10(1/alpha)
    pFa=numpy.log10(1/alpha)
    b=0.1+(0.2/n**2)*(1.0-numpy.exp(-(thetar/(thetas-thetar))**2))
    
    Sad=1+(1.0/(pFa-pF0))*(numpy.log10(numpy.abs(h))-pFa+b*numpy.log(1.0+numpy.exp((pFa-numpy.log10(numpy.abs(h)))/b)))
    
    Sad=PositiveH_Correction(h,Sad,1.0)
    Sad=GenericReplacement(Sad<0,Sad,0.0)
    
    return Sad


def h2theta_PDIRet(h,Par):
    if len(Par)==4:
        thetar,thetas,alpha,n=Par
        pF0=6.8
    elif len(Par)==5:
        thetar,thetas,alpha,n,pF0=Par
    h0=-10**pF0
    SadPar=[thetar,thetas,alpha,n,pF0]
    Sad=h2theta_PDI_AdsorptiveSaturation(h,SadPar)
    
#    Scap=(h2theta_VanGenuchten4(h,[0.0,1.0,alpha,n])-h2theta_VanGenuchten4(h0,[0.0,1.0,alpha,n])) \
#        /(1.0-h2theta_VanGenuchten4(h0,[0.0,1.0,alpha,n])) # Scaled to be 0 at pF0
    Scap=(h2theta_VanGenuchten4(h,[0.0,1.0,alpha,n])-h2theta_VanGenuchten4(h0,[0.0,1.0,alpha,n])+1e-9) \
        /(1.0-h2theta_VanGenuchten4(h0,[0.0,1.0,alpha,n])+1e-9) # Scaled to be 0 at pF0 
    
    theta=thetar*Sad+(thetas-thetar)*Scap

    theta=PositiveH_Correction(h,theta,thetas)
    theta=GenericReplacement(theta<0,theta,0.0)
    return theta


def h2theta_BimodalPDIRet(h,Par):
    if len(Par)==7:
        thetar,thetas,w1,alpha1,n1,alpha2,n2=Par
        pF0=6.8
    elif len(Par)==8:
        thetar,thetas,w1,alpha1,n1,alpha2,n2,pF0=Par
    w2=1.0-w1
    h0=-10**pF0
    
#    SadPar=[thetar,thetas,[alpha1,alpha2][numpy.argmax([alpha1,alpha2])],[n1,n2][numpy.argmax([alpha1,alpha2])],pF0]
    if alpha1.size>1:
        Sad_alpha=numpy.zeros(alpha1.size)
        Sad_n=numpy.zeros(alpha1.size)
        for i in range(alpha1.size):
            Sad_alpha[i]=[alpha1[i],alpha2[i]][numpy.argmax([alpha1[i],alpha2[i]])]
            Sad_n[i]=[n1[i],n2[i]][numpy.argmax([alpha1[i],alpha2[i]])]
        SadPar=[thetar,thetas,Sad_alpha,Sad_n,pF0]
    else:
        SadPar=[thetar,thetas,[alpha1,alpha2][numpy.argmax([alpha1,alpha2])],[n1,n2][numpy.argmax([alpha1,alpha2])],pF0]

    
    Sad=h2theta_PDI_AdsorptiveSaturation(h,SadPar)
    
#    Scap=((w1*h2theta_VanGenuchten4(h,[0.0,1.0,alpha1,n1])+w2*h2theta_VanGenuchten4(h,[0.0,1.0,alpha2,n2]))-(w1*h2theta_VanGenuchten4(h0,[0.0,1.0,alpha1,n1])+w2*h2theta_VanGenuchten4(h0,[0.0,1.0,alpha2,n2]))) \
#        /(1-(w1*h2theta_VanGenuchten4(h0,[0.0,1.0,alpha1,n1])+w2*h2theta_VanGenuchten4(h0,[0.0,1.0,alpha2,n2]))) # Scaled to be 0 at pF0
    Scap=((w1*h2theta_VanGenuchten4(h,[0.0,1.0,alpha1,n1])+w2*h2theta_VanGenuchten4(h,[0.0,1.0,alpha2,n2]))-(w1*h2theta_VanGenuchten4(h0,[0.0,1.0,alpha1,n1])+w2*h2theta_VanGenuchten4(h0,[0.0,1.0,alpha2,n2]))+1e-9) \
        /(1-(w1*h2theta_VanGenuchten4(h0,[0.0,1.0,alpha1,n1])+w2*h2theta_VanGenuchten4(h0,[0.0,1.0,alpha2,n2]))+1e-9) # Scaled to be 0 at pF0
    
    theta=thetar*Sad+(thetas-thetar)*Scap

    theta=PositiveH_Correction(h,theta,thetas)
    theta=GenericReplacement(theta<0,theta,0.0)
    return theta


def h2theta_BimodalPDIRet_FixedAlpha(h,Par):
    if len(Par)==7:
        thetar,thetas,w1,alpha1,n1,alpha2,n2=Par
        pF0=6.8
    elif len(Par)==8:
        thetar,thetas,w1,alpha1,n1,alpha2,n2,pF0=Par
    w2=1.0-w1
    h0=-10**pF0
    
    SadPar=[thetar,thetas,alpha1,n1,pF0]
    
    Sad=h2theta_PDI_AdsorptiveSaturation(h,SadPar)
    
#    Scap=((w1*h2theta_VanGenuchten4(h,[0.0,1.0,alpha1,n1])+w2*h2theta_VanGenuchten4(h,[0.0,1.0,alpha2,n2]))-(w1*h2theta_VanGenuchten4(h0,[0.0,1.0,alpha1,n1])+w2*h2theta_VanGenuchten4(h0,[0.0,1.0,alpha2,n2]))) \
#        /(1-(w1*h2theta_VanGenuchten4(h0,[0.0,1.0,alpha1,n1])+w2*h2theta_VanGenuchten4(h0,[0.0,1.0,alpha2,n2]))) # Scaled to be 0 at pF0
    Scap=((w1*h2theta_VanGenuchten4(h,[0.0,1.0,alpha1,n1])+w2*h2theta_VanGenuchten4(h,[0.0,1.0,alpha2,n2]))-(w1*h2theta_VanGenuchten4(h0,[0.0,1.0,alpha1,n1])+w2*h2theta_VanGenuchten4(h0,[0.0,1.0,alpha2,n2]))+1e-9) \
        /(1-(w1*h2theta_VanGenuchten4(h0,[0.0,1.0,alpha1,n1])+w2*h2theta_VanGenuchten4(h0,[0.0,1.0,alpha2,n2]))+1e-9) # Scaled to be 0 at pF0
    
    theta=thetar*Sad+(thetas-thetar)*Scap

    theta=PositiveH_Correction(h,theta,thetas)
    theta=GenericReplacement(theta<0,theta,0.0)
    return theta

## INVERSE

def theta2h_VanGenuchten(theta,m,n,alpha,thetar,thetas):
    h=((((thetas-thetar)/(theta-thetar))**(1/m)-1)**(1/n))/alpha
    if numpy.size(theta)>1:
        h[numpy.where(theta>thetas)]=0
    else:
        if theta >= thetas:
            h=0
    return h

#==============================================================================
# Conductivity
#==============================================================================

def h2K_Mualem4(h,Par):
    Ks,l,alpha,n=Par
    m=1.0-1.0/n
    
    K=Ks*h2theta_VanGenuchten5(h,[0.0,1.0,alpha,n,m])**l * (alpha*(1.0-(1.0-h2theta_VanGenuchten5(h,[0.0,1.0,alpha,n,m])**(1.0/m))**m))**2/(alpha**2)
    
    K=PositiveH_Correction(h,K,Ks)
    return K
    
def h2K_Mualem5(h,Par):
    Ks,l,alpha,n,m=Par
    
    K=Ks*h2theta_VanGenuchten5(h,[0.0,1.0,alpha,n,m])**l * (alpha*(1.0-(1.0-h2theta_VanGenuchten5(h,[0.0,1.0,alpha,n,m])**(1.0/m))**m))**2/(alpha**2)

    K=PositiveH_Correction(h,K,Ks)
    return K

def h2K_MualemBimodal(h,Par):
    Ks,l,w1,alpha1,n1,alpha2,n2=Par
    m1=1.0-1.0/n1
    m2=1.0-1.0/n2
    w2=1.0-w1
    K=Ks*(h2theta_Durner(h,[0.0,1.0,w1,alpha1,n1,alpha2,n2])**l \
        *(w1*alpha1 * (1.0 - (alpha1*numpy.abs(h))**(n1-1) * h2theta_VanGenuchten5(h,[0.0,1.0,alpha1,n1,m1])) \
        + w2*alpha2 * (1.0 - (alpha2*numpy.abs(h))**(n2-1) * h2theta_VanGenuchten5(h,[0.0,1.0,alpha2,n2,m2])))**2 \
        /(w1*alpha1+w2*alpha2)**2)

    K=PositiveH_Correction(h,K,Ks)
    return K
    
def h2K_PDI_Kfilm(h,Par):
    if len(Par)==7:
        thetar,thetas,Ks,l,omega1,alpha,n=Par
        pF0=6.8
        a=-1.5
    elif len(Par)==9:
        thetar,thetas,Ks,l,omega1,alpha,n,pF0,a=Par
    
    ha=-1.0/numpy.max(alpha)
    h0=-10**pF0
    
    SadPar=[0.0,1.0,alpha,n,pF0]
    Sad=h2theta_PDI_AdsorptiveSaturation(h,SadPar)
    Kfilm=(h0/ha)**(a*(1.0-Sad))
    
    return Kfilm
    
def h2K_PDIK(h,Par,VaporConductivity=False):
    if len(Par)==7:
        thetar,thetas,Ks,l,omega1,alpha,n=Par
        pF0=6.8
        a=-1.5
    elif len(Par)==9:
        thetar,thetas,Ks,l,omega1,alpha,n,pF0,a=Par
    
    m=1.0-1.0/n
    omega2=1.0-omega1
    ha=-1.0/alpha
    pFa=numpy.log10(ha)
    h0=-10**pF0
    Par_Ret=[thetar,thetas,alpha,n,pF0]
    
    Kfilm=h2K_PDI_Kfilm(h,Par)
    
#    Scap=(h2theta_VanGenuchten4(h,[0.0,1.0,alpha,n])-h2theta_VanGenuchten4(h0,[0.0,1.0,alpha,n])) \
#        /(1-h2theta_VanGenuchten4(h0,[0.0,1.0,alpha,n])) # Scaled to be 0 at pF0
    Scap=(h2theta_VanGenuchten4(h,[0.0,1.0,alpha,n])-h2theta_VanGenuchten4(h0,[0.0,1.0,alpha,n])+1e-9) \
        /(1-h2theta_VanGenuchten4(h0,[0.0,1.0,alpha,n])+1e-9) # Scaled to be 0 at pF0
    
    Scap=GenericReplacement(Scap<0,Scap,0.0)
    
    Kcap=Scap**l \
        *(1.0 -
        (1.0 - h2theta_VanGenuchten5(h,[0.0,1.0,alpha,n,m])**(1/m)+1e-9)**m \
        /(1.0 - h2theta_VanGenuchten5(h0,[0.0,1.0,alpha,n,m])**(1/m)+1e-9)**m \
        )**2
    
    Kcap=PositiveH_Correction(h,Kcap,1.0)
    Kcap=GenericReplacement(Kcap<0,Kcap,0.0)
    
    if VaporConductivity:
        Kvap=h2K_PDI_Kvap(h,Par_Ret,SWR_Model=h2theta_PDIRet)
        K=Ks*(omega2*Kcap+omega1*Kfilm)+Kvap
    else:
        K=Ks*(omega2*Kcap+omega1*Kfilm)
    
    return K

def h2K_BimodalPDI_Kfilm(h,Par,FixedAlpha=False):
    if len(Par)==10:
        thetar,thetas,Ks,l,w1,omega1,alpha1,n1,alpha2,n2=Par
        pF0=6.8
        a=-1.5
    elif len(Par)==12:
        thetar,thetas,Ks,l,w1,omega1,alpha1,n1,alpha2,n2,pF0,a=Par
    ha=-1.0/numpy.max([alpha1,alpha2])
    h0=-10**pF0
    
    if FixedAlpha:
        SadPar=[0.0,1.0,alpha1,n1,pF0]
    else:
        if alpha1.size>1:
            Sad_alpha=numpy.zeros(alpha1.size)
            Sad_n=numpy.zeros(alpha1.size)
            for i in range(alpha1.size):
                Sad_alpha[i]=[alpha1[i],alpha2[i]][numpy.argmax([alpha1[i],alpha2[i]])]
                Sad_n[i]=[n1[i],n2[i]][numpy.argmax([alpha1[i],alpha2[i]])]
            SadPar=[0.0,1.0,Sad_alpha,Sad_n,pF0]
        else:
            SadPar=[0.0,1.0,[alpha1,alpha2][numpy.argmax([alpha1,alpha2])],[n1,n2][numpy.argmax([alpha1,alpha2])],pF0]
    Sad=h2theta_PDI_AdsorptiveSaturation(h,SadPar)
    Kfilm=(h0/ha)**(a*(1.0-Sad))
    
    return Kfilm

def h2K_PDI_Kvap(h_in,Par,SWR_Model=h2theta_BimodalPDIRet):
    import numpy
    import copy
    h=copy.deepcopy(h_in)
    thetas=Par[1]
    rho_w=1e3 # kg/m3
    M=0.018015 # kg/mol
    g=9.81 # m/s2
    R=8.134 #J/mol/kg
    T=300 # K

    h=PositiveH_Correction(h,h,0.0)
    
    theta_A=thetas-SWR_Model(h,Par)
    Tort=theta_A**(7.0/3)/thetas**2
    D_A=2.14e-5*(T/273.15)**2
    D=theta_A*D_A*Tort
    
    rho_sv=1e-3*numpy.exp(31.3716-6014.76/T-7.924951e-3*T)/T
    Hr=numpy.exp(-numpy.abs(h*100.0)*M*g/(R*T)) # Convert h to meters!
    
    Kvap=rho_sv/rho_w*Hr*D*M*g/(R*T)
    
    return Kvap

def h2K_BimodalPDIK(h,Par,VaporConductivity=False):
    if len(Par)==10:
        thetar,thetas,Ks,l,w1,omega1,alpha1,n1,alpha2,n2=Par
        pF0=6.8
        a=-1.5
    elif len(Par)==12:
        thetar,thetas,Ks,l,w1,omega1,alpha1,n1,alpha2,n2,pF0,a=Par
    m1=1.0-1.0/n1
    m2=1.0-1.0/n2
    w2=1.0-w1
    omega2=1.0-omega1
    
    h0=-10**pF0
    Par_Ret=[thetar,thetas,w1,alpha1,n1,alpha2,n2,pF0]
    
    Kfilm=h2K_BimodalPDI_Kfilm(h,Par)
    Scap=((w1*h2theta_VanGenuchten4(h,[0.0,1.0,alpha1,n1])+w2*h2theta_VanGenuchten4(h,[0.0,1.0,alpha2,n2]))-(w1*h2theta_VanGenuchten4(h0,[0.0,1.0,alpha1,n1])+w2*h2theta_VanGenuchten4(h0,[0.0,1.0,alpha2,n2]))+1e-9) \
        /(1-(w1*h2theta_VanGenuchten4(h0,[0.0,1.0,alpha1,n1])+w2*h2theta_VanGenuchten4(h0,[0.0,1.0,alpha2,n2]))+1e-9) # Scaled to be 0 at pF0
        
    Scap=GenericReplacement(Scap<0,Scap,0.0)
    
    Kcap=Scap**l \
        *(1.0 -
        (w1*alpha1 * (1.0 - h2theta_VanGenuchten5(h,[0.0,1.0,alpha1,n1,m1])**(1/m1))**m1 \
        + w2*alpha2 * (1.0 - h2theta_VanGenuchten5(h,[0.0,1.0,alpha2,n2,m2])**(1/m2))**m2+1e-9) \
        /(w1*alpha1 * (1.0 - h2theta_VanGenuchten5(h0,[0.0,1.0,alpha1,n1,m1])**(1/m1))**m1 \
        + w2*alpha2 * (1.0 - h2theta_VanGenuchten5(h0,[0.0,1.0,alpha2,n2,m2])**(1/m2))**m2+1e-9) )**2
    
    Kcap=PositiveH_Correction(h,Kcap,1.0)
    Kcap=GenericReplacement(Kcap<0,Kcap,0.0)
    
    if VaporConductivity:
        Kvap=h2K_PDI_Kvap(h,Par_Ret,SWR_Model=h2theta_BimodalPDIRet)
        K=Ks*(omega2*Kcap+omega1*Kfilm)+Kvap
    else:
        K=Ks*(omega2*Kcap+omega1*Kfilm)
    
    return K

def h2K_BimodalPDIK_FixedAlpha(h,Par,VaporConductivity=False):
    if len(Par)==10:
        thetar,thetas,Ks,l,w1,omega1,alpha1,n1,alpha2,n2=Par
        pF0=6.8
        a=-1.5
    elif len(Par)==12:
        thetar,thetas,Ks,l,w1,omega1,alpha1,n1,alpha2,n2,pF0,a=Par
    m1=1.0-1.0/n1
    m2=1.0-1.0/n2
    w2=1.0-w1
    omega2=1.0-omega1
    
    h0=-10**pF0
    Par_Ret=[thetar,thetas,w1,alpha1,n1,alpha2,n2,pF0]
    
    Kfilm=h2K_BimodalPDI_Kfilm(h,Par,FixedAlpha=True)
    Scap=((w1*h2theta_VanGenuchten4(h,[0.0,1.0,alpha1,n1])+w2*h2theta_VanGenuchten4(h,[0.0,1.0,alpha2,n2]))-(w1*h2theta_VanGenuchten4(h0,[0.0,1.0,alpha1,n1])+w2*h2theta_VanGenuchten4(h0,[0.0,1.0,alpha2,n2]))+1e-9) \
        /(1-(w1*h2theta_VanGenuchten4(h0,[0.0,1.0,alpha1,n1])+w2*h2theta_VanGenuchten4(h0,[0.0,1.0,alpha2,n2]))+1e-9) # Scaled to be 0 at pF0
        
    Scap=GenericReplacement(Scap<0,Scap,0.0)
    
    Kcap=Scap**l \
        *(1.0 -
        (w1*alpha1 * (1.0 - h2theta_VanGenuchten5(h,[0.0,1.0,alpha1,n1,m1])**(1/m1))**m1 \
        + w2*alpha2 * (1.0 - h2theta_VanGenuchten5(h,[0.0,1.0,alpha2,n2,m2])**(1/m2))**m2+1e-9) \
        /(w1*alpha1 * (1.0 - h2theta_VanGenuchten5(h0,[0.0,1.0,alpha1,n1,m1])**(1/m1))**m1 \
        + w2*alpha2 * (1.0 - h2theta_VanGenuchten5(h0,[0.0,1.0,alpha2,n2,m2])**(1/m2))**m2+1e-9) )**2
    
    Kcap=PositiveH_Correction(h,Kcap,1.0)
    Kcap=GenericReplacement(Kcap<0,Kcap,0.0)
    
    if VaporConductivity:
        Kvap=h2K_PDI_Kvap(h,Par_Ret,SWR_Model=h2theta_BimodalPDIRet)
        K=Ks*(omega2*Kcap+omega1*Kfilm)+Kvap
    else:
        K=Ks*(omega2*Kcap+omega1*Kfilm)
    
    return K


#==============================================================================
# Capacity
#==============================================================================

def h2C_VanGenuchten4(h,Par):
    thetar,thetas,alpha,n=Par
    m=1.0-1.0/n
    
    C=-(thetas-thetar)*(n*alpha**n*abs(h)**(n-1))*(-m)*(1+alpha**n*abs(h)**n)**(-m-1)

    C=PositiveH_Correction(h,C,0.0)
    
    return C
    
def h2C_VanGenuchten5(h,Par):
    thetar,thetas,alpha,n,m=Par
    
    C=-(thetas-thetar)*(n*alpha**n*abs(h)**(n-1))*(-m)*(1+alpha**n*abs(h)**n)**(-m-1)

    C=PositiveH_Correction(h,C,0.0)
    return C


def h2C_PDICap(h,Par):
    if len(Par)==4:
        thetar,thetas,alpha,n=Par
        pF0=6.8
    elif len(Par)==5:
        thetar,thetas,alpha,n,pF0=Par
    h0=-10**pF0
    pFa=numpy.log10(1/alpha)
    pF=numpy.log10(numpy.abs(h))
    m=1.0-1.0/n
    
    b=0.1+(0.2/n**2)*(1.0-numpy.exp(-(thetar/(thetas-thetar))**2))
    
    Cap_AD=(abs(h)*numpy.log(10.0)*(pFa-pF0))**(-1.0) * (1.0 - numpy.exp((pFa-pF)/b)/(1.0+numpy.exp((pFa-pF)/b)))
    
    dSdh=-alpha*n*m*(alpha*abs(h))**(n-1.0)*(1+(alpha*abs(h))**n)**(-m-1.0)
    
    C=(thetas-thetar+1e-9)/(1.0-h2theta_VanGenuchten5(h0,[0.0,1.0,alpha,n,m])+1e-9) * dSdh + thetar*Cap_AD
    C=-C
    
    C=PositiveH_Correction(h,C,0.0)

    return C

def h2C_BimodalPDICap(h,Par):
    if len(Par)==7:
        thetar,thetas,w1,alpha1,n1,alpha2,n2=Par
        pF0=6.8
    elif len(Par)==8:
        thetar,thetas,w1,alpha1,n1,alpha2,n2,pF0=Par
    w2=1.0-w1
    m1=1.0-1.0/n1
    m2=1.0-1.0/n2
    h0=-10**pF0
    
    if alpha1.size>1:
        Sad_alpha=numpy.zeros(alpha1.size)
        Sad_n=numpy.zeros(alpha1.size)
        for i in range(alpha1.size):
            Sad_alpha[i]=[alpha1[i],alpha2[i]][numpy.argmax([alpha1[i],alpha2[i]])]
            Sad_n[i]=[n1[i],n2[i]][numpy.argmax([alpha1[i],alpha2[i]])]
        SadPar=[thetar,thetas,Sad_alpha,Sad_n,pF0]
    else:
        SadPar=[thetar,thetas,[alpha1,alpha2][numpy.argmax([alpha1,alpha2])],[n1,n2][numpy.argmax([alpha1,alpha2])],pF0]
        
    pFa=numpy.log10(1.0/SadPar[2])
    pF=numpy.log10(numpy.abs(h))
    b=0.1+(0.2/SadPar[3]**2)*(1.0-numpy.exp(-(thetar/(thetas-thetar))**2))
    
    Cap_AD=(abs(h)*numpy.log(10.0)*(pFa-pF0))**(-1.0) * (1.0 - numpy.exp((pFa-pF)/b)/(1.0+numpy.exp((pFa-pF)/b)))
    
    dSdh= w1 * (-alpha1*n1*m1*(alpha1*abs(h))**(n1-1.0)*(1+(alpha1*abs(h))**n1)**(-m1-1.0)) \
        + w2 * (-alpha2*n2*m2*(alpha2*abs(h))**(n2-1.0)*(1+(alpha2*abs(h))**n2)**(-m2-1.0))
    
    C=(thetas-thetar+1e-9)/(1.0-(w1*h2theta_VanGenuchten4(h0,[0.0,1.0,alpha1,n1])+w2*h2theta_VanGenuchten4(h0,[0.0,1.0,alpha2,n2]))+1e-9) * dSdh + thetar*Cap_AD
    C=-C
    
    C=PositiveH_Correction(h,C,0.0)

    return C

def h2C_BimodalPDICap_FixedAlpha(h,Par):
    if len(Par)==7:
        thetar,thetas,w1,alpha1,n1,alpha2,n2=Par
        pF0=6.8
    elif len(Par)==8:
        thetar,thetas,w1,alpha1,n1,alpha2,n2,pF0=Par
    w2=1.0-w1
    m1=1.0-1.0/n1
    m2=1.0-1.0/n2
    h0=-10**pF0
    
    SadPar=[thetar,thetas,alpha1,n1,pF0]
        
    pFa=numpy.log10(1.0/SadPar[2])
    pF=numpy.log10(numpy.abs(h))
    b=0.1+(0.2/SadPar[3]**2)*(1.0-numpy.exp(-(thetar/(thetas-thetar))**2))
    
    Cap_AD=(abs(h)*numpy.log(10.0)*(pFa-pF0))**(-1.0) * (1.0 - numpy.exp((pFa-pF)/b)/(1.0+numpy.exp((pFa-pF)/b)))
    
    dSdh= w1 * (-alpha1*n1*m1*(alpha1*abs(h))**(n1-1.0)*(1+(alpha1*abs(h))**n1)**(-m1-1.0)) \
        + w2 * (-alpha2*n2*m2*(alpha2*abs(h))**(n2-1.0)*(1+(alpha2*abs(h))**n2)**(-m2-1.0))
    
    C=(thetas-thetar+1e-9)/(1.0-(w1*h2theta_VanGenuchten4(h0,[0.0,1.0,alpha1,n1])+w2*h2theta_VanGenuchten4(h0,[0.0,1.0,alpha2,n2]))+1e-9) * dSdh + thetar*Cap_AD
    C=-C
    
    C=PositiveH_Correction(h,C,0.0)

    return C

