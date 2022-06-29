# -*- coding: utf-8 -*-
"""
Created on Mon Oct 26 13:37:36 2015

@author: Administrator
"""

def RunAlbonHydrus(FC,SFC,HC,K):
#    FC='Z'
#    SFC='REF'
#    HC='C'
#    K='KSK0'
#    K='KS'
#    K='K0'
    
    import os
    import pylab
    import numpy
    import copy 
    
    #from Hydrus_Funky import *
    from DIY_Hydrus_Funky import DoHydrus,gausselim,F_j1k,P_j1k,BC_DiricheletTop,BC_DiricheletBottom,BC_NeumannTop,BC_NeumannBottom,BC_FreeDrainageBottom,ha2Hr_Feddes,Hr2ha_Feddes,PennmanMonteith,Flux_1D,Flux_Top,Flux_Bottom,Flux_2D,Flux,Storage,dStorage_dt,MassBalance,MassBalance_RelErr
    from RetentionConductivityCapacity_Funky import PositiveH_Correction,GenericReplacement,h2theta_VanGenuchten4,h2theta_VanGenuchten5,h2theta_Durner,h2theta_PDI_AdsorptiveSaturation,h2theta_PDIRet,h2theta_BimodalPDIRet,h2theta_BimodalPDIRet_FixedAlpha,theta2h_VanGenuchten,h2K_Mualem4,h2K_Mualem5,h2K_MualemBimodal,h2K_PDI_Kfilm,h2K_PDIK,h2K_BimodalPDI_Kfilm,h2K_PDI_Kvap,h2K_BimodalPDIK,h2K_BimodalPDIK_FixedAlpha,h2C_VanGenuchten4,h2C_VanGenuchten5,h2C_PDICap,h2C_BimodalPDICap,h2C_BimodalPDICap_FixedAlpha
    
    
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
    # Open Control
    #==============================================================================
    
    #ProjectFolder='C:\\Users\\Administrator\\Documents\\Arbeid\\Nevenactiviteiten\\Hydrus_Matlab\\test_MakeItRain'
    #
    #TLevelOutName=os.path.join(ProjectFolder,'T_Level.out')
    #TlevelData,TLevel_Header=ReadTLevelOut(TLevelOutName)
    #
    #NodInfOutName=os.path.join(ProjectFolder,'Nod_Inf.out')
    #NodInf_T,NodInf_Data,NodInf_Header=ReadNodInfOut(NodInfOutName)
    #
    #ObsOutName=os.path.join(ProjectFolder,'Obs_Node.out')
    #ObsNode_Nodes,ObsNode_t,ObsNode_L=ReadObsNodeOut(ObsOutName)
    
    #t_hydrus=TlevelData[:,0]
    #dT_Hydrus=numpy.append(t_hydrus[0],t_hydrus[1:]-t_hydrus[:-1])
    
    #==============================================================================
    # Open Boundary Conditions
    #==============================================================================
    
    #BCfilename='C:\\Users\\Administrator\\Documents\\GitHub\\Hydrus_DIY\\DataLanduyt\\Landuyt_BC.txt'
    BCfilename='/home/supersoil/Documents/Jan/git/Hydrus_DIY/DataLanduyt/Landuyt_BC.txt'
    
    fID=open(BCfilename)
    BC_Header=numpy.array(fID.readline().replace('\n','').split('\t'))
    fID.close()
    
    BC_All=numpy.genfromtxt(BCfilename,delimiter='\t',skip_header=1,dtype='str')
    
    BC_T=BC_All[:,1].astype('float')*3600*24
    BC_GWT=BC_All[:,2].astype('float')
    BC_P=BC_All[:,3].astype('float')*1e-1/3600 # cm / s
    BC_Rn=BC_All[:,4].astype('float')*1e-6*3600*24 # MJ/m2/d
    BC_Temp=BC_All[:,5].astype('float') # DegreeC
    BC_Hr=BC_All[:,6].astype('float')*1e-2
    BC_Uz=BC_All[:,7].astype('float')
    
    BC_Ha=Hr2ha_Feddes(BC_Hr,T=BC_Temp+273.0)
    BC_ET0=PennmanMonteith(BC_Temp,BC_Hr,BC_Rn,BC_Uz,zm=1.75,zh=1.75) / 24.0 / 3600 # cm/s
    
    fig=pylab.figure()
    ax=fig.add_subplot(611)
    ax.plot(BC_T,BC_GWT)
    ax.set_ylabel('GWT (cm)')
    ax=fig.add_subplot(612)
    ax.plot(BC_T,BC_P)
    ax.set_ylabel('Precipitation (cm/s)')
    ax=fig.add_subplot(613)
    ax.plot(BC_T,BC_Rn)
    ax.set_ylabel('Net Radiation (MJ/m2/d)')
    ax=fig.add_subplot(614)
    ax.plot(BC_T,BC_Temp)
    ax.set_ylabel('Temperature (K)')
    ax=fig.add_subplot(615)
    ax.plot(BC_T,BC_Hr)
    ax.set_ylabel('Relative Humidity (%)')
    ax=fig.add_subplot(616)
    ax.plot(BC_T,BC_Uz)
    #ax.set_yscale('log')
    ax.set_ylabel('Wind Speed (m/s)')
    ax.set_xlabel('Time (h)')
    figlist.append(fig)
    
    fig=pylab.figure()
    ax=fig.add_subplot(411)
    ax.plot(BC_T,BC_GWT)
    ax.set_ylabel('GWT (cm)')
    ax=fig.add_subplot(412)
    ax.plot(BC_T,BC_P)
    ax.set_ylabel('Precipitation (cm/s)')
    ax=fig.add_subplot(413)
    ax.plot(BC_T,BC_ET0)
    ax.set_ylabel('ET0 (cm/s)')
    ax=fig.add_subplot(414)
    ax.plot(BC_T,-BC_Ha)
    #ax.set_yscale('log')
    ax.set_ylabel('h air (Feddes) (cm)')
    figlist.append(fig)
    
    BC_q=-BC_P+BC_ET0
    
    #==============================================================================
    # Initial
    #==============================================================================
    
    nz=261
    #nz=131
    zmin=0
    zmax=-130
    z=numpy.linspace(zmin,zmax,nz)    # cm // choose z-axis upwards or downwards!
    #dz=(z[2:]-z[:-2])/2
    
    t0=0.
    #t_end=1*24*60*60.                 # s
    #t_end=2e4                # s
    t_end=BC_T[-1]
    
    dt_init=1e-1                 # s
    dt_BC=10.0              # s
    dt_min=1e-4                # s
    dt_max=1e3                 # s
    
    h0=numpy.ones(nz)*-100      # cm
    #h0[0]=5                     # cm
    #h0[0]=0                     # cm
    #h0[-1]=-100                  # cm
    h0=BC_GWT[0]-z
    
    
    #t_rain_L=[5e4,2*24*60*60.]
    #q_rain_L=[-4e-4,3e-4] # cm/s
    #t_rain_L=[5e4,8e4,2*24*60*60.]
    #q_rain_L=[-4e-4,3e-4,-8e-4] # cm/s
    t_rain_L=list(BC_T) # s
    q_rain_L=list(BC_q) # cm/s
    
    
    #q_bottom=0  # cm/s
    
    #==============================================================================
    # Iteration parameters
    #==============================================================================
    
    theta_tol=0.001     # m3/m3
    h_tol=0.1             # cm
    q_tol=1e-6
    pF_tol=0.1
    
    #theta_tol=0.001     # m3/m3
    #h_tol=0.1             # cm
    #q_tol=1e5
    #pF_tol=1e5
    
    #theta_tol_rel=1e-2     # m3/m3
    #h_tol_rel=1e-2             # cm
    #q_tol_rel=1e-2
    
    nit_decr_dt=6
    decr_dt=0.8
    nit_incr_dt=2
    incr_dt=1.2
    
    nit_term=20
    dt_term=0.4
    
    maxHesitation_N=500
    maxHesitation_T=2e1
    
    #==============================================================================
    # Soil properties
    #==============================================================================
    
    #incl=numpy.pi/2      # inclinatie v z-as tov verticale
    incl=0
    
    ## Homogeneous Loam
    
    #n=numpy.ones(nz)*1.56
    ##m=numpy.ones(nz)*0.6
    #m=1-1/n
    #alpha=numpy.ones(nz)*0.036  # 1/cm
    #thetar=numpy.ones(nz)*0.078
    #thetas=numpy.ones(nz)*0.43
    #l=numpy.ones(nz)*0.5
    #Ks=numpy.ones(nz)*2.88889e-4   # cm/s
    #
    #ModelSWRC=h2theta_VanGenuchten5
    #ParSWRC=[thetar,thetas,alpha,n,m]
    #ModelKh=h2K_Mualem5
    #ParKh=[Ks,l,alpha,n,m]
    #ModelCap=h2C_VanGenuchten5
    #ParCap=[thetar,thetas,alpha,n,m]
    
    ## Albon 3 Layer
    
    #Albonfilename='Z:\\shares\\BW12\\OG Bodemfysica\\ALBON soil compaction\\albon data\\Albon for future research\\Albon_AllHydraulic.txt'
    Albonfilename='/home/supersoil/Documents/Jan/git/Hydrus_DIY/DataLanduyt/Albon_AllHydraulic.txt'
    fID=open(Albonfilename)
    Albon_Header=numpy.array(fID.readline().replace('\n','').split('\t'))
    fID.close()
    
    Albon_All=numpy.genfromtxt(Albonfilename,delimiter='\t',skip_header=1,dtype='str')
    
#    FC='Z'
#    SFC='REF'
#    HC='C'
#    K='KSK0'
#    K='KS'
#    K='K0' 
    
    Layers=['TOP','CSUB','SUB']
    
    
    ## KSK0 or K0
    
    if (K=='K0') or (K=='KSK0'):
    
        thetar_WR=numpy.ones(nz)
        thetas_WR=numpy.ones(nz)
        w1_WR=numpy.ones(nz)
        alpha1_WR=numpy.ones(nz)
        n1_WR=numpy.ones(nz)
        alpha2_WR=numpy.ones(nz)
        n2_WR=numpy.ones(nz)
        thetar_K=numpy.ones(nz)
        thetas_K=numpy.ones(nz)
        Ks_K=numpy.ones(nz)
        l_K=numpy.ones(nz)
        w1_K=numpy.ones(nz)
        omega1_K=numpy.ones(nz)
        alpha1_K=numpy.ones(nz)
        n1_K=numpy.ones(nz)
        alpha2_K=numpy.ones(nz)
        n2_K=numpy.ones(nz)
        SoilProp=[thetar_WR,thetas_WR,w1_WR,alpha1_WR,n1_WR,alpha2_WR,n2_WR,
                  thetar_K,thetas_K,Ks_K,l_K,w1_K,omega1_K,alpha1_K,n1_K,alpha2_K,n2_K]
        
        if K == 'K0':
            SoilProp_Keys=['K0_PDI_Bimodal_SWR thetar (m3/m3)','K0_PDI_Bimodal_SWR thetas (m3/m3)', 'K0_PDI_Bimodal_SWR w1 (-)', 'K0_PDI_Bimodal_SWR alpha1 (1/cm)', 'K0_PDI_Bimodal_SWR n1 (-)','K0_PDI_Bimodal_SWR alpha2 (1/cm)', 'K0_PDI_Bimodal_SWR n2 (-)',
                   'K0_PDI_Bimodal_K thetar (m3/m3)', 'K0_PDI_Bimodal_K thetas (m3/m3)', 'K0_PDI_Bimodal_K Ks (cm/d)', 'K0_PDI_Bimodal_K l (-)', 'K0_PDI_Bimodal_K w1 (-)', 'K0_PDI_Bimodal_K omega1 (-)', 'K0_PDI_Bimodal_K alpha1 (1/cm)','K0_PDI_Bimodal_K n1 (-)', 'K0_PDI_Bimodal_K alpha2 (1/cm)', 'K0_PDI_Bimodal_K n2 (-)',]
        elif K == 'KSK0':
            SoilProp_Keys=['KS_PDI_Bimodal_SWR thetar (m3/m3)','KS_PDI_Bimodal_SWR thetas (m3/m3)', 'KS_PDI_Bimodal_SWR w1 (-)', 'KS_PDI_Bimodal_SWR alpha1 (1/cm)', 'KS_PDI_Bimodal_SWR n1 (-)','KS_PDI_Bimodal_SWR alpha2 (1/cm)', 'KS_PDI_Bimodal_SWR n2 (-)',
               'KS_PDI_Bimodal_K thetar (m3/m3)', 'KS_PDI_Bimodal_K thetas (m3/m3)', 'KS_PDI_Bimodal_K Ks (cm/d)', 'KS_PDI_Bimodal_K l (-)', 'KS_PDI_Bimodal_K w1 (-)', 'KS_PDI_Bimodal_K omega1 (-)', 'KS_PDI_Bimodal_K alpha1 (1/cm)','KS_PDI_Bimodal_K n1 (-)', 'KS_PDI_Bimodal_K alpha2 (1/cm)', 'KS_PDI_Bimodal_K n2 (-)',]
        
        nSP=len(SoilProp)
        
        Z1_L=[]
        for L in Layers:
            filt=(Albon_All[:,0]==FC) & (Albon_All[:,1]==SFC) & (Albon_All[:,2]==HC) & (Albon_All[:,3]==L)
            Data_Layer=Albon_All[filt,:][0]
            Z1=float(Data_Layer[Albon_Header=='Z1 (cm)'][0])
            Z2=float(Data_Layer[Albon_Header=='Z2 (cm)'][0])
            if L == 'SUB':
                Z2=-zmax
            print('Set layer %s : %s cm - %s cm'%(L,Z1,Z2))
            for iL in range(nSP):
                SoilProp[iL][(z<=-Z1) & (z>=-Z2)] = float(Data_Layer[Albon_Header==SoilProp_Keys[iL]][0])
            Z1_L.append(Z1)
        
        [thetar_WR,thetas_WR,w1_WR,alpha1_WR,n1_WR,alpha2_WR,n2_WR,thetar_K,thetas_K,Ks_K,l_K,w1_K,omega1_K,alpha1_K,n1_K,alpha2_K,n2_K] = SoilProp
        Ks_K=Ks_K/(24*3600) # cm/d to cm/s
        
        # Fix alpha
        filter_fixAlpha=alpha1_WR < alpha2_WR
        alpha2_WR_T=copy.deepcopy(alpha2_WR)
        n2_WR_T=copy.deepcopy(n2_WR)
        alpha2_WR[filter_fixAlpha]=alpha1_WR[filter_fixAlpha]
        alpha1_WR[filter_fixAlpha]=alpha2_WR_T[filter_fixAlpha]
        n2_WR[filter_fixAlpha]=n1_WR[filter_fixAlpha]
        n1_WR[filter_fixAlpha]=n2_WR_T[filter_fixAlpha]
        w1_WR[filter_fixAlpha]=1-w1_WR[filter_fixAlpha]
        
        filter_fixAlpha=alpha1_K < alpha2_K
        alpha2_K_T=copy.deepcopy(alpha2_K)
        n2_K_T=copy.deepcopy(n2_K)
        alpha2_K[filter_fixAlpha]=alpha1_K[filter_fixAlpha]
        alpha1_K[filter_fixAlpha]=alpha2_K_T[filter_fixAlpha]
        n2_K[filter_fixAlpha]=n1_K[filter_fixAlpha]
        n1_K[filter_fixAlpha]=n2_K_T[filter_fixAlpha]
        w1_K[filter_fixAlpha]=1-w1_K[filter_fixAlpha]
        
        ModelSWRC=h2theta_BimodalPDIRet_FixedAlpha
        ParSWRC=[thetar_WR,thetas_WR,w1_WR,alpha1_WR,n1_WR,alpha2_WR,n2_WR]
        ModelKh=h2K_BimodalPDIK_FixedAlpha
        ParKh=[thetar_K,thetas_K,Ks_K,l_K,w1_K,omega1_K,alpha1_K,n1_K,alpha2_K,n2_K]
        ModelCap=h2C_BimodalPDICap_FixedAlpha
        ParCap=[thetar_WR,thetas_WR,w1_WR,alpha1_WR,n1_WR,alpha2_WR,n2_WR]
    
    
    ## KS
    
    elif K=='KS':
        thetar=numpy.ones(nz)
        thetas=numpy.ones(nz)
        alpha=numpy.ones(nz)
        n=numpy.ones(nz)
        m=numpy.ones(nz)
        Ks=numpy.ones(nz)
        SoilProp=[thetar,thetas,alpha,n,m,Ks]
        
        SoilProp_Keys=['K0_VanGenuchten thetar (m3/m3)', 'K0_VanGenuchten thetas (m3/m3)', 'K0_VanGenuchten alpha (1/cm)', 'K0_VanGenuchten n (-)', 'K0_VanGenuchten m (-)','Ks - geometrisch (m/s) ']
        nSP=len(SoilProp)
        
        Z1_L=[]
        for L in Layers:
            filt=(Albon_All[:,0]==FC) & (Albon_All[:,1]==SFC) & (Albon_All[:,2]==HC) & (Albon_All[:,3]==L)
            Data_Layer=Albon_All[filt,:][0]
            Z1=float(Data_Layer[Albon_Header=='Z1 (cm)'][0])
            Z2=float(Data_Layer[Albon_Header=='Z2 (cm)'][0])
            if L == 'SUB':
                Z2=-zmax
            print('Set layer %s : %s cm - %s cm'%(L,Z1,Z2))
            for iL in range(nSP):
                SoilProp[iL][(z<=-Z1) & (z>=-Z2)] = float(Data_Layer[Albon_Header==SoilProp_Keys[iL]][0])
            Z1_L.append(Z1)
        
        [thetar,thetas,alpha,n,m,Ks] = SoilProp
#        Ks=Ks/(24*3600) # cm/d to cm/s
        Ks=Ks*100.0 # m/s to cm/s
        l=numpy.ones(nz)*0.5
        
        ModelSWRC=h2theta_VanGenuchten5
        ParSWRC=[thetar,thetas,alpha,n,m]
        ModelKh=h2K_Mualem5
        ParKh=[Ks,l,alpha,n,m]
        ModelCap=h2C_VanGenuchten5
        ParCap=[thetar,thetas,alpha,n,m]
    
    else:
        raise NameError('K_NotCorrectlyDefined')
        
    ## PLOT
    
    Tensionplot=-numpy.logspace(-1,7,1000)
    fig=pylab.figure()
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
    #ax21.set_xlim([1e-2,1e5])
    #ax21.set_ylim([0.0,0.5])
    #ax22.set_xlim([1e-2,1e5])
    #ax22.set_ylim([1e-6,1e4])
    #ax23.set_xlim([1e-2,1e5])
    #ax23.set_ylim([1e-6,1e4])
    ax21.set_xscale('log')
    ax22.set_xscale('log')
    ax23.set_xscale('log')
    ax22.set_yscale('log')
    nL=len(Layers)
    for iL in range(nL):
        color=cm(float(iL)/(nL-0.999))
        ax21.plot(-Tensionplot ,ModelSWRC(Tensionplot,map(lambda x : x[z==-Z1_L[iL]],ParSWRC)),'-',color=color,label='%s'%Layers[iL])
        ax22.plot(-Tensionplot ,ModelKh(Tensionplot,map(lambda x : x[z==-Z1_L[iL]],ParKh)),'-',color=color,label='%s: Ks = %1.1e cm/s'%(Layers[iL],ParKh[2][z==-Z1_L[iL]][0]))
        ax23.plot(-Tensionplot ,ModelCap(Tensionplot,map(lambda x : x[z==-Z1_L[iL]],ParCap)),'-',color=color,label='%s'%Layers[iL])
    ax21.legend(loc=1)
    ax22.legend(loc=1)
    ax23.legend(loc=1)
    fig.suptitle('%s_%s_%s'%(FC,SFC,HC))
    figlist.append(fig)
    
#    ParSWRC_2D=map(lambda x : x[None,:],ParSWRC)
#    ParKh_2D=map(lambda x : x[None,:],ParKh)
    
    #==============================================================================
    # Iteration
    #==============================================================================
    
    #t_request=[4000,80000] # list of forced t
    t_request=[] # list of forced t
    t_request=t_request+t_rain_L
    t_request=numpy.sort(numpy.unique(t_request))
    
    dt=[]       # List with dt (s)
    t=[t0]      # List with t (s)
    h=[h0]      # List with h (m)
    nit=[]      # List with n it
    q_BCtop=[]
    q_runoff=[]
    BC_top=[]
    h_SL=[0,]
    
    DoDiricheletAnyway = False
    dt_now=dt_init
    end_T=0
    iT=0
    end_IT=0
    
    BC_switch_prev=0
    h_SurfaceLayer=numpy.array([0.,])
    h_SurfaceLayer_MAX=10.0
    #h_SurfaceLayer_MAX=0.01
    h_SurfaceLayer_MIN=-1e5
    
    #while end_T==0:
    while (end_T==0) & (dt_now > dt_min):
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
        z_gwt=BC_GWT[numpy.array(t_rain_L)>t[iT]][0]
        
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
    #    if any(t[iT] == numpy.array(t_rain_L)):
    #        ResetdT=numpy.array([ True,])
        
    #    if ((iT>0) & (BC_switch_prev!=BC_switch) & end_IT==1)|(ResetdT):
        if ((iT>0) & (BC_switch_prev!=BC_switch)  & end_IT==1)|(ResetdT):
            if dt_now>dt_init:
                dt_now=dt_init
        #        dt_now=dt_BC
                t_now=t[iT]+dt_now
    
        # Override everything and use Hydrus time
    #    dt_now=dT_Hydrus[iT]    # use Hydrus dt
    #    t_now=t[iT]+dt_now     # use Hydrus dt
        
    
    #    print('Time: %1.1f / %1.1f - BC: %s - q_rain = %1.2e - h_SL = %1.2e' %(t_now,t_end,['Dirichelet','Neumann'][BC_switch[0]],q_rain,h_SurfaceLayer))
        print('Time: %1.2e / %1.2e - dT = %1.2e -  BC: %s - q_rain = %1.2e - h_SL = %1.2e' %(t_now,t_end,dt_now,['Dirichelet','Neumann'][BC_switch[0]],q_rain,h_SurfaceLayer))
        
        
        ## initialize h_prevIT
        
        if iT > 0:
            dt_prev=dt[iT-1]
            h_prevIT=h[iT]+dt_now/dt_prev*(h[iT]-h[iT-1]) ## Linear extrapolation of previous solution
        else:
            h_prevIT=copy.deepcopy(h[iT])
    #    h_prevIT=copy.deepcopy(h[iT])
        
        ## Experimental
        if (h_prevIT[0] > -1e-2) & (q_rain>=0) & (BC_switch==1):
            print('experimental stabilisation correction')
            h_prevIT[0] = -1e-2
        
        ## Start Iteration    
        
        nit_now=0
        end_IT=0
        
        while (end_IT==0) & (nit_now<=nit_term):
    #        print(nit_now)
            
            
            # GET P & F WITHOUT BOUNDARY CONDITIONS
            P=P_j1k(h_prevIT,z,dt_now,incl,ModelSWRC,ParSWRC,ModelKh,ParKh,ModelCap,ParCap)
            F=F_j1k(h_prevIT,h_prevT,z,dt_now,incl,ModelSWRC,ParSWRC,ModelKh,ParKh,ModelCap,ParCap)
    
            # TOP BOUNDARY CONDITION
    
            if (BC_switch) & (~DoDiricheletAnyway) :
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
            h_BCbottom=z_gwt-z[-1]
            if h_BCbottom>=0:
                P,F=BC_DiricheletBottom(P,F,h_prevIT,h_prevT,dt_now,h_BCbottom,z,ModelKh,ParKh)
            else:
                P,F=BC_FreeDrainageBottom(P,F,h_prevIT,h_prevT,dt_now,z,incl,ModelSWRC,ParSWRC,ModelKh,ParKh,ModelCap,ParCap)
    #        P,F=BC_FreeDrainageBottom(P,F,h_prevIT,h_prevT,dt_now,z,incl,ModelSWRC,ParSWRC,ModelKh,ParKh,ModelCap,ParCap)
    
    #        P,F=BC_NeumannBottom(P,F,h_prevIT,h_prevT,dt_now,q_bottom,z,incl,ModelSWRC,ParSWRC,ModelKh,ParKh,ModelCap,ParCap)
            
            
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
        
        if (h_now[0]>0) & BC_switch & (q_rain<0):
            end_IT = 0
            DoDiricheletAnyway=True
        
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
            q_BCtop.append(q_rain)
            q_runoff.append((h_runoff)/dt_now)
            
            # New time
            DoDiricheletAnyway=False
            iT+=1
            if (nit_now <= nit_incr_dt) & (dt_now<dt_max):
                print('Increasing time step')
                dt_now=dt_now*incr_dt
            elif (nit_now >= nit_decr_dt) & (dt_now>dt_min):
                print('Decreasing time step')
                dt_now=dt_now*decr_dt
                
        elif(dt_now>dt_min):
            if ~DoDiricheletAnyway:
                print('Decreasing time step drastically')
                dt_now=dt_now*dt_term
                print dt_now
        
        
        if dt_now < dt_min :
            print('Not converging, minimal time step reached, ending calculation')
            end_T=1
        
        if len(t)>maxHesitation_N:
            if t[-1] - t[-maxHesitation_N] < maxHesitation_T:
                print('Hesitating too long, ending calculation')
                end_T=1
        
        ## Test stability of the convergence
        # ?
        
    h=numpy.array(h)
    #theta=ModelSWRC(h,ParSWRC_2D)
    theta=numpy.zeros_like(h)
    for iT in range(len(t)):
        theta[iT,:] = ModelSWRC(h[iT,:],ParSWRC)
    t=numpy.array(t)
    dt=numpy.array(dt)
    nit=numpy.array(nit)
    H=h+(-z*numpy.cos(incl))
    q_BCtop=numpy.array(q_BCtop)
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
    
    #q_top_hydrus=TlevelData[:,3]
    
    q_top=Flux_Top(z,h,dt,incl,ModelSWRC,ParSWRC,ModelKh,ParKh)
    
    
    #q_bottom_hydrus=TlevelData[:,5]
    
    q_bottom=Flux_Bottom(z,h,dt,incl,ModelSWRC,ParSWRC,ModelKh,ParKh)
    
    
    #q_runoff_hydrus=TlevelData[:,14]
    
    #==============================================================================
    # Check
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
    ax=fig.add_subplot(211)
    ax.plot(nit,'k-',label='Replicate')
    ax.set_ylim([0,None])
    ax.set_title('number of iterations')
    ax=fig.add_subplot(212)
    ax.plot(dt,'k-',label='Replicate')
    ax.set_yscale('log')
    ax.set_title('dt')
    figlist.append(fig)
    
    fig=pylab.figure()
    ax=fig.add_subplot(111)
    ax.plot(t[1:],q_top,'k-',label='Replicate')
    ax.plot(t[1:],q_BCtop,'b-',label='BC')
    #ax.plot(t_hydrus,q_top_hydrus,'r.-',label='Hydrus')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('q_top (cm/s)')
    ax.legend()
    ax.set_title('q_top (cm/s)')
    figlist.append(fig)
    
    fig=pylab.figure()
    ax=fig.add_subplot(111)
    ax.plot(t[1:],q_top,'k-',label='Replicate')
    ax.plot(t[1:],q_BCtop,'b-',label='BC')
    #ax.plot(t_hydrus,q_top_hydrus,'r.-',label='Hydrus')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('q_top (cm/s)')
    ax.set_ylim([-2e-4,5e-5])
    ax.legend()
    ax.set_title('q_top (cm/s)')
    figlist.append(fig)
    
    fig=pylab.figure()
    ax=fig.add_subplot(111)
    ax.plot(t[1:],q_bottom,'k-',label='Replicate')
    #ax.plot(t_hydrus,q_bottom_hydrus,'r.-',label='Hydrus')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('q_bottom (cm/s)')
    ax.legend()
    ax.set_title('q_bottom (cm/s)')
    figlist.append(fig)
    
    
    fig=pylab.figure()
    ax=fig.add_subplot(111)
    ax.plot(t[1:],q_runoff,'k-',label='Replicate')
    #ax.plot(t_hydrus,q_runoff_hydrus,'r.-',label='Hydrus')
    ax.set_xlabel('Time (s)')
    ax.set_ylabel('q_runoff (cm/s)')
    ax.legend()
    ax.set_title('q_runoff (cm/s)')
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
    
#    nO=3
#    for iO in range(nO):
#        fig=pylab.figure(figsize=figsize1)
#        ax=fig.add_subplot(111)
#        iZ=ObsNode_Nodes[iO]-1
#        ax.plot(t,h[:,iZ],'k-',label='Replicate')
#        #ax.plot(ObsNode_t,ObsNode_L[iO][:,0],'r-',label='Hydrus')
#        ax.set_xlabel('Time')
#        ax.set_ylabel('Head')
#        ax.legend(loc=4)
#        #ax.set_xlim([0,.1])
#        #ax.set_ylim([-1000,0])
#        ax.set_title('Observation node z = %s cm' %z[iZ])
#        figlist.append(fig)
    
    
    fig=pylab.figure()
    ax=fig.add_subplot(111)
    ax.axhline(0,color='k')
    ax.plot(t[1:][BC_top=='Dirichelet'],h_SL[:-1][BC_top=='Dirichelet'],'r.',label='Dirichelet')
    ax.plot(t[1:][BC_top=='Neumann'],h_SL[:-1][BC_top=='Neumann'],'b.',label='Neumann')
    ax.plot(t,h[:,0],'k--')
    ax.set_xlabel('h Surface Layer (cm)')
    ax.set_ylabel('Head')
    ax.set_ylim([-1e-2,1e-1])
    ax.legend()
    figlist.append(fig)



    pylab.close('all')
    
    #==============================================================================
    # Save
    #==============================================================================
    
    ## Save data
    import sys
    basename=os.path.basename(sys.argv[0])[:-3]
    
    #Root='C:\Users\Administrator\Documents\Arbeid\Nevenactiviteiten\Hydrus_MakeItRain'
    Root='/home/supersoil/Documents/Jan/AlbonHydrus_DoTheWalk'
#    Root='/home/supersoil/Documents/Jan/AlbonHydrus_DoTheWalk/FirstTest'
    postfix='_%s_%s_%s_%s'%(FC,SFC,HC,K)
    
    npzname=os.path.join(Root,basename+postfix)
    numpy.savez(npzname,h=h,theta=theta,t=t,dt=dt,nit=nit,H=H,z=z,
                q_runoff=q_runoff,q_top=q_top,q_bottom=q_bottom,
                h_SL=h_SL,BC_top=BC_top,q_BCtop=q_BCtop,h_BCbottom=h_BCbottom,
                ModelSWRC=ModelSWRC,ParSWRC=ParSWRC,ModelKh=ModelKh,ParKh=ParKh,ModelCap=ModelCap,ParCap=ParCap)
    
    
    ## save plots to pdf
    
    pdfname = os.path.join(Root, basename + postfix + '.pdf')
    pp = PdfPages(pdfname)
    for fig in figlist:
        pp.savefig(fig)
    
    pp.close()
    
    extension='.png'
    for iF in range(len(figlist)):
        fig=figlist[iF]
        figname = os.path.join(Root, basename+'_%s'%iF + postfix + extension)
        fig.savefig(figname, bbox_inches='tight')
    
    #pylab.show()
    #os.startfile(pdfname)
