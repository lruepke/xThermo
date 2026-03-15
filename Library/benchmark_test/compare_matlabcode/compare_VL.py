# -*- coding: utf-8 -*-
import os
import numpy as np
import time
import linecache
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import patches
import matplotlib.ticker as ticker
from matplotlib.ticker import MultipleLocator
from tabulate import tabulate
import copy
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
# 3d plot
import helpfunc
mpl.rcParams['font.family'] = 'Arial'  # default font family
mpl.rcParams['mathtext.fontset'] = 'cm'  # font for math
fmt_figs = ['pdf']  # ['svg','pdf']
figpath = '.'
phaseBoundary='VL'
result_path='VL'
def savefig(figname):
    for fmt_fig in fmt_figs:
        figname_full = '%s/%s.%s' % (figpath, figname, fmt_fig)
        plt.savefig(figname_full, bbox_inches='tight')
        print('figure saved: ', figname_full)
compare = lambda a,b : float(str('%.6e'%(a)))-float(str('%.6e'%(b)))
# Import package of xThermal
from xThermal import H2O
from xThermal import NaCl
from xThermal import H2ONaCl
sw_84 = H2ONaCl.cH2ONaCl("IAPS84")
sw_95 = H2ONaCl.cH2ONaCl("IAPWS95")

def plot_3d_matlab():
    fig=plt.figure()
    ax = fig.add_subplot(111,projection='3d',facecolor='None')
    helpfunc.set_axis_diagram_3D(ax)
    # read result calculated by matlab code
    TT=np.loadtxt('%s/TT_%s.txt'%(result_path,phaseBoundary)) #K
    PP=np.loadtxt('%s/PP_%s.txt'%(result_path,phaseBoundary)) #Pa
    # plot
    for fname,color in zip(['Xl','Xv'],['green','orange']):
        XX=np.loadtxt('%s/%s_%s.txt'%(result_path,fname,phaseBoundary))
        ax.plot_wireframe(XX*100,TT-273.15,PP/1E5, color=color,lw=0.1)
    savefig('PhaseBoundary_%s_3D'%(phaseBoundary))
def contourf_phase(ax,TT,pp,phase,phase_name,ax_cb=None):
    cmap = plt.get_cmap("Dark2")
    # customize cmap
    colors=list(copy.deepcopy(cmap.colors))
    colors[0:8]=['lightblue','red','lightgreen','lightgray','violet','yellow','lightcyan','lightcyan']
    cmap.colors=tuple(colors)
    CS=ax.contourf(TT,pp,phase, cmap=cmap,vmin=phase.min()-0.5, vmax=phase.max()+0.5, levels=np.linspace(phase.min()-0.5,phase.max()+0.5,len(phase_name)+1))
    if(ax_cb is None): ax_cb = ax.inset_axes([0,1.03,1,0.05])
    cb=plt.colorbar(CS, cax=ax_cb, orientation='horizontal',ticklocation='top',ticks=np.arange(phase.min(),phase.max()+1))
    cb.ax.set_xticklabels(phase_name)
    return CS,ax_cb,cb
def plot_props_VLH(sw,water,mmc5='../Driesner2007a/1-s2.0-S0016703707002943-mmc5.txt',mmc3='../Driesner2007b/1-s2.0-S0016703707002955-mmc3.txt'):
    if(not os.path.exists(mmc5)):
        print('Please set correct mmc1 file path: %s'%(mmc5))
        exit()
    data=np.loadtxt(mmc5, skiprows=7)
    T0_5,P0_5,XV0_5,XL0_5=data[:,0]+273.15,data[:,1]*1E5,data[:,2],data[:,3]
    if(not os.path.exists(mmc3)):
        print('Please set correct mmc1 file path: %s'%(mmc3))
        exit()
    data=np.loadtxt(mmc3, skiprows=5)
    T0,P0,XV0,rhoV0,hV0,XL0,rhoL0,hL0=data[:,0]+273.15,data[:,1]*1E5,data[:,2],data[:,3],data[:,4],data[:,5],data[:,6],data[:,7]
    # calculate
    T = np.linspace(H2ONaCl.T_MIN_VLH, H2ONaCl.T_MAX_VLH, 500)
    P = np.array(sw.P_VLH(T))
    # X_haliteLiquidus=np.array(sw.X_HaliteLiquidus(T,P))
    XL_,XV_ = np.array(sw.X_VLH(T,P))
    rhoV_, rhoL_ = np.array(sw.Rho_phase(T, P, XV_, H2ONaCl.Vapor)), np.array(sw.Rho_phase(T, P, XL_, H2ONaCl.Liquid))
    hV_, hL_     = np.array(sw.H_phase(T, P, XV_, H2ONaCl.Vapor)), np.array(sw.H_phase(T, P, XL_, H2ONaCl.Liquid))

    # plot
    fig,axes = plt.subplots(1,6,figsize=(30,4),gridspec_kw={'wspace':0.1},sharey=True)
    # 1. liquid salinity
    ax=axes[0]
    line=helpfunc.plot_coloredline(ax,np.array(sw.Wt2Mol(XL_)), T-273.15, P/1E5,cmap='rainbow')
    ax_cb = ax.inset_axes([0.8,0.15,0.02,0.6])
    fig.colorbar(line,cax = ax_cb,label='Pressure (bar)')
    ax.plot(XL0_5, T0_5-273.15,color='gray',ls='dashed', marker='.',markersize=5,mec='w',mfc='k',markeredgewidth=0.3,label='Driesner & Heinrich(2007)')
    ax.xaxis.set_major_locator(MultipleLocator(0.2))
    ax.xaxis.set_minor_locator(MultipleLocator(0.04))
    ax.set_xlabel('Liquid composition X$_{\mathregular{NaCl}}$ (mole fraction)')
    # 2. vapor salinity
    ax=axes[1]
    helpfunc.plot_coloredline(ax,np.array(sw.Wt2Mol(XV_)), T-273.15, P/1E5,cmap='rainbow')
    ax.set_xscale('log')
    ax.plot(XV0_5, T0_5-273.15,color='gray',ls='dashed', marker='.',markersize=5,mec='w',mfc='k',markeredgewidth=0.3,label='Driesner & Heinrich(2007)')
    ax.set_xlabel('Vapor composition X$_{\mathregular{NaCl}}$ (mole fraction)')
    # 3. liquid density
    ax=axes[2]
    helpfunc.plot_coloredline(ax,rhoL_, T-273.15, P/1E5,cmap='rainbow')
    ax.plot(rhoL0, T0-273.15,color='gray',ls='dashed', marker='.',markersize=5,mec='w',mfc='k',markeredgewidth=0.3,label='Driesner(2007)')
    ax.xaxis.set_major_locator(MultipleLocator(100))
    ax.xaxis.set_minor_locator(MultipleLocator(20))
    ax.set_xlabel('Liquid density (kg/m$^{\mathregular{3}}$)')
    # 4. Vapor density
    ax=axes[3]
    helpfunc.plot_coloredline(ax,rhoV_, T-273.15, P/1E5,cmap='rainbow')
    ax.plot(rhoV0, T0-273.15,color='gray',ls='dashed', marker='.',markersize=5,mec='w',mfc='k',markeredgewidth=0.3,label='Driesner(2007)')
    ax.xaxis.set_major_locator(MultipleLocator(20))
    ax.xaxis.set_minor_locator(MultipleLocator(4))
    ax.set_xlabel('Vapor density (kg/m$^{\mathregular{3}}$)')
    # 5. liquid enthalpy
    ax=axes[4]
    helpfunc.plot_coloredline(ax,hL_/1E6, T-273.15, P/1E5,cmap='rainbow')
    # ax.plot(hL_/1E6, T-273.15,'.')
    ax.plot(hL0/1E6, T0-273.15,color='gray',ls='dashed', marker='.',markersize=5,mec='w',mfc='k',markeredgewidth=0.3,label='Driesner(2007)')
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.set_xlabel('Liquid enthalpy (MJ/kg)')
    # 6. vapor enthalpy
    ax=axes[5]
    helpfunc.plot_coloredline(ax,hV_/1E6, T-273.15, P/1E5,cmap='rainbow')
    ax.plot(hV0/1E6, T0-273.15,color='gray',ls='dashed', marker='.',markersize=5,mec='w',mfc='k',markeredgewidth=0.3,label='Driesner(2007)')
    ax.xaxis.set_major_locator(MultipleLocator(0.5))
    ax.xaxis.set_minor_locator(MultipleLocator(0.1))
    ax.set_xlabel('Vapor enthalpy (MJ/kg)')

    for ax in axes:
        ax.yaxis.set_major_locator(MultipleLocator(100))
        ax.yaxis.set_minor_locator(MultipleLocator(20))
        ax.grid(which='major',lw=0.04,color='k')
        ax.grid(which='minor',lw=0.04,color='gray')
        ax.axhline(H2ONaCl.T_MIN_VLH-273.15,label='T=%.2f $^{\circ}$C, p=1 bar'%(H2ONaCl.T_MIN_VLH-273.15),color='r',ls='dotted')
        ax.legend()
        # ax.set_ylim(300,T.max()-273.15 + 20)
    axes[0].set_ylabel('Temperature ($^{\circ}$C)')
    savefig('VLH_props_%s'%(sw.name_backend()))

    # let's checkout what happens in the (low T, low p) and (high T, low p) region for the liquid enthalpy
    q1,q2 = np.array(sw.q1q2_Tstar_H(P, XL_))
    q1_v,q2_v = np.array(sw.q1q2_Tstar_H(P, XV_))
    Tstar_H,Tstar_H_v = q1 + q2*(T-273.15) + 273.15, q1_v + q2_v*(T-273.15) + 273.15
    n1,n2 = np.array(sw.n1n2_Tstar_V(P, XL_))
    n1_v,n2_v = np.array(sw.n1n2_Tstar_V(P, XV_))
    Tstar_V,Tstar_V_v = n1 + n2*(T-273.15) + 273.15, n1_v + n2_v*(T-273.15) + 273.15
    T_water,p_water=np.linspace(np.array([T.min(),Tstar_H.min(),Tstar_H_v.min(),Tstar_V.min(),Tstar_V_v.min()]).min(),np.array([T.max(),Tstar_H.max(),Tstar_H_v.max(),Tstar_V.max(),Tstar_V_v.max()]).max()+100,1000), np.linspace(np.log10(P.min()),np.log10(P.max())+0.1,500)
    TT,pp=np.meshgrid(T_water,10**p_water)
    phase,rho,h = np.zeros_like(TT),np.zeros_like(TT),np.zeros_like(TT)
    for i in range(0,TT.shape[0]):
        for j in range(0,TT.shape[1]):
            water.UpdateState_TP(TT[i][j], pp[i][j])
            phase[i][j]=water.phase()
            rho[i][j]=water.rhomass()
            h[i][j]=water.hmass()
    # get phase names
    phase_unique = np.sort(np.unique(phase))
    phase_name = ['']*len(phase_unique)
    for i,phase0 in enumerate(phase_unique):
        phase[phase==phase0]=i+phase_unique.max()+10
        phase_name[i]=water.phase_name(int(phase0))
    fig,axes=plt.subplots(1,3,figsize=(21,5),gridspec_kw={'wspace':0.05},sharey=True)
    axes[0].set_ylabel('Pressure (bar)')
    l_L,l_V=[],[]
    for ax,prop, cmap,label in zip(axes, [phase, rho, h/1E6],['Paired','YlGnBu_r','RdBu'],['Phase','Density (kg/m$^{\mathregular{3}}$)','Specific enthalpy (MJ/kg)']):
        CS=[]
        if(ax==axes[0]):
            CS,ax_cb,cb=contourf_phase(ax,TT-273.15,pp/1E5,prop,phase_name,ax.inset_axes([0,1.03,1,0.03]))
        else:
            CS=ax.contourf(TT-273.15,pp/1E5,prop,levels=50,cmap=cmap)
            ax_cb=ax.inset_axes([0,1.02,1,0.03])
            plt.colorbar(CS,cax=ax_cb,label=label,orientation='horizontal')
            ax_cb.xaxis.set_label_position('top')
            ax_cb.xaxis.set_ticks_position('top')
        l_L,=ax.plot(Tstar_H-273.15, P/1E5,lw=2,label='$T^*_h$: liquid')
        l_V,=ax.plot(Tstar_H_v-273.15, P/1E5,lw=2,label='$T^*_h$: vapor')
        ax.plot(Tstar_V-273.15, P/1E5,lw=2,ls='dashed',label='$T^*_V$: liquid')
        ax.plot(Tstar_V_v-273.15, P/1E5,lw=2,ls='dashed',label='$T^*_V$: vapor')
        ax.plot(T-273.15, P/1E5,lw=0.8,marker='.',markevery=20,mec='w',mew=0.5,ms=10,label='$T_{VLH}$')
        ax.set_yscale('log')
        ax.set_xscale('log')
        ax.set_xlabel('Temperature ($^{\circ}$C)')
    axes[0].legend(loc='center left') #,bbox_to_anchor=[1.01,0]
    savefig('Tstar_VLH_%s'%(sw.name_backend()))
def plot_err(ax,x,y,data,label='',cmap='rainbow',scale_data='log',vmin=1E-6,vmax=1,s=1):
    # plot difference between xThermal and Driesner(2007b)
    norm = None
    if(scale_data=='log'):
        norm = mpl.colors.LogNorm(vmin=vmin, vmax=vmax)
    CS=ax.scatter(x, y,c=data,cmap=cmap,norm=norm,s=s,zorder=3)
    if(label!=''):
        ax_cb=ax.inset_axes([0,1.02,1,0.05])
        plt.colorbar(CS,cax=ax_cb,label=label,orientation='horizontal',extend='both')
        ax_cb.xaxis.set_label_position('top')
        ax_cb.xaxis.set_ticks_position('top')
def benchmark_VL(sw,mmc4='../Driesner2007b/1-s2.0-S0016703707002955-mmc4.txt'):
    # compare
    if(not os.path.exists(mmc4)):
        print('Please set correct mmc1 file path: %s'%(mmc4))
        exit()
    data=np.loadtxt(mmc4, skiprows=7)
    T0,P0,XV0,rhoV0,hV0,XL0,rhoL0,hL0=data[:,0]+273.15,data[:,1]*1E5,data[:,2],data[:,3],data[:,4],data[:,5],data[:,6],data[:,7]
    XL0_wt,XV0_wt = np.array(sw.Mol2Wt(XL0))*100, np.array(sw.Mol2Wt(XV0))*100
    # ind=(T0>sw.Tmin_VLH())
    # # only compare the result in valid range of pressure: >1bar
    # T0,P0,XV0,rhoV0,hV0,XL0,rhoL0,hL0 = T0[ind],P0[ind],XV0[ind],rhoV0[ind],hV0[ind],XL0[ind],rhoL0[ind],hL0[ind]
    # 1. calculate halite liquidus
    XL_,XV_ = np.array(sw.XL_VL(T0,P0)), np.array(sw.XV_VL(T0,P0))
    XL_mol_,XV_mol_ = np.array(sw.Wt2Mol(XL_)), np.array(sw.Wt2Mol(XV_))
    # 2. calculate saturated liquid density and vapor density
    rhoV_, rhoL_ = np.array(sw.Rho_phase(T0, P0, XV_, H2ONaCl.Vapor)), np.array(sw.Rho_phase(T0, P0, XL_, H2ONaCl.Liquid))
    hV_, hL_     = np.array(sw.H_phase(T0, P0, XV_, H2ONaCl.Vapor)), np.array(sw.H_phase(T0, P0, XL_, H2ONaCl.Liquid))
    # compare result dict
    Data0 = {'XV':XV0,'rhoV':rhoV0,'hV':hV0,'XL':XL0,'rhoL':rhoL0,'hL':hL0}
    Data_ = {'XV':XV_mol_,'rhoV':rhoV_,'hV':hV_,'XL':XL_mol_,'rhoL':rhoL_,'hL':hL_}
    Err,RErr={},{}
    for key in Data0.keys(): Err[key],RErr[key] = Data0[key]-Data_[key], np.abs(Data0[key]-Data_[key])/(Data0[key])*100.0
    # print to file
    fpout = open('%s/mmc4_%s.csv'%(result_path,sw.name_backend()),'w')
    fpout.write('T[C],P[bar],XV(Driesner)[mol],XV(xThermal)[mol],XV(diff)[mol],RhoV(Driesner)[kg/m3],RhoV(xThermal),RhoV(err),HV(Driesner)[J/kg],HV(xThermal),HV(err),XL(Driesner)[mol],XL(xThermal)[mol],XL(diff)[mol],RhoL(Driesner)[kg/m3],RhoL(xThermal),RhoL(err),HL(Driesner)[J/kg],HL(xThermal),HL(err)\n')
    for i in range(0,len(T0)):
        fpout.write('%.6e,%.6e'%(T0[i]-273.15,P0[i]/1E5))
        for key in Data0.keys():
            fpout.write(',%.6e,%.6e,%.6e'%(Data0[key][i], Data_[key][i],compare(Data0[key][i],Data_[key][i])))
        fpout.write('\n')
    fpout.close()

    # plot difference
    T_crit=np.linspace(T0.min(), T0.max(), 200)
    P_crit,X_crit=np.array(sw.P_X_Critical(T_crit))
    fig,axes2=plt.subplots(2,6,figsize=(30,10),gridspec_kw={'wspace':0.05,'hspace':0.05})
    # X-T space
    axes=axes2[0,:]
    for ax in [axes[1],axes[3],axes[5]]: ax.set_xscale('log')
    for ax in axes[0:2]: plot_err(ax, np.append(XV0_wt,XL0_wt), np.append(T0-273.15, T0-273.15), np.append(RErr['XV'],RErr['XL']),'Relative difference: Composition (%)')
    for ax in axes[2:4]: plot_err(ax, np.append(XV0_wt,XL0_wt), np.append(T0-273.15, T0-273.15), np.append(RErr['rhoV'],RErr['rhoL']),'Relative difference: Density (%)')
    for ax in axes[4:]: plot_err(ax, np.append(XV0_wt,XL0_wt), np.append(T0-273.15, T0-273.15), np.append(RErr['hV'],RErr['hL']),'Relative difference: Specific enthalpy (%)')
    for ax in axes: ax.plot(X_crit*100, T_crit-273.15,color='k',lw=2)
    axes[0].set_ylabel('Temperature ($^{\circ}$C)')
    for ax in axes[1:]: ax.yaxis.set_ticklabels([])
    # X-p space
    axes=axes2[1,:]
    for ax in [axes[1],axes[3],axes[5]]: ax.set_xscale('log')
    for ax in axes[0:2]: plot_err(ax, np.append(XV0_wt,XL0_wt), np.append(P0/1E5,P0/1E5), np.append(RErr['XV'],RErr['XL']))
    for ax in axes[2:4]: plot_err(ax, np.append(XV0_wt,XL0_wt), np.append(P0/1E5,P0/1E5), np.append(RErr['rhoV'],RErr['rhoL']))
    for ax in axes[4:]: plot_err(ax, np.append(XV0_wt,XL0_wt), np.append(P0/1E5,P0/1E5), np.append(RErr['hV'],RErr['hL']))
    for ax in axes: ax.plot(X_crit*100, P_crit/1E5,color='k',lw=2)
    axes[0].set_ylabel('Pressure (bar)')
    for ax in axes[1:]: ax.yaxis.set_ticklabels([])
    for i in range(axes2.shape[1]):
        for ax in axes2[:,i]:
            ax.grid(which='major',lw=0.04,color='k')
            ax.grid(which='minor',lw=0.04,color='gray')
        axes2[1,i].set_xlabel('wt.% NaCl')
    for ax in axes2[0,:]: ax.xaxis.set_ticklabels([])
    savefig('diff_VL')
    # statistics of the difference
    table=[]
    for key,name in zip(list(Err.keys()),['XV (mole fraction)','RhoV (kg/m3)','HV (J/kg)','XL (mole fraction)','RhoL (kg/m3)','HL (J/kg)']):
        RErr[key] = RErr[key][~(np.isnan(RErr[key]) | np.isinf(RErr[key]))]
        table.append([name,Err[key].min(),Err[key].max(),RErr[key].min(),RErr[key].max()])
    print(tabulate(table, headers=['Critical property', 'Err. Min', 'Err. Max','RE. Min(%)','RE. Max(%)'],numalign="right",floatfmt=".6f"))

def calculate(sw=sw_84):
    # read result calculated by matlab code
    TT=np.loadtxt('%s/TT_%s.txt'%(result_path,phaseBoundary)) #K
    PP=np.loadtxt('%s/PP_%s.txt'%(result_path,phaseBoundary)) #Pa
    XL_,XV_ = np.zeros_like(TT), np.zeros_like(TT)
    for i in range(0,TT.shape[0]):
        XL_[i,:],XV_[i,:] = np.array(sw.XL_VL(TT[i,:],PP[i,:])), np.array(sw.XV_VL(TT[i,:],PP[i,:]))
    # save result
    for fname,XX in zip(['Xl','Xv'],[XL_,XV_]):
        fpout=open('%s/%s_%s_xThermal.txt'%(result_path,fname,phaseBoundary),'w')
        for i in range(0,TT.shape[0]):
            for j in range(0,TT.shape[1]):
                fpout.write('%.8E '%(XX[i][j]))
            fpout.write('\n')
        fpout.close()
    
    fig=plt.figure()
    ax = fig.add_subplot(111,projection='3d',facecolor='None')
    helpfunc.set_axis_diagram_3D(ax)
    for XX,color in zip([XL_, XV_],['green','orange']):
        ax.plot_wireframe(XX*100,TT-273.15,PP/1E5, color=color,lw=0.1)
    savefig('PhaseBoundary_%s_3D_xThermal'%(phaseBoundary))
def plot_diff():
    Err, RErr = {},{}
    TT=np.loadtxt('%s/TT_%s.txt'%(result_path,phaseBoundary)) #K
    PP=np.loadtxt('%s/PP_%s.txt'%(result_path,phaseBoundary)) #Pa
    for dataname,color in zip(['Xl','Xv'],['green','orange']):
        XX=np.loadtxt('%s/%s_%s.txt'%(result_path,dataname,phaseBoundary))
        XX_=np.loadtxt('%s/%s_%s_xThermal.txt'%(result_path,dataname,phaseBoundary))
        Err[dataname]=XX-XX_
        RErr[dataname]=Err[dataname]/XX*100
    fig,axes=plt.subplots(1,2,figsize=(14,6),gridspec_kw={'wspace':0.1,'hspace':0.1})
    norm = mpl.colors.LogNorm()
    norm = mpl.colors.SymLogNorm(linthresh=1E-9, linscale=1, vmin=-1, vmax=1, base=10)
    for dataname,ax in zip(['Xl','Xv'],axes):
        ind_nan=(np.isnan(Err[dataname]) | np.isinf(Err[dataname]) | (Err[dataname]==0))
        Err[dataname][ind_nan]=0.001
        CS=ax.contourf(TT-273.15,PP/1E5,np.log10(np.abs(Err[dataname])),cmap='rainbow',levels=50)
        ax_cb=ax.inset_axes([0,1.02,1,0.05])
        plt.colorbar(CS,cax=ax_cb,label=dataname,orientation='horizontal',extend='both')
        ax_cb.xaxis.set_label_position('top')
        ax_cb.xaxis.set_ticks_position('top')
        # ax_cb.set_xscale('log')
        # print(np.log10(np.abs(Err[dataname])).min())
    savefig('diff_%s'%(phaseBoundary))
    # statistics of the difference
    table=[]
    for key,name in zip(list(Err.keys()),['XV (mole fraction)','XL (mole fraction)']):
        RErr[key] = RErr[key][~(np.isnan(RErr[key]) | np.isinf(RErr[key]))]
        ind_nan=(np.isnan(Err[key]) | np.isinf(Err[key]))
        Err[key] = Err[key][~ind_nan]
        T,P=TT[~ind_nan],PP[~ind_nan]
        ind_max=(Err[key]==Err[key].max())
        table.append([name,Err[key].min(),Err[key].max()])
    print(tabulate(table, headers=['Critical property', 'Err. Min', 'Err. Max'],numalign="right",floatfmt=".6f"))
    

# plot_3d_matlab()
# calculate()
plot_diff()