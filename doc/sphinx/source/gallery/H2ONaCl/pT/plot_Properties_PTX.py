# -*- coding: utf-8 -*-
"""
6. Properties
==================================================
.. include:: /include.rst_
Calculate and plot phase diagram.
"""
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
from matplotlib.patches import Patch
import copy
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)
# 3d plot
import helpfunc
mpl.rcParams['font.family'] = 'Arial'  # default font family
mpl.rcParams['mathtext.fontset'] = 'cm'  # font for math
dpi=100
fmt_figs = ['pdf']  # ['svg','pdf']
figpath = '.'
result_path='../../../gallery_H2ONaCl/pT'
def savefig(figname):
    for fmt_fig in fmt_figs:
        figname_full = '%s/%s.%s' % (figpath, figname, fmt_fig)
        plt.savefig(figname_full, bbox_inches='tight',dpi=dpi)
        print('figure saved: ', figname_full)
compare = lambda a,b : float(str('%.6e'%(a)))-float(str('%.6e'%(b)))
# Import package of xThermal
from xThermal import H2O
from xThermal import NaCl
from xThermal import H2ONaCl
sw_84 = H2ONaCl.cH2ONaCl("IAPS84")
sw_95 = H2ONaCl.cH2ONaCl("IAPWS95")
#sw_95_CoolProp = H2ONaCl.cH2ONaCl("IAPWS95_CoolProp")

def plot_Phase(ax, XX,YY,Phase,sw, cmap="Dark2"):
    cmap = plt.get_cmap(cmap)
    # customize cmap
    phase_unique = np.sort(np.unique(Phase))
    phase_name = ['']*len(phase_unique)
    for i,phase0 in enumerate(phase_unique):
        Phase[Phase==phase0]=i+phase_unique.max()+10
        phase_name[i]=sw.phase_name(int(phase0))
    colors=list(copy.deepcopy(cmap.colors))
    colors[0:8]=['red','lightblue','lightgreen','lightgray','violet','yellow','lightcyan','lightcyan']
    cmap.colors=tuple(colors)
    CS=ax.contourf(XX,YY, Phase,cmap=cmap,vmin=Phase.min()-0.5, vmax=Phase.max()+0.5, levels=np.linspace(Phase.min()-0.5,Phase.max()+0.5,len(phase_name)+1))
    ax_cb = ax.inset_axes([0,1.01,1,0.05])
    cb=plt.colorbar(CS, cax=ax_cb, orientation='horizontal',ticklocation='top',ticks=np.arange(Phase.min(),Phase.max()+1))
    cb.ax.set_xticklabels(phase_name)

def plot_prop(ax, XX,YY,prop, prop_name='prop', prop_unit='unit', cmap='rainbow',levels=50):
    CS=ax.contourf(XX,YY,prop, levels=levels, cmap=cmap)
    ax_cb = ax.inset_axes([0,1.01,1,0.05])
    cb=plt.colorbar(CS, cax=ax_cb, orientation='horizontal',ticklocation='top',label='%s (%s)'%(prop_name, prop_unit))

constT, constP, constX = 400+273.15, 250E5, 0.1

# %%
# Isobaric section
# --------------------------

def plot_props_P0(P0, sw):
    T = np.linspace(sw.Tmin(),sw.Tmax(),150)
    X = np.linspace(1E-6,1,150)
    TT,XX = np.meshgrid(T,X)
    PP = TT*0 + P0
    state = sw.UpdateState_TPX(TT.reshape(-1,), PP.reshape(-1,), XX.reshape(-1,))
    # plot
    fig, axes=plt.subplots(2,4,figsize=(28,12),sharey=True, sharex=True,gridspec_kw={'wspace':0.1,'hspace':0.2})
    Rho = np.array(state.Rho).reshape(TT.shape)
    Phase = np.array(state.phase).reshape(TT.shape)
    H = np.array(state.H).reshape(TT.shape)
    Cp = np.array(state.Cp).reshape(TT.shape)
    Mu = np.array(state.Mu).reshape(TT.shape)
    Mu_l = np.array(state.Mu_l).reshape(TT.shape)
    Mu_v = np.array(state.Mu_v).reshape(TT.shape)

    xx,yy = XX*100, TT-273.15
    # 1. phase
    ax = axes[0][0]
    plot_Phase(ax, xx,yy, Phase, sw)

    # 2. Rho
    ax = axes[0][1]
    plot_prop(ax, xx,yy, Rho,prop_name="Density", prop_unit="kg/m$^\mathregular{3}$")
    # 3. H
    ax = axes[0][2]
    plot_prop(ax, xx,yy, H/1E6,prop_name="Specific enthalpy", prop_unit="MJ/kg",cmap='GnBu')
    # 3. Cp
    ax = axes[0][3]
    plot_prop(ax, xx,yy, Cp,prop_name="Isobaric specific heat", prop_unit="J/kg/K",cmap='YlOrRd')
    # 1. Mu_l
    ax=axes[1][0]
    plot_prop(ax, xx,yy, np.log10(Mu_l),prop_name="log$_{10}$ $\mu_l$", prop_unit="$\mu Pa \cdot s$",cmap='GnBu')

    # 2. Mu_v
    ax=axes[1][1]
    plot_prop(ax, xx,yy, np.log10(Mu_v),prop_name="log$_{10}$ $\mu_v$", prop_unit="$\mu Pa \cdot s$",cmap='BuGn')

    # kappa: Isothermal compressibility
    ax=axes[1][2]
    kappa = np.array(state.IsothermalCompressibility).reshape(TT.shape)
    plot_prop(ax, xx,yy, np.log10(kappa), prop_name="log$_{10}$ $\kappa$", prop_unit="Pa$^{-1}$",cmap='BuGn')

    # beta: isobaric expansivity
    ax=axes[1][3]
    beta = np.array(state.IsobaricExpansivity).reshape(TT.shape)
    plot_prop(ax, xx,yy, np.log10((beta)), prop_name="log$_{10}$ $\\beta$", prop_unit="T$^{-1}$",cmap='BuGn')
    print(beta.min(),beta.max())
    for ax in axes.reshape(-1,):
        ax.axvline(constX*100, ls='dashed', color='k')
        ax.axhline(constT - 273.15, ls='dashdot', color='k')

    # set axis
    for ax in axes[1,:]: ax.set_xlabel('Bulk salinity (wt.% NaCl)')
    for ax in axes[:,0]: ax.set_ylabel('Temperature ($^{\circ}$C)')

    savefig('Props_P%.0fbar'%(P0/1E5))

plot_props_P0(constP, sw_84)

# %%
# Isothermal section
# --------------------------

def plot_props_T0(T0, sw):
    X = np.linspace(1E-6,1,400)
    P = np.linspace(1, 1000, 400)*1E5
    PP,XX = np.meshgrid(P,X)
    TT = PP*0 + T0
    state = sw.UpdateState_TPX(TT.reshape(-1,), PP.reshape(-1,), XX.reshape(-1,))
    # plot
    fig, axes=plt.subplots(2,4,figsize=(28,12),sharey=True, sharex=True,gridspec_kw={'wspace':0.1,'hspace':0.2})
    Rho = np.array(state.Rho).reshape(TT.shape)
    Phase = np.array(state.phase).reshape(TT.shape)
    H = np.array(state.H).reshape(TT.shape)
    Cp = np.array(state.Cp).reshape(TT.shape)
    Mu = np.array(state.Mu).reshape(TT.shape)
    Mu_l = np.array(state.Mu_l).reshape(TT.shape)
    Mu_v = np.array(state.Mu_v).reshape(TT.shape)

    xx,yy = XX*100, PP/1E5
    # 1. phase
    ax = axes[0][0]
    plot_Phase(ax, xx,yy, Phase, sw)
    # 2. Rho
    ax = axes[0][1]
    plot_prop(ax, xx,yy, Rho,prop_name="Density", prop_unit="kg/m$^\mathregular{3}$")
    # 3. H
    ax = axes[0][2]
    plot_prop(ax, xx,yy, H/1E6,prop_name="Specific enthalpy", prop_unit="MJ/kg",cmap='GnBu')
    # 3. Cp
    ax = axes[0][3]
    plot_prop(ax, xx,yy, Cp,prop_name="Isobaric specific heat", prop_unit="J/kg/K",cmap='YlOrRd')
    # 1. Mu_l
    ax=axes[1][0]
    plot_prop(ax, xx,yy, np.log10(Mu_l),prop_name="log$_{10}$ $\mu_l$", prop_unit="$\mu Pa \cdot s$",cmap='GnBu')
    # ax.grid()
    # 2. Mu_v
    ax=axes[1][1]
    plot_prop(ax, xx,yy, np.log10(Mu_v),prop_name="log$_{10}$ $\mu_v$", prop_unit="$\mu Pa \cdot s$",cmap='BuGn')

    # kappa: Isothermal compressibility
    ax=axes[1][2]
    kappa = np.array(state.IsothermalCompressibility).reshape(TT.shape)
    plot_prop(ax, xx,yy, np.log10(kappa), prop_name="log$_{10}$ $\kappa$", prop_unit="Pa$^{-1}$",cmap='BuGn')

    # beta: isobaric expansivity
    ax=axes[1][3]
    beta = np.array(state.IsobaricExpansivity).reshape(TT.shape)
    plot_prop(ax, xx,yy, np.log10((beta)), prop_name="log$_{10}$ $\\beta$", prop_unit="T$^{-1}$",cmap='BuGn')
    print(beta.min(),beta.max())
    for ax in axes.reshape(-1,):
        ax.axvline(constX*100, ls='dashed', color='k')
        ax.axhline(constP/1E5, ls='dashdot', color='k')

    # set axis
    for ax in axes[1,:]: ax.set_xlabel('Bulk salinity (wt.% NaCl)')
    for ax in axes[:,0]: ax.set_ylabel('Pressure (bar)')

    savefig('Props_T%.0fC'%(T0-273.15))

plot_props_T0(constT, sw_84)

# %%
# Constant X
# --------------------------
def plot_props_X0(X0,sw):
    dT = 1
    T = np.linspace(sw.Tmin(),sw.Tmax()-dT,400)
    P = np.linspace(1E5, 1000E5, 400)
    # # close to critical curve
    # T = np.linspace(550+273.15, 650+273.15, 100)
    # P = np.linspace(800E5, 1000E5, 100)

    TT,PP = np.meshgrid(T,P)
    XX = TT*0 + X0
    state = sw.UpdateState_TPX(TT.reshape(-1,), PP.reshape(-1,), XX.reshape(-1,))
    # state2 = sw.UpdateState_TPX(TT.reshape(-1,)+dT, PP.reshape(-1,), XX.reshape(-1,))
    # plot
    fig, axes=plt.subplots(2,4,figsize=(28,12),sharey=True, sharex=True,gridspec_kw={'wspace':0.1,'hspace':0.2})
    Rho = np.array(state.Rho).reshape(TT.shape)
    Phase = np.array(state.phase).reshape(TT.shape)
    H = np.array(state.H).reshape(TT.shape)
    # H2 = np.array(state2.H).reshape(TT.shape)/1E6
    # Cp2 = (H2-H)/dT
    Cp = np.array(state.Cp).reshape(TT.shape)
    Mu = np.array(state.Mu).reshape(TT.shape)
    Mu_l = np.array(state.Mu_l).reshape(TT.shape)
    Mu_v = np.array(state.Mu_v).reshape(TT.shape)

    xx,yy = TT-273.15,PP/1E5
    # 1. phase
    ax = axes[0][0]
    plot_Phase(ax, xx,yy, Phase, sw)
    # 2. Rho
    ax = axes[0][1]
    plot_prop(ax, xx,yy, Rho,prop_name="Density", prop_unit="kg/m$^\mathregular{3}$")
    # 3. H
    ax = axes[0][2]
    plot_prop(ax, xx,yy, H/1E6,prop_name="Specific enthalpy", prop_unit="MJ/kg",cmap='GnBu')
    # 3. Cp
    ax = axes[0][3]
    plot_prop(ax, xx,yy, Cp,prop_name="Isobaric specific heat", prop_unit="J/kg/K",cmap='YlOrRd')

    # 1. Mu_l
    ax=axes[1][0]
    plot_prop(ax, xx,yy, np.log10(Mu_l),prop_name="log$_{10}$ $\mu_l$", prop_unit="$\mu Pa \cdot s$",cmap='GnBu')
    # ax.grid()
    # 2. Mu_v
    ax=axes[1][1]
    plot_prop(ax, xx,yy, np.log10(Mu_v),prop_name="log$_{10}$ $\mu_v$", prop_unit="$\mu Pa \cdot s$",cmap='BuGn')

    # kappa: Isothermal compressibility
    ax=axes[1][2]
    kappa = np.array(state.IsothermalCompressibility).reshape(TT.shape)
    plot_prop(ax, xx,yy, np.log10(kappa), prop_name="log$_{10}$ $\kappa$", prop_unit="Pa$^{-1}$",cmap='BuGn')

    # beta: isobaric expansivity
    ax=axes[1][3]
    beta = np.array(state.IsobaricExpansivity).reshape(TT.shape)
    plot_prop(ax, xx,yy, np.log10((beta)), prop_name="log$_{10}$ $\\beta$", prop_unit="T$^{-1}$",cmap='BuGn')
    print(beta.min(),beta.max())
    for ax in axes.reshape(-1,):
        ax.axvline(constT - 273.15, ls='dashed', color='k')
        ax.axhline(constP/1E5, ls='dashdot', color='k')

    # set axis
    for ax in axes[1,:]: ax.set_xlabel('Temperature ($^{\circ}$C)')
    for ax in axes[:,0]: ax.set_ylabel('Pressure (bar)')

    # ax.set_xscale('log')
    savefig('Props_X%.0fwt'%(X0*100))

plot_props_X0(constX, sw_84)

# %%
# Profile along T-axis
# ---------------------------
def plot_Phase_Profile(ax, sw,x, Phase, showPhaseName=False):
    colors_phase={'Liquid':'lightblue', 'V+L': 'lightgreen', 'Vapor':'purple','Unknown':'lightgray', 'L+H':'violet','V+H':'orange'}
    start,end = 0, 0
    for i in range(0,len(Phase)):
        if(Phase[i]!=Phase[start]):
            end = i-1
            phasename = sw.phase_name(int(Phase[start]))
            ax.axvspan(x[start], x[end], color=colors_phase[sw.phase_name(int(Phase[start]))], label=phasename, alpha=0.4)
            if(showPhaseName): ax.text(x[start:end+1].mean(), np.array(ax.get_ylim()).mean(), phasename, va='center', ha='center')
            start = i
    # the last setment
    phasename = sw.phase_name(int(Phase[start]))
    ax.axvspan(x[end], x[-1], color=colors_phase[phasename])
    if(showPhaseName): ax.text(x[end:].mean(), np.array(ax.get_ylim()).mean(), phasename)
    ax.legend(ncol=4, loc='lower left', bbox_to_anchor=[0, 1])
def plot_profile_T(sw, P0=constP, X0=constX):
    # sw.showProgressBar(False)
    T_K = np.linspace(0.1, 999, 400) + 273.15
    P = T_K*0 + P0
    X = T_K*0 + X0
    state = sw.UpdateState_TPX(T_K, P, X)

    # print(state)
    x=T_K- 273.15
    # plot
    fig, axes = plt.subplots(2,4,figsize=(28,12), sharex=True,gridspec_kw={'wspace':0.3,'hspace':0.2})
    for ax in axes[-1,:]:
        ax.set_xlim(x.min(), x.max())
        ax.xaxis.set_major_locator(MultipleLocator(200))
        ax.xaxis.set_minor_locator(MultipleLocator(40))
        ax.set_xlabel('Temperature ($^{\circ}$C)')
    # Phase
    ax=axes[0][0]
    Phase = np.array(state.phase)
    plot_Phase_Profile(ax, sw, x, Phase, True)

    # Density
    ax=axes[0][1]
    ax.plot(x, np.array(state.Rho))
    plot_Phase_Profile(ax, sw, x, Phase)
    ax.set_ylabel('Density $\\rho$ (kg/m$^3$)')

    # Enthalpy
    ax=axes[0][2]
    ax.plot(x, np.array(state.H)/1E6)
    plot_Phase_Profile(ax, sw, x, Phase)
    ax.set_ylabel('Enthalpy H (MJ/kg)')

    # Cp
    ax = axes[0][3]
    ax.plot(x, np.array(state.Cp)/1E6, label='Analytical approach')
    plot_Phase_Profile(ax, sw, x, Phase)
    ax.set_ylabel('Isobaric specific heat $C_p$ (MJ/kg/K)')
    # test finite difference method result
    # dT = 0.1
    # state_T2 = sw.UpdateState_TPX(T_K+dT, P, X)
    # axes[0][3].plot(T, (np.array(state_T2.H) - np.array(state.H))/dT/1E6, label='Diff approach: dT=%.3f $^{\circ}$C'%(dT))
    # axes[0][3].legend()

    # Mu
    ax=axes[1][0]
    ax.set_yscale('log')
    ax.plot(x, np.array(state.Mu_l))
    plot_Phase_Profile(ax, sw, x, Phase)
    ax.set_ylabel('Dynamic viscosity $\mu_l$ (Pa s)')

    ax = axes[1][1]
    ax.set_yscale('log')
    ax.plot(x, np.array(state.Mu_v))
    plot_Phase_Profile(ax, sw, x, Phase)
    ax.set_ylabel('Dynamic viscosity $\mu_l$ (Pa s)')

    # Isothermal compressibility
    # print(state.dRhodT)
    ax=axes[1][2]
    ax.set_yscale('log')
    l,=ax.plot(x, -1.0/np.array(state.Rho_l)*np.array(state.dRhodT_l), label='Liquid', lw=2)
    v,=ax.plot(x, -1.0/np.array(state.Rho_v)*np.array(state.dRhodT_v), label='Vapor', lw=2)
    b,=ax.plot(x, -1.0/np.array(state.Rho)*np.array(state.dRhodT), label='Bulk', lw=1)
    # axes[1][2].plot(T, np.array(state.IsobaricExpansivity), label='Analytical approach')
    # axes[1][2].plot(T, -1.0/np.array(state.Rho_l)*(np.array(state_T2.Rho_l) - np.array(state.Rho_l))/dT, color=l.get_color(), ls='dashed', lw=2, label='Diff approach: Liquid: dT=%.3f $^{\circ}$C'%(dT))
    # axes[1][2].plot(T, -1.0/np.array(state.Rho_v)*(np.array(state_T2.Rho_v) - np.array(state.Rho_v))/dT, color=v.get_color(), ls='dashed', lw=2, label='Diff approach: Liquid: dT=%.3f $^{\circ}$C'%(dT))
    # axes[1][2].plot(T, -1.0/np.array(state.Rho)*(np.array(state_T2.Rho) - np.array(state.Rho))/dT, color=b.get_color(), ls='dashed', lw=2, label='Diff approach: Liquid: dT=%.3f $^{\circ}$C'%(dT))
    # ax.legend()
    plot_Phase_Profile(ax, sw, x, Phase)
    ax.set_ylabel('Isothermal compressibility $\kappa$ (1/Pa)')

    # Isobaric expansivity
    # print(state.dRhodT)
    ax=axes[1][3]
    ax.set_yscale('log')
    ax.plot(x, 1.0/np.array(state.Rho_l)*np.array(state.dRhodP_l), label='Liquid')
    ax.plot(x, 1.0/np.array(state.Rho_v)*np.array(state.dRhodP_v), label='Vapor')
    ax.plot(x, 1.0/np.array(state.Rho)*np.array(state.dRhodP), label='Bulk')
    # axes[1][2].plot(T, np.array(state.IsobaricExpansivity), label='Analytical approach')
    # axes[1][2].plot(T, -1.0/np.array(state.Rho)*(np.array(state_T2.Rho) - np.array(state.Rho))/dT, label='Diff approach: dT=%.3f $^{\circ}$C'%(dT))
    # ax.legend()
    plot_Phase_Profile(ax, sw, x, Phase)
    ax.set_ylabel('Isobaric expansivity $\\beta$ (1/K)')

    savefig('Props_Tprofile_P%.0fbar_X%.0fwt'%(P0/1E5, X0*100))
plot_profile_T(sw_84)


# %%
# Profile along P-axis
# ---------------------------
def plot_profile_P(sw, T0=constT, X0=constX):
    # sw.showProgressBar(False)
    P = np.linspace(0.1, 1000, 400) *1E5
    T_K = P*0 + T0
    X = P*0 + X0
    state = sw.UpdateState_TPX(T_K, P, X)
    # print(state)
    x = P/1E5
    # plot
    fig, axes = plt.subplots(2,4,figsize=(28,12), sharex=True,gridspec_kw={'wspace':0.3,'hspace':0.2})
    for ax in axes[-1,:]:
        ax.set_xlim(x.min(), x.max())
        ax.xaxis.set_major_locator(MultipleLocator(200))
        ax.xaxis.set_minor_locator(MultipleLocator(40))
        ax.set_xlabel('Pressure (bar)')

    # Phase
    ax=axes[0][0]
    Phase = np.array(state.phase)
    plot_Phase_Profile(ax, sw, x, Phase, True)

    # Density
    ax=axes[0][1]
    ax.plot(x, np.array(state.Rho))
    plot_Phase_Profile(ax, sw, x, Phase)
    ax.set_ylabel('Density $\\rho$ (kg/m$^3$)')

    # Enthalpy
    ax=axes[0][2]
    ax.plot(x, np.array(state.H)/1E6)
    plot_Phase_Profile(ax, sw, x, Phase)
    ax.set_ylabel('Enthalpy H (MJ/kg)')

    # Cp
    ax = axes[0][3]
    ax.plot(x, np.array(state.Cp)/1E6, label='Analytical approach')
    plot_Phase_Profile(ax, sw, x, Phase)
    ax.set_ylabel('Isobaric specific heat $C_p$ (MJ/kg/K)')
    # test finite difference method result
    # dT = 0.1
    # state_T2 = sw.UpdateState_TPX(T_K+dT, P, X)
    # axes[0][3].plot(T, (np.array(state_T2.H) - np.array(state.H))/dT/1E6, label='Diff approach: dT=%.3f $^{\circ}$C'%(dT))
    # axes[0][3].legend()

    # Mu
    ax=axes[1][0]
    ax.set_yscale('log')
    ax.plot(x, np.array(state.Mu_l))
    plot_Phase_Profile(ax, sw, x, Phase)
    ax.set_ylabel('Dynamic viscosity $\mu_l$ (Pa s)')

    ax = axes[1][1]
    ax.set_yscale('log')
    ax.plot(x, np.array(state.Mu_v))
    plot_Phase_Profile(ax, sw, x, Phase)
    ax.set_ylabel('Dynamic viscosity $\mu_l$ (Pa s)')

    # Isothermal compressibility
    # print(state.dRhodT)
    ax=axes[1][2]
    ax.set_yscale('log')
    l,=ax.plot(x, -1.0/np.array(state.Rho_l)*np.array(state.dRhodT_l), label='Liquid', lw=2)
    v,=ax.plot(x, -1.0/np.array(state.Rho_v)*np.array(state.dRhodT_v), label='Vapor', lw=2)
    b,=ax.plot(x, -1.0/np.array(state.Rho)*np.array(state.dRhodT), label='Bulk', lw=1)
    # axes[1][2].plot(T, np.array(state.IsobaricExpansivity), label='Analytical approach')
    # axes[1][2].plot(T, -1.0/np.array(state.Rho_l)*(np.array(state_T2.Rho_l) - np.array(state.Rho_l))/dT, color=l.get_color(), ls='dashed', lw=2, label='Diff approach: Liquid: dT=%.3f $^{\circ}$C'%(dT))
    # axes[1][2].plot(T, -1.0/np.array(state.Rho_v)*(np.array(state_T2.Rho_v) - np.array(state.Rho_v))/dT, color=v.get_color(), ls='dashed', lw=2, label='Diff approach: Liquid: dT=%.3f $^{\circ}$C'%(dT))
    # axes[1][2].plot(T, -1.0/np.array(state.Rho)*(np.array(state_T2.Rho) - np.array(state.Rho))/dT, color=b.get_color(), ls='dashed', lw=2, label='Diff approach: Liquid: dT=%.3f $^{\circ}$C'%(dT))
    # ax.legend()
    plot_Phase_Profile(ax, sw, x, Phase)
    ax.set_ylabel('Isothermal compressibility $\kappa$ (1/Pa)')

    # Isobaric expansivity
    # print(state.dRhodT)
    ax=axes[1][3]
    ax.set_yscale('log')
    ax.plot(x, 1.0/np.array(state.Rho_l)*np.array(state.dRhodP_l), label='Liquid')
    ax.plot(x, 1.0/np.array(state.Rho_v)*np.array(state.dRhodP_v), label='Vapor')
    ax.plot(x, 1.0/np.array(state.Rho)*np.array(state.dRhodP), label='Bulk')
    # axes[1][2].plot(T, np.array(state.IsobaricExpansivity), label='Analytical approach')
    # axes[1][2].plot(T, -1.0/np.array(state.Rho)*(np.array(state_T2.Rho) - np.array(state.Rho))/dT, label='Diff approach: dT=%.3f $^{\circ}$C'%(dT))
    # ax.legend()
    plot_Phase_Profile(ax, sw, x, Phase)
    ax.set_ylabel('Isobaric expansivity $\\beta$ (1/K)')

    savefig('Props_Pprofile_T%.0fdegC%.0fwt'%(T0 - 273.15, X0*100))

plot_profile_P(sw_84)

# %%
# Profile along X-axis
# ---------------------------
def plot_profile_X(sw, T0=constT, P0=constP):
    # sw.showProgressBar(False)
    X = np.linspace(1E-4, 1, 400)
    T_K = X*0 + T0
    P = X*0 + P0
    state = sw.UpdateState_TPX(T_K, P, X)
    # print(state)
    x = X*100
    # plot
    fig, axes = plt.subplots(2,4,figsize=(28,12), sharex=True,gridspec_kw={'wspace':0.3,'hspace':0.2})
    for ax in axes[-1,:]:
        ax.set_xlim(x.min(), x.max())
        ax.xaxis.set_major_locator(MultipleLocator(10))
        ax.xaxis.set_minor_locator(MultipleLocator(2))
        ax.set_xlabel('Bulk salinity (wt.% NaCl)')

    # Phase
    ax=axes[0][0]
    Phase = np.array(state.phase)
    plot_Phase_Profile(ax, sw, x, Phase, True)

    # Density
    ax=axes[0][1]
    ax.plot(x, np.array(state.Rho))
    plot_Phase_Profile(ax, sw, x, Phase)
    ax.set_ylabel('Density $\\rho$ (kg/m$^3$)')

    # Enthalpy
    ax=axes[0][2]
    ax.plot(x, np.array(state.H)/1E6)
    plot_Phase_Profile(ax, sw, x, Phase)
    ax.set_ylabel('Enthalpy H (MJ/kg)')

    # Cp
    ax = axes[0][3]
    ax.plot(x, np.array(state.Cp)/1E6, label='Analytical approach')
    plot_Phase_Profile(ax, sw, x, Phase)
    ax.set_ylabel('Isobaric specific heat $C_p$ (MJ/kg/K)')
    # test finite difference method result
    # dT = 0.1
    # state_T2 = sw.UpdateState_TPX(T_K+dT, P, X)
    # axes[0][3].plot(T, (np.array(state_T2.H) - np.array(state.H))/dT/1E6, label='Diff approach: dT=%.3f $^{\circ}$C'%(dT))
    # axes[0][3].legend()

    # Mu
    ax=axes[1][0]
    ax.set_yscale('log')
    ax.plot(x, np.array(state.Mu_l))
    plot_Phase_Profile(ax, sw, x, Phase)
    ax.set_ylabel('Dynamic viscosity $\mu_l$ (Pa s)')

    ax = axes[1][1]
    ax.set_yscale('log')
    ax.plot(x, np.array(state.Mu_v))
    plot_Phase_Profile(ax, sw, x, Phase)
    ax.set_ylabel('Dynamic viscosity $\mu_l$ (Pa s)')

    # Isothermal compressibility
    # print(state.dRhodT)
    ax=axes[1][2]
    ax.set_yscale('log')
    l,=ax.plot(x, -1.0/np.array(state.Rho_l)*np.array(state.dRhodT_l), label='Liquid', lw=2)
    v,=ax.plot(x, -1.0/np.array(state.Rho_v)*np.array(state.dRhodT_v), label='Vapor', lw=2)
    b,=ax.plot(x, -1.0/np.array(state.Rho)*np.array(state.dRhodT), label='Bulk', lw=1)
    # axes[1][2].plot(T, np.array(state.IsobaricExpansivity), label='Analytical approach')
    # axes[1][2].plot(T, -1.0/np.array(state.Rho_l)*(np.array(state_T2.Rho_l) - np.array(state.Rho_l))/dT, color=l.get_color(), ls='dashed', lw=2, label='Diff approach: Liquid: dT=%.3f $^{\circ}$C'%(dT))
    # axes[1][2].plot(T, -1.0/np.array(state.Rho_v)*(np.array(state_T2.Rho_v) - np.array(state.Rho_v))/dT, color=v.get_color(), ls='dashed', lw=2, label='Diff approach: Liquid: dT=%.3f $^{\circ}$C'%(dT))
    # axes[1][2].plot(T, -1.0/np.array(state.Rho)*(np.array(state_T2.Rho) - np.array(state.Rho))/dT, color=b.get_color(), ls='dashed', lw=2, label='Diff approach: Liquid: dT=%.3f $^{\circ}$C'%(dT))
    # ax.legend()
    plot_Phase_Profile(ax, sw, x, Phase)
    ax.set_ylabel('Isothermal compressibility $\kappa$ (1/Pa)')

    # Isobaric expansivity
    # print(state.dRhodT)
    ax=axes[1][3]
    ax.set_yscale('log')
    ax.plot(x, 1.0/np.array(state.Rho_l)*np.array(state.dRhodP_l), label='Liquid')
    ax.plot(x, 1.0/np.array(state.Rho_v)*np.array(state.dRhodP_v), label='Vapor')
    ax.plot(x, 1.0/np.array(state.Rho)*np.array(state.dRhodP), label='Bulk')
    # axes[1][2].plot(T, np.array(state.IsobaricExpansivity), label='Analytical approach')
    # axes[1][2].plot(T, -1.0/np.array(state.Rho)*(np.array(state_T2.Rho) - np.array(state.Rho))/dT, label='Diff approach: dT=%.3f $^{\circ}$C'%(dT))
    # ax.legend()
    plot_Phase_Profile(ax, sw, x, Phase)
    ax.set_ylabel('Isobaric expansivity $\\beta$ (1/K)')

    savefig('Props_Xprofile_T%.0fdegCP%.0fbar'%(T0 - 273.15, P0/1E5))

plot_profile_X(sw_84)
