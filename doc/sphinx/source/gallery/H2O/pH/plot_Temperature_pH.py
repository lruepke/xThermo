# -*- coding: utf-8 -*-
"""
2. Temperature
=========================
"""

import numpy as np 
import time
import linecache
import matplotlib.pyplot as plt
import matplotlib as mpl
from matplotlib import patches
import matplotlib.ticker as ticker
from matplotlib.ticker import MultipleLocator
mpl.rcParams['font.family'] = 'Arial'  #default font family
mpl.rcParams['mathtext.fontset'] = 'cm' #font for math
fmt_figs=['pdf'] #['svg','pdf']
figpath='.'
def savefig(figname):
    for fmt_fig in fmt_figs:
        figname_full = '%s/%s.%s'%(figpath,figname,fmt_fig)
        plt.savefig(figname_full, bbox_inches='tight')
        print('figure saved: ',figname_full)

# Import package of xThermal
from xThermal import H2O
iaps84 = H2O.cIAPS84()
# iapws95_CoolProp = H2O.cIAPWS95_CoolProp(). #coolprop not exposed to python
iapws95 = H2O.cIAPWS95()

# Calculate and plot result
def plot_prop(water, ax=None,name_prop='Temperature',unit_prop='K',H=[],p=[],cmap='GnBu'):
    if(len(H)==0): H = np.linspace(0.1E6, 4.5E6,150)
    if(len(p)==0): p = np.linspace(1E5,1000e5,150)
    HH,PP = np.meshgrid(H,p)
    prop = np.zeros_like(HH)
    for i in range(0,HH.shape[0]):
        for j in range(0,HH.shape[1]):
            props = water.UpdateState_HPX(HH[i][j], PP[i][j])
            prop[i][j] = props.T
    # plot
    isSaveFig=False
    if(ax is None): isSaveFig = True
    if(ax==None):
        fig=plt.figure(figsize=(7,7))
        ax=plt.gca()
    CS = ax.contourf(HH/1E6,PP/1E5, prop, levels=50, cmap=cmap)
    ax_cb = ax.inset_axes([1.01, 0, 0.05, 1])
    plt.colorbar(CS, cax=ax_cb, orientation='vertical',label='%s (%s)'%(name_prop,unit_prop))
    # labels
    ax.text(0.98,0.98,water.name(),ha='right',va='top',bbox={'fc':'w','ec':'gray'}, transform=ax.transAxes)
    ax.set_xlabel('Specific enthalpy (MJ/kg)')
    ax.set_ylabel('Pressure (bar)')
    if(isSaveFig): savefig('H2O_%s_%s'%(name_prop, water.name()))
    return HH/1E6, PP/1E5,prop

def plot_error(ax, XX,YY,ZZ, eosA, eosB,cmap='RdBu',unit='K'):
    CS=ax.contourf(XX,YY,ZZ,levels=50,cmap=cmap, norm=mpl.colors.CenteredNorm())
    ax_cb = ax.inset_axes([1.01, 0, 0.05, 1])
    plt.colorbar(CS, cax=ax_cb, orientation='vertical',label='Diff (%s)'%(unit))
    ax.text(0.98,0.98,'%s-%s'%(eosA.name(),eosB.name()),ha='right',va='top',bbox={'fc':'w','ec':'gray'}, transform=ax.transAxes)
    ax.set_xlabel('Specific enthalpy (MJ/kg)')

# %%
# IAPS84 EOS
# -------------------------
TT,PP,prop=plot_prop(iaps84)

# %%
# Comparison: IAPS84 and IAPWS95
# --------------------------------------------------

fig, axes = plt.subplots(1,5,figsize=(30,5),gridspec_kw={'wspace':0.35})
TT,PP,prop84 = plot_prop(iaps84,ax=axes[0])
TT,PP,prop95 = plot_prop(iapws95,ax=axes[1])
#TT,PP,prop95_CoolProp = plot_prop(iapws95_CoolProp,ax=axes[2])
plot_error(axes[3],TT,PP,prop95 - prop84,eosA=iapws95, eosB=iaps84)
#plot_error(axes[4],TT,PP,prop95 - prop95_CoolProp,eosA=iapws95, eosB=iapws95_CoolProp)
for ax in axes[1:]: ax.set_ylabel(None)
savefig('T_diff_84_95')
