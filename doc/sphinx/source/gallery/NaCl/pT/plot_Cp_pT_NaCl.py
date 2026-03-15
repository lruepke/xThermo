# -*- coding: utf-8 -*-
"""
7. Isobaric specific heat capacity
======================================
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
from xThermal import NaCl
salt = NaCl.cNaCl('IAPS84')

# Calculate and plot result
def plot_prop(salt, ax=None,name_prop='Cp',unit_prop='J/kg/K',T=[],p=[],cmap='GnBu'):
    if(len(T)==0): T = np.linspace(274,1273,400)
    if(len(p)==0): p = np.linspace(1E5,1000e5,400)
    TT,PP = np.meshgrid(T,p)
    prop = np.zeros_like(TT)
    for i in range(0,TT.shape[0]):
        for j in range(0,TT.shape[1]):
            props = salt.UpdateState_TPX(TT[i][j], PP[i][j])
            prop[i][j] = props.Cp
    # plot
    isSaveFig=False
    if(ax is None): isSaveFig = True
    if(ax==None):
        fig=plt.figure(figsize=(7,7))
        ax=plt.gca()
    CS = ax.contourf(TT-273.15,PP/1E5, prop, levels=50, cmap=cmap)
    ax_cb = ax.inset_axes([1.01, 0, 0.05, 1])
    plt.colorbar(CS, cax=ax_cb, orientation='vertical',label='%s (%s)'%(name_prop,unit_prop))
    # plot melting curve
    T_melting = np.array(salt.Melting_T(p))
    ax.plot(T_melting-273.15, p/1E5, color='orange')
    # labels
    ax.text(0.98,0.98,salt.name(),ha='right',va='top',bbox={'fc':'w','ec':'gray'}, transform=ax.transAxes)
    ax.set_xlabel('Temperature ($^{\circ}$C)')
    ax.set_ylabel('Pressure (bar)')
    if(isSaveFig): savefig('NaCl_%s_%s'%(name_prop, salt.name()))
    return TT-273.15,PP/1E5,prop

# %%
# IAPS84 EOS
# -------------------------
TT,PP,prop=plot_prop(salt)
