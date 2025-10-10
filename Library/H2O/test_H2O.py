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
# from iapws import IAPWS95_CoolProp 
# from iapws import IAPWS97
# import iapws.iapws95 as iapws

from xThermal import H2O
iaps84 = H2O.cIAPS84()
iapws95 = H2O.cIAPWS95_CoolProp()

fmt_figs=['svg','pdf']
figpath='.'

def savefig(figname):
    for fmt_fig in fmt_figs:
        figname_full = '%s/%s.%s'%(figpath,figname,fmt_fig)
        plt.savefig(figname_full, bbox_inches='tight')
        print('figure saved: ',figname_full)

def test_constants(water):
    print("\nTest %s"%(water.name()))
    print('Tmin [K]: %f '%(water.Tmin()))
    print('Tmax [K]: %f '%(water.Tmax()))
    print('pmin [Pa]: %f '%(water.pmin()))
    print('pmax [Pa]: %f '%(water.pmax()))
    print('T_triple [K]: %f '%(water.Ttriple()))
    print('T_crit [K]: %f '%(water.T_critical()))
    print('p_crit [Pa]: %f '%(water.p_critical()))
    print('rho_crit [kg/m^3]: %f '%(water.rhomass_critical()))
    # print('rho_crit [mol/m^3]: %f '%(water.rhomolar_critical()))
    print('Molar mass [kg/mol]: %f '%(water.molar_mass()))

def test_prop(water,name_prop='rho',T=None,p=None):
    if(T==None): T = np.linspace(274,1273,100)
    if(p==None): p = np.linspace(1E5,1000e5,100)
    TT,PP = np.meshgrid(T,p)
    prop = np.zeros_like(TT)
    for i in range(0,TT.shape[0]):
        for j in range(0,TT.shape[1]):
            water.UpdateState_TP(TT[i][j], PP[i][j])
            prop[i][j] = water.rhomass()
    # plot
    fig=plt.figure(figsize=(7,7))
    ax=plt.gca()
    CS = ax.contourf(TT-273.15,PP/1E5, prop, levels=50, cmap='GnBu')
    ax_cb = ax.inset_axes([1.01, 0, 0.05, 1])
    plt.colorbar(CS, cax=ax_cb, orientation='vertical',label=name_prop)
    # labels
    ax.text(0.98,0.98,water.name(),ha='right',va='top',bbox={'fc':'w','ec':'gray'}, transform=ax.transAxes)
    ax.set_xlabel('Temperature ($^{\circ}$C)')
    ax.set_ylabel('Pressure (bar)')

    savefig('H2O_%s_%s'%(name_prop, water.name()))
# 1.
# test_constants(iaps84)
# test_constants(iapws95)
# 2. 
for water in [iaps84, iapws95]: test_prop(water)
