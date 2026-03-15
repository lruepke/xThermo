# Heat a fluid with initial composition(bulk salinity) by increasing bulk enthalpy,
# visualize the phase change path
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
from matplotlib.ticker import MultipleLocator
import matplotlib.patches as mpatches
from matplotlib.collections import PatchCollection
mpl.rcParams['font.family'] = 'Arial'  # default font family
mpl.rcParams['mathtext.fontset'] = 'cm'  # font for math
from xThermal import H2ONaCl
sw = H2ONaCl.cH2ONaCl("IAPS84")

# make log-linear hybrid axes for isobaric slice
def makeLogLinearAxes(ax,width_logaxis=0.4,xlim_log=(2E-3,1),xmax=100,xloc_major=20,xloc_minor=4,yloc='left',color_log=(127/255,0,1),xlabel='Bulk salinity (wt % NaCl)',bkcolor='None',ylim=None, yloc_major=None,yloc_minor=None,ylabel=None):
    # log segment
    ax_log=ax.inset_axes([0,0,width_logaxis,1],facecolor=bkcolor)
    ax_log.set_xscale('log')
    ax_log.set_xlim(xlim_log)
    #     ax_log.xaxis.set_minor_locator(mpl.ticker.LogLocator(base=10.0, subs=(0.1,0.2,0.3,0.4, 0.5, 0.6, 0.7, 0.8, 0.9),numticks=10))
    ax_log.spines['right'].set_visible(False)
    ax_log.spines['bottom'].set_color(color_log)
    ax_log.spines['top'].set_color(color_log)
    ax_log.tick_params(axis='x', which='both', colors=color_log)
    ax_log.annotate("", xy=(0.03,-0.11), xycoords='axes fraction', xytext=(1,-0.11), textcoords='axes fraction', arrowprops=dict(arrowstyle="<->",connectionstyle="arc3",color=color_log),clip_on=False,zorder=0)
    ax_log.set_xlabel('Log-scale',color=color_log,bbox={'fc':'w','ec':'None'})

    # linear segment
    ax_linear=ax.inset_axes([width_logaxis,0, 1-width_logaxis, 1],facecolor=bkcolor)
    ax_linear.yaxis.set_label_position('right')
    ax_linear.yaxis.set_ticks_position('right')
    ax_linear.set_xlim(xlim_log[1], xmax+1)
    ax_linear.xaxis.set_major_locator(MultipleLocator(xloc_major))
    ax_linear.xaxis.set_minor_locator(MultipleLocator(xloc_minor))
    #ax_linear.spines['left'].set_visible(False)
    ax_linear.spines['left'].set_linewidth(0.5)
    ax_linear.spines['left'].set_linestyle('solid')
    ax_linear.spines['left'].set_color(color_log)
    ax_linear.spines['left'].set_capstyle("butt")
    ax_linear.annotate("", xy=(-0.03,-0.11), xycoords='axes fraction', xytext=(1,-0.11), textcoords='axes fraction', arrowprops=dict(arrowstyle="<->",connectionstyle="arc3"),clip_on=False,zorder=0)
    ax_linear.set_xlabel('Linear-scale',labelpad=10,bbox={'fc':'w','ec':'None'})
    # ytick location
    ax_linear.yaxis.set_ticklabels([]) if(yloc=='left') else ax_log.yaxis.set_ticklabels([])
    # share y axis
    # ax_linear.get_shared_y_axes().join(ax_log, ax_linear)
    for axis in [ax_log,ax_linear]:
        if(ylim!=None): axis.set_ylim(ylim)
        if(yloc_major!=None): axis.yaxis.set_major_locator(MultipleLocator(yloc_major))
        if(yloc_minor!=None): axis.yaxis.set_minor_locator(MultipleLocator(yloc_minor))
    if(ylabel!=None): ax_log.set_ylabel(ylabel) if(yloc=='left') else ax_linear.set_ylabel(ylabel)
    # hide base axes
    for spine in ax.spines: ax.spines[spine].set_visible(False)
    ax.xaxis.set_ticks([])
    ax.yaxis.set_ticks([])
    if(xlabel!=None): ax.set_xlabel(xlabel,labelpad=40)

    return ax_log,ax_linear

def plot_PhaseRegions(PhaseRegions, ax_linear_TPX, ax_log_TPX, ax_linear_HPX, ax_log_HPX, lw):
    phaseRegion_color = {}
    # plot regions
    for key in PhaseRegions.regions.keys():
        for i in range(0, len(PhaseRegions.regions[key])):
            for ax in [ax_linear_TPX, ax_log_TPX]:
                ax.fill(np.array(PhaseRegions.regions[key][i].X)*100, np.array(PhaseRegions.regions[key][i].T)-273.15, fc=PhaseRegions.regions[key][i].fc, ec='None')
            for ax in [ax_linear_HPX, ax_log_HPX]:
                ax.fill(np.array(PhaseRegions.regions[key][i].X)*100, np.array(PhaseRegions.regions[key][i].H)/1E6, fc=PhaseRegions.regions[key][i].fc, ec='None')
                phaseRegion_color[key] = PhaseRegions.regions[key][i].fc
            if(i==0): ax_linear_HPX.axhspan(-1000, -1000, fc=PhaseRegions.regions[key][i].fc, ec='k', linewidth=0.5, label=key.replace(': ', '\n'))
    # plot lines
    for key in PhaseRegions.lines.keys():
        for i in range(0, len(PhaseRegions.lines[key])):
            label = None #key if(i==0) else None
            for ax in [ax_linear_TPX, ax_log_TPX]:
                ax.plot(np.array(PhaseRegions.lines[key][i].X)*100, np.array(PhaseRegions.lines[key][i].T)-273.15, color=PhaseRegions.lines[key][i].color, label=label, lw=lw)
            for ax in [ax_linear_HPX, ax_log_HPX]:
                ax.plot(np.array(PhaseRegions.lines[key][i].X)*100, np.array(PhaseRegions.lines[key][i].H)/1E6, color=PhaseRegions.lines[key][i].color, label=label, lw=lw)
    # plot points
    for key in PhaseRegions.points.keys():
        for i in range(0, len(PhaseRegions.points[key])):
            label = None #key if(i==0) else None
            for ax in [ax_linear_TPX, ax_log_TPX]:
                X = np.array(PhaseRegions.points[key][i].X)*100
                if((ax == ax_linear_TPX) & (X.max() < 1)): continue
                ax.plot(X, np.array(PhaseRegions.points[key][i].T)-273.15, 'o', mfc=PhaseRegions.points[key][i].mfc, mec=PhaseRegions.points[key][i].mec, label=label, mew=lw)
            for ax in [ax_linear_HPX, ax_log_HPX]:
                X = np.array(PhaseRegions.points[key][i].X)*100
                if((ax == ax_linear_HPX) & (X.max() < 1)): continue
                ax.plot(X, np.array(PhaseRegions.points[key][i].H)/1E6, 'o', mfc=PhaseRegions.points[key][i].mfc, mec=PhaseRegions.points[key][i].mec, label=label, mew=lw)
    # legend
    ax_linear_HPX.axhspan(-1000, -1000, fc='w', ec='k', label='L', linewidth=0.5)
    ax_linear_HPX.legend(ncol=2, loc='upper right', frameon=True, columnspacing=1, framealpha=0.8, handlelength=1)
    phaseRegion_color['Liquid'] = 'w'
    return phaseRegion_color
def PhaseDiagramSlice_constP(P0=250E5,X0=0.4, Tmin=200+273.15, Tmax=H2ONaCl.T_MAX_VLH+50,fmt='jpeg',dpi=900,figpath='.',x0_log_linear=0.01,lw=0.5, points=30):
    # calculate phase regions
    PhaseRegions = sw.Slice_constP(P0)
    print(PhaseRegions)
    # calculate properties at const P and const X condition
    state = sw.UpdateState_TPX(Tmin, P0, X0)
    H_min = state.H
    state = sw.UpdateState_TPX(Tmax, P0, X0)
    H_max = state.H
    # calculate profile along T and H
    T, H = [], []
    if(P0<=H2ONaCl.P_Peak_VLH):
        TminTmax = sw.T_VLH_P0(P0)
        # ------------ refine H in VLH zone
        state1_Tmin = sw.UpdateState_TPX(TminTmax[0]-1, P0, X0)
        state2_Tmin = sw.UpdateState_TPX(TminTmax[0]+1, P0, X0)
        H = np.linspace(H_min, state1_Tmin.H, points)
        H = np.append(H,np.linspace(state1_Tmin.H, state2_Tmin.H, points))
        state1_Tmax = sw.UpdateState_TPX(TminTmax[1]-1, P0, X0)
        state2_Tmax = sw.UpdateState_TPX(TminTmax[1]+1, P0, X0)
        H = np.append(H,np.linspace(state2_Tmin.H, state1_Tmax.H, points))
        H = np.append(H,np.linspace(state1_Tmax.H, state2_Tmax.H, points))
        H = np.append(H,np.linspace(state2_Tmax.H, H_max, points))
    else:
        H = np.linspace(H_min, H_max, 200)
    P_, X_ = H*0 + P0, H*0 + X0
    H_ = H
    state_HPX = sw.UpdateState_HPX(H_, P_, X_)
    T_, H_l_, H_v_, H_h_ = np.array(state_HPX.T), np.array(state_HPX.H_l), np.array(state_HPX.H_v), np.array(state_HPX.H_h)
    phase_, X_l_, X_v_, S_l_, S_v_, S_h_ = np.array(state_HPX.phase), np.array(state_HPX.X_l), np.array(state_HPX.X_v), np.array(state_HPX.S_l), np.array(state_HPX.S_v), np.array(state_HPX.S_h)
    Rho_l_, Rho_v_, Rho_h_ = np.array(state_HPX.Rho_l), np.array(state_HPX.Rho_v), np.array(state_HPX.Rho_h)
    phase_unique = np.unique(phase_)
    # for pu in phase_unique: print(sw.phase_name(int(pu)))


    # plot
    for i in range(0, len(H_)): # len(H_)
        # i= 5
        # Plot phase diagram
        plt.close()
        fig,axes=plt.subplots(1,2,figsize=(16,5),gridspec_kw={'wspace':0.05})
        xlim_log=(2E-4,x0_log_linear*100) # wt.% NaCl
        # make log-linear hybrid axis
        ax_log_TPX,ax_linear_TPX=makeLogLinearAxes(axes[0],xlim_log=xlim_log,ylabel='Temperature ($^{\circ}$C)',ylim=(1,1000),yloc_major=200,yloc_minor=40)
        ax_log_HPX,ax_linear_HPX=makeLogLinearAxes(axes[1],xlim_log=xlim_log,yloc='right',ylabel='Specific enthalpy (MJ/kg)',ylim=(0,np.array(PhaseRegions.lines['V+L:V'][0].H).max()/1E6),yloc_major=1, yloc_minor=0.2)
        # ax_water_TPX=createWaterAxes(ax_log_TPX)
        # ax_water_HPX=createWaterAxes(ax_log_HPX)
        for ax_log in [ax_log_TPX, ax_log_HPX]:
            ax_log.xaxis.set_major_locator(mpl.ticker.LogLocator(base=10.0, subs=(1.0,),numticks=10))
            ax_log.xaxis.set_minor_locator(mpl.ticker.LogLocator(base=10.0, subs=(0.1,0.2,0.3,0.4, 0.5, 0.6, 0.7, 0.8, 0.9),numticks=10))
        # ----- plot profile along T and H
        ax_profile_TPX = axes[0].inset_axes([0, 1.05, 1, 0.3])
        ax_profile_HPX = axes[1].inset_axes([0, 1.05, 1, 0.3])
        for ax, label, ytickpos, xlim in zip([ax_profile_HPX, ax_profile_TPX], [ax_linear_HPX.get_ylabel(), ax_log_TPX.get_ylabel()], ['right', 'left'], [ax_linear_HPX.get_ylim(), ax_linear_TPX.get_ylim()]):
            ax.xaxis.set_label_position('top')
            ax.xaxis.set_ticks_position('top')
            ax.set_xlabel(label)
            ax.yaxis.set_ticks_position(ytickpos)
            ax.set_xlim(xlim)

        T, H, H_l, H_v, H_h, P, X, phase = T_[0:i+1], H_[0:i+1], H_l_[0:i+1], H_v_[0:i+1], H_h_[0:i+1], P_[0:i+1], X_[0:i+1], phase_[0:i+1]
        X_l, X_v, S_l, S_v, S_h, Rho_l, Rho_v, Rho_h = X_l_[0:i+1], X_v_[0:i+1], S_l_[0:i+1], S_v_[0:i+1], S_h_[0:i+1], Rho_l_[0:i+1], Rho_v_[0:i+1], Rho_h_[0:i+1]
        # plot regions
        phaseRegion_color = plot_PhaseRegions(PhaseRegions, ax_linear_TPX, ax_log_TPX, ax_linear_HPX, ax_log_HPX, lw)
        # plot profile on phase region
        ax_path_TPX = ax_log_TPX if(X0<x0_log_linear) else ax_linear_TPX
        ax_path_HPX = ax_log_HPX if(X0<x0_log_linear) else ax_linear_HPX
        ms = 4
        mfc_l, mfc_v, mfc_h = 'r', 'b', 'purple'
        # L+H
        ind = (phase == H2ONaCl.TwoPhase_LH)
        ax_path_TPX.plot(X_l[ind]*100, T[ind]-273.15, '.', ms = ms, mfc=mfc_l, mec='None')
        ax_path_HPX.plot(X_l[ind]*100, H_l[ind]/1E6, '.', ms = ms, mfc=mfc_l, mec='None')
        for ax1, ax2 in zip([ax_log_TPX, ax_linear_TPX], [ax_log_HPX, ax_linear_HPX]):
            ax1.plot(X_v[ind]*0 + 100, T[ind]-273.15, '.', ms = ms, mfc=mfc_h, mec='None')
            ax2.plot(X_v[ind]*0 + 100, H_h[ind]/1E6, '.', ms = ms, mfc=mfc_h, mec='None')
        # Liquid
        ind = (phase == H2ONaCl.SinglePhase_L)
        ax_path_TPX.plot(X_l[ind]*100, T[ind]-273.15, '.', ms = ms, mfc=mfc_l, mec='None')
        ax_path_HPX.plot(X_l[ind]*100, H[ind]/1E6, '.', ms = ms, mfc=mfc_l, mec='None')
        # V+L
        ind = (phase == H2ONaCl.TwoPhase_VL)
        ax_path_TPX.plot(X_l[ind]*100, T[ind]-273.15, '.', ms = ms, mfc=mfc_l, mec='None')
        ax_path_HPX.plot(X_l[ind]*100, H_l[ind]/1E6, '.', ms = ms, mfc=mfc_l, mec='None')
        for ax1, ax2 in zip([ax_log_TPX, ax_linear_TPX], [ax_log_HPX, ax_linear_HPX]):
            ax1.plot(X_v[ind]*100, T[ind]-273.15, '.', ms = ms, mfc=mfc_v, mec='None')
            ax2.plot(X_v[ind]*100, H_v[ind]/1E6, '.', ms = ms, mfc=mfc_v, mec='None')
        # V+H
        ind = (phase == H2ONaCl.TwoPhase_VH)
        for ax1, ax2 in zip([ax_log_TPX, ax_linear_TPX], [ax_log_HPX, ax_linear_HPX]):
            ax1.plot(X_v[ind]*100, T[ind]-273.15, '.', ms = ms, mfc=mfc_v, mec='None')
            ax2.plot(X_v[ind]*100, H_v[ind]/1E6, '.', ms = ms, mfc=mfc_v, mec='None')
        ax_linear_TPX.plot(X_v[ind]*0 + 100, T[ind]-273.15, '.', ms = ms, mfc=mfc_h, mec='None', zorder=20)
        ax_linear_HPX.plot(X_v[ind]*0 + 100, H_h[ind]/1E6, '.', ms = ms, mfc=mfc_h, mec='None', zorder=20)
        # plot bulk salinity path
        ax_path_TPX.plot(T*0 + X0*100, T-273.15, color='gray', linewidth=0.5, ls='dashed')
        ax_path_HPX.plot(T*0 + X0*100, H/1E6, color='gray', linewidth=0.5, ls='dashed')
        # plot bulk salinity with saturation and phase salinity
        rhom=S_l[-1]*Rho_l[-1] + S_v[-1]*Rho_v[-1] + S_h[-1]*Rho_h[-1]
        X_f = (S_l[-1]*Rho_l[-1]*X_l[-1] + S_v[-1]*Rho_v[-1]*X_v[-1])/rhom # mass fraction of salt in fluid
        X_s = (S_h[-1]*Rho_h[-1])/rhom # mass fraction of salt in solid
        Xm = X_f + X_s
        if(phase[-1] == H2ONaCl.SinglePhase_L):
            ax_path_TPX.plot(X0*100, T[-1] - 273.15, 'o', mfc=mfc_l, mec='k')
            ax_path_HPX.plot(X0*100, H[-1]/1E6, 'o', mfc=mfc_l, mec='k')
        else:
            r1,r2=1,0.4
            w_ax,h_ax=1.5, 0.105
            min_display,max_display,data_display=ax_path_HPX.transAxes.transform((0,0)),ax_path_HPX.transAxes.transform((1,1)),ax_path_HPX.transData.transform((X0*100,H[-1]/1E6))
            x0y0_frac=(data_display-min_display)/(max_display-min_display)
            ax_pie = ax_path_HPX.inset_axes([x0y0_frac[0]-w_ax/2, x0y0_frac[1]-h_ax/2, w_ax,h_ax], transform=ax_path_HPX.transAxes, zorder=10)
            ax_pie.pie([S_l[-1], S_v[-1], S_h[-1]],colors=[mfc_l,mfc_v,mfc_h], shadow=False,radius=r1, wedgeprops=dict(width=r2, edgecolor='None',lw=0.1, alpha=0.8),startangle=0)
            ax_pie.pie([X_f/Xm, X_s/Xm],colors=['c',mfc_h], shadow=False,radius=r1-r2, wedgeprops=dict(width=r1-r2-0.05, edgecolor='w',lw=0.1, alpha=0.8),startangle=0)
            # TPX space
            min_display,max_display,data_display=ax_path_TPX.transAxes.transform((0,0)),ax_path_TPX.transAxes.transform((1,1)),ax_path_TPX.transData.transform((X0*100,T[-1]-273.15))
            x0y0_frac=(data_display-min_display)/(max_display-min_display)
            ax_pie = ax_path_TPX.inset_axes([x0y0_frac[0]-w_ax/2, x0y0_frac[1]-h_ax/2, w_ax,h_ax], transform=ax_path_TPX.transAxes, zorder=10)
            ax_pie.pie([S_l[-1], S_v[-1], S_h[-1]],colors=[mfc_l,mfc_v,mfc_h], shadow=False,radius=r1, wedgeprops=dict(width=r2, edgecolor='None',lw=0.1, alpha=0.8),startangle=0)
            ax_pie.pie([X_f/Xm, X_s/Xm],colors=['c',mfc_h], shadow=False,radius=r1-r2, wedgeprops=dict(width=r1-r2-0.05, edgecolor='w',lw=0.1, alpha=0.8),startangle=0)
        # ----- plot profile along T and H
        ax_profile_TPX_xv = ax_profile_TPX.twinx()
        ax_profile_HPX_xv = ax_profile_HPX.twinx()
        ax_profile_TPX_xv.set_ylim(0, X_v_.max()*100)
        ax_profile_TPX.set_ylim(0, X_l_.max()*100)
        ax_profile_HPX.yaxis.set_ticks([])
        ax_profile_TPX_xv.spines["right"].set_position(("axes", 2.05))
        ax_profile_TPX_xv.spines['right'].set_color(mfc_v)
        ax_profile_TPX_xv.tick_params(axis='y', colors=mfc_v)
        ax_profile_TPX_xv.set_ylabel('$X_{v}$ (wt. % NaCl)', color=mfc_v)
        for spine in ['left', 'top', 'bottom']: ax_profile_TPX_xv.spines[spine].set_color('None')
        ax_profile_TPX.spines['left'].set_color(mfc_l)
        ax_profile_TPX.tick_params(axis='y', colors=mfc_l)
        ax_profile_TPX.set_ylabel('$X_{l}$ (wt. % NaCl)', color=mfc_l)
        ax_profile_TPX.xaxis.set_minor_locator(MultipleLocator(40))
        ax_profile_HPX.xaxis.set_minor_locator(MultipleLocator(0.2))
        # plot phase regions along profile
        ind_start, end_start = 0, 0
        for i_phase in range(0,len(phase)):
            phaseName = sw.phase_name(int(phase[ind_start]))
            if(phaseName=='V+L+H'): phaseName = 'V+L+H: T=%.3f$^{\circ}$C'%(T[ind_start]-273.15)
            fc = phaseRegion_color[phaseName]
            if(phase[ind_start]!=phase[i_phase]):
                ax_profile_TPX.axvspan(T[ind_start]-273.15, T[i_phase-1]-273.15, fc=fc, alpha=0.5)
                ax_profile_HPX.axvspan(H[ind_start]/1E6, H[i_phase-1]/1E6, fc=fc, alpha=0.5)
                ind_start = i_phase-1
            if(i_phase==(len(phase)-1)):
                ax_profile_TPX.axvspan(T[ind_start]-273.15, T[-1]-273.15, fc=fc, alpha=0.5)
                ax_profile_HPX.axvspan(H[ind_start]/1E6, H[-1]/1E6, fc=fc, alpha=0.5)
        ax_profile_TPX.plot(T-273.15, X_l*100, color=mfc_l)
        ax_profile_TPX_xv.plot(T-273.15, X_v*100, color=mfc_v)
        ax_profile_HPX.plot(H/1E6, X_l*100, color=mfc_l)
        ax_profile_HPX_xv.plot(H/1E6, X_v*100, color=mfc_v)
        ax_profile_HPX_xv.axis('off')
        ax_profile_HPX_xv.set_ylim(ax_profile_TPX_xv.get_ylim())
        ax_profile_HPX.set_ylim(ax_profile_TPX.get_ylim())
        # plot saturation
        ax_profile_HPX_Sh = ax_profile_HPX.twinx()
        ax_profile_HPX_Svl = ax_profile_HPX.twinx()
        for ax_s, c,label, ylim in zip([ax_profile_HPX_Sh, ax_profile_HPX_Svl], ['purple', 'c'], ['$S_h$ (%)', '$S_l, S_v$ (%)'], [(-0.05, S_h_.max()*100*1.1), (-5, 105)]):
            ax_s.spines['right'].set_color(c)
            ax_s.tick_params(axis='y', colors=c)
            ax_s.set_ylabel(label, color=c)
            for spine in ['left','top','bottom']: ax_s.spines[spine].set_color('None')
            ax_s.set_ylim(ylim)
        # ax_profile_HPX_S.set_ylim(-10,110)
        ax_profile_HPX_Sh.tick_params(axis='y', labelleft=True, labelright=False)
        ax_profile_HPX_Sh.plot(H/1E6, S_h*100, ls='dashed', color='purple', lw=1, label='$S_h$')
        ax_profile_HPX_Svl.plot(H/1E6, S_l*100, ls='solid', color='c', lw=1, label='$S_l$')
        ax_profile_HPX_Svl.plot(H/1E6, S_v*100, ls='dashed', color='c', lw=1, label='$S_v$')
        ax_profile_HPX_Sh.spines["right"].set_position(("axes", 0))
        ax_profile_HPX_Svl.spines["right"].set_position(("axes", 0.06))
        ax_profile_HPX_Svl.legend()

        # plot text
        ax_log_HPX.text(-0.02, -0.18, '(Guo et al., 2023, CG)', ha='center', transform=ax_log_HPX.transAxes, fontsize=12, fontweight='bold', color='blue')
        for ax, TorH, str_TorH, unit_TorH in zip([ax_log_TPX, ax_log_HPX], [T[-1]-273.15, H[-1]/1E6], ['T', 'H'], ['$^{\circ}$C', 'MJ/kg']):
            ax.text(0.02, 0.02, 'P = %.0f MPa'%(P0/1E6), ha='left', transform=ax.transAxes, fontsize=12)
            ax.text(0.02, 0.08, '%s = %.2f %s'%(str_TorH, TorH, unit_TorH), ha='left', transform=ax.transAxes, fontsize=12, fontweight='bold')
        # plot phase saturation
        w_rect,h_rect,x0,y0=0.8,0.06,0.02,0.14
        rect_l = mpatches.Rectangle((x0, y0),w_rect*S_l[-1],h_rect,fc=mfc_l, ec="None",transform=ax_log_HPX.transAxes)
        ax_log_HPX.add_patch(rect_l)
        rect_v = mpatches.Rectangle((x0 + w_rect*S_l[-1], y0),w_rect*S_v[-1],h_rect,fc=mfc_v, ec="None",transform=ax_log_HPX.transAxes)
        ax_log_HPX.add_patch(rect_v)
        rect_h = mpatches.Rectangle((x0 + w_rect*S_l[-1] + w_rect*S_v[-1], y0),w_rect*S_h[-1],h_rect,fc=mfc_h, ec="None",transform=ax_log_HPX.transAxes)
        ax_log_HPX.add_patch(rect_h)
        if(S_l[-1]!=0):
            ax_log_HPX.text(x0, y0+h_rect+0.01, 'Liquid', ha='left', va='bottom', color=mfc_l, transform=ax_log_HPX.transAxes)
            if(S_l[-1]>0.15): ax_log_HPX.text(x0 + w_rect*S_l[-1]/2.0, y0+h_rect/2.0, '%.0f %%'%(S_l[-1]*100), ha='center', va='center', color='w', transform=ax_log_HPX.transAxes, fontweight='bold')
        if(S_v[-1]!=0):
            ax_log_HPX.text(x0 + w_rect/2.0, y0+h_rect+0.01, 'Vapor', ha='center', va='bottom', color=mfc_v, transform=ax_log_HPX.transAxes)
            if(S_v[-1]>0.15):
                if(S_v[-1]>0.9):
                    ax_log_HPX.text(x0 + w_rect*S_l[-1] + w_rect*S_v[-1]/2.0, y0+h_rect/2.0, '%.3f %%'%(S_v[-1]*100), ha='center', va='center', color='w', transform=ax_log_HPX.transAxes, fontweight='bold')
                else:
                    ax_log_HPX.text(x0 + w_rect*S_l[-1] + w_rect*S_v[-1]/2.0, y0+h_rect/2.0, '%.0f %%'%(S_v[-1]*100), ha='center', va='center', color='w', transform=ax_log_HPX.transAxes, fontweight='bold')
        if(S_h[-1]!=0):
            ax_log_HPX.text(x0 + w_rect, y0+h_rect+0.01, 'Halite', ha='right', va='bottom', color=mfc_h, transform=ax_log_HPX.transAxes)
            if(S_h[-1]>0.15): ax_log_HPX.text(x0 + w_rect*S_l[-1] + w_rect*S_v[-1] + w_rect*S_h[-1]/2.0, y0+h_rect/2.0, '%.0f %%'%(S_h[-1]*100), ha='center', va='center', color='w', transform=ax_log_HPX.transAxes, fontweight='bold')
        ax_log_HPX.text(0.02, y0+h_rect+0.06, 'Phase saturation (%)', ha='left', va='bottom', color='gray', transform=ax_log_HPX.transAxes)
        # mass fraction of salt in fluid and halite: display how much salt is precipitated
        y0 = 0.34
        rect_f = mpatches.Rectangle((x0, y0),w_rect*X_f/Xm,h_rect,fc='c', ec="None",transform=ax_log_HPX.transAxes)
        ax_log_HPX.add_patch(rect_f)
        rect_s = mpatches.Rectangle((x0 + w_rect*X_f/Xm, y0),w_rect*X_s/Xm,h_rect,fc=mfc_h, ec="None",transform=ax_log_HPX.transAxes)
        ax_log_HPX.add_patch(rect_s)
        if(X_s!=0):
            ax_log_HPX.text(x0+w_rect, y0+h_rect+0.01, 'Solid', ha='right', va='bottom', color=mfc_h, transform=ax_log_HPX.transAxes)
            if((X_s/Xm)>0.15):
                if((X_s/Xm)>0.9):
                    ax_log_HPX.text(x0 + w_rect*X_f/Xm + w_rect*X_s/Xm/2.0, y0+h_rect/2.0, '%.3f %%'%(X_s/Xm*100), ha='center', va='center', color='w', transform=ax_log_HPX.transAxes, fontweight='bold')
                else:
                    ax_log_HPX.text(x0 + w_rect*X_f/Xm + w_rect*X_s/Xm/2.0, y0+h_rect/2.0, '%.0f %%'%(X_s/Xm*100), ha='center', va='center', color='w', transform=ax_log_HPX.transAxes, fontweight='bold')
        if(X_f!=0):
            ax_log_HPX.text(x0, y0+h_rect+0.01, 'Fluid', ha='left', va='bottom', color='c', transform=ax_log_HPX.transAxes)
            if((X_f/Xm)>0.15): ax_log_HPX.text(x0 + w_rect*X_f/Xm/2.0, y0+h_rect/2.0, '%.0f %%'%(X_f/Xm*100), ha='center', va='center', color='w', transform=ax_log_HPX.transAxes, fontweight='bold')
        ax_log_HPX.text(0.02, y0+h_rect+0.06, 'Mass fraction of salt (%)', ha='left', va='bottom', color='gray', transform=ax_log_HPX.transAxes)
        plt.savefig('%s/%04d.%s'%(figpath,i,fmt),bbox_inches='tight', dpi=dpi)
        print('%d/%d complete'%(i, len(H_)))
    plt.savefig('slice.pdf', bbox_inches = 'tight')
    print('Run command: ffmpeg -r 10 -f image2 -i %%4d.%s  -vcodec libx264 -vf "crop=trunc(iw/2)*2:trunc(ih/2)*2"  -pix_fmt yuv420p animation_X%.1f.mp4'%(fmt, X0*100))

PhaseDiagramSlice_constP()