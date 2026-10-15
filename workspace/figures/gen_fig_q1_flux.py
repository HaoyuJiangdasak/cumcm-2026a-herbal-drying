"""问题 1 表面通量与累积失水：面积图 + 驱动力对照。

a：热流密度 |q| 与水分通量 j 随时间变化，面积填充示意瞬时强度；
b：累积失水量与驱动温差 (T_a-T_R)、驱动浓度差 (C_R-C_a) 的同步演化。
全部取自 problem_1_results.json 的 flux_series / mass_series（1800 个 1 s 输出）。
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _figcommon import *  # noqa: F401,F403

d1 = load('problem_1_results')
fs = d1['flux_series']
ms = d1['mass_series']
dg = d1['diagnostics']
t = np.asarray(fs['times_s'], dtype=float) / 60.0
qh = np.abs(np.asarray(fs['q_heat_W_m2'], dtype=float))
qm = np.asarray(fs['q_mass_kg_m2s'], dtype=float) * 1e6
Ta = np.asarray(fs['Ta'], dtype=float)
Ts = np.asarray(fs['T_surface'], dtype=float)
Ca = np.asarray(fs['Ca'], dtype=float)
Cs = np.asarray(fs['C_surface'], dtype=float)
cum = np.asarray(ms['flux_cum'], dtype=float)

fig = plt.figure(figsize=(6.4, 4.3), layout='constrained')
gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 1.0], hspace=0.05)
ax = fig.add_subplot(gs[0, 0])
ax2 = ax.twinx()
axb = fig.add_subplot(gs[1, 0], sharex=ax)
axb2 = axb.twinx()

ax.fill_between(t, 0, qh, color=C_HEAT, alpha=0.20, lw=0, zorder=1)
lh, = ax.plot(t, qh, color=C_HEAT, lw=2.0, zorder=3, label='热流密度 $|q_h|$')
lm, = ax2.plot(t, qm, color=C_MOIST, lw=2.0, ls='--', zorder=3,
               label='水分通量 $j_m$')
ax.set_ylabel('$|q_h|$ (W/m$^2$)', color=C_HEAT)
ax2.set_ylabel('$j_m$ ($10^{-6}$ kg/(m$^2\\cdot$s))', color=C_MOIST)
ax.tick_params(axis='y', colors=C_HEAT)
ax2.tick_params(axis='y', colors=C_MOIST)
ax.set_xlim(0, t[-1])
ax.set_ylim(0, 152)
ax2.set_ylim(1.05, 2.30)
ax2.grid(False)
ax.tick_params(labelbottom=False)
ax.legend(handles=[lh, lm], loc='upper left', frameon=False,
          fontsize=FS_S, labelspacing=0.35, handlelength=1.6)
ax.annotate('%.2f' % qh[-1], xy=(t[-1], qh[-1]), xytext=(24.0, 132.0),
            fontsize=FS_S, color=C_TEXT, ha='center',
            arrowprops=dict(arrowstyle='-|>', color=C_LINE, lw=1.0,
                            shrinkA=2, shrinkB=3))
panel(ax, 'a')

axb.fill_between(t, 0, cum, color=C_MEAN, alpha=0.22, lw=0, zorder=1)
lc, = axb.plot(t, cum, color=C_MEAN, lw=2.0, zorder=3, label='累积失水 $M_{cum}$')
ld, = axb2.plot(t, Ta - Ts, color=C_HEAT, lw=1.7, ls='-.', zorder=3,
                label='驱动温差 $T_a-T_R$')
le, = axb2.plot(t, Cs - Ca, color=C_MOIST, lw=1.7, ls=':', zorder=3,
                label='驱动浓差 $C_R-C_a$')
axb.set_xlabel('时间 $t$ (min)')
axb.set_ylabel('累积失水 (kg/m$^2$)', color=C_MEAN_T)
axb2.set_ylabel('驱动温差 (℃) / 浓差 (kg/kg)')
axb.tick_params(axis='y', colors=C_MEAN_T)
axb.set_xlim(0, t[-1])
axb.set_ylim(0, 0.196)
axb2.set_ylim(0, 6.2)
axb2.grid(False)
axb.legend(handles=[lc, ld, le], loc='upper left', ncol=1, frameon=False,
           fontsize=FS_S, labelspacing=0.35, handlelength=1.6)
axb.text(20.6, 0.028, '终值 %.5f' % dg['flux_total'], fontsize=FS_S,
         color=ink(C_MEAN, on=C_MEAN, alpha=0.22), ha='center', va='center')
panel(axb, 'b')

out(fig, 'fig_q1_flux')
