"""问题 3 干燥终点判据：max_r C 逼近阈值的全过程 + 终点邻域放大。

a：max_r C(r,t)、表面 C(R,t) 与平均 C̄ 三条曲线随 57.5 h 全程下降，水平线为
阈值 0.15；b：终点附近 1 min 分辨率的括号区间，显示判据被单调穿越一次。
数据来源 problem_3_results.json 的 criterion_series / t_dry / result_grid。
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _figcommon import *  # noqa: F401,F403

d3 = load('problem_3_results')
cs = d3['criterion_series']
td = d3['t_dry']
th = np.asarray(cs['times_h'], dtype=float)
Cmax = np.asarray(cs['C_max'], dtype=float)
Csur = np.asarray(cs['C_surface'], dtype=float)
Cbar = np.asarray(cs['C_bar'], dtype=float)
t_dry_h = float(td['hours'])

fig = plt.figure(figsize=(6.4, 3.2), layout='constrained')
gs = fig.add_gridspec(1, 2, width_ratios=[1.42, 1.0], wspace=0.10)
ax = fig.add_subplot(gs[0, 0])
axz = fig.add_subplot(gs[0, 1])

ax.axhspan(0, C_TH, color=C_GRN1, alpha=0.45, lw=0, zorder=0)
ax.plot(th, Cmax, color=C_HEAT, lw=2.2, zorder=4, label='判据量 $C_{\\max}$')
ax.plot(th, Cbar, color=C_MEAN, lw=1.8, ls='-.', zorder=3, label='平均 $\\bar{C}$')
ax.plot(th, Csur, color=C_SURF, lw=1.8, ls='--', zorder=3, label='表面 $C(R)$')
ax.axhline(C_TH, color=C_LINE, lw=1.3, ls=':', zorder=2)
ax.axvline(t_dry_h, color=C_ACC2, lw=1.4, zorder=2)
ax.set_xlabel('时间 $t$ (h)')
ax.set_ylabel('干基含水率 $C$ (kg/kg)')
ax.set_xlim(0, th[-1])
ax.set_ylim(0, 2.72)
ax.legend(loc='upper right', frameon=False, fontsize=FS_S,
          labelspacing=0.35, handlelength=1.7)
ax.annotate('$t_{dry}$=%.4f h' % t_dry_h, xy=(t_dry_h, C_TH),
            xytext=(30.0, 0.74), fontsize=FS_S, color=C_TEXT, ha='center',
            arrowprops=dict(arrowstyle='-|>', color=C_LINE, lw=1.0,
                            shrinkA=2, shrinkB=3))
panel(ax, 'a')

w = 1.25
sel = th >= t_dry_h - w
tz, Cz = th[sel], Cmax[sel]
axz.axhspan(float(Cz.min()) - 2e-5, C_TH, color=C_GRN1, alpha=0.45, lw=0, zorder=0)
axz.plot(tz, Cz, color=C_HEAT, lw=2.0, zorder=4, label='$C_{\\max}$ 实算输出')
axz.scatter(tz, Cz, s=16, color=C_HEAT, marker=MK[2], zorder=5,
            edgecolor='white', linewidth=0.4)
axz.axhline(C_TH, color=C_LINE, lw=1.3, ls=':', zorder=2)
axz.axvline(t_dry_h, color=C_ACC2, lw=1.4, zorder=3)
axz.set_xlabel('时间 $t$ (h)')
axz.set_ylabel('$C_{\\max}$ (kg/kg)')
axz.set_xlim(t_dry_h - w, t_dry_h + 0.22)
axz.set_ylim(0.14994, float(Cz.max()) + 6e-5)
axz.yaxis.set_major_formatter(FuncFormatter(lambda v, _: '%.4f' % v))
axz.legend(loc='upper right', frameon=False, fontsize=FS_S,
           labelspacing=0.35, handlelength=1.5)
axz.annotate('%.7f' % float(td['C_max_at_t_dry']), xy=(tz[-1], Cz[-1]),
             xytext=(t_dry_h - 0.66, 0.15012), fontsize=FS_S, color=C_TEXT,
             ha='center',
             arrowprops=dict(arrowstyle='-|>', color=C_LINE, lw=1.0,
                             shrinkA=2, shrinkB=3))
panel(axz, 'b')

out(fig, 'fig_q3_criterion')
