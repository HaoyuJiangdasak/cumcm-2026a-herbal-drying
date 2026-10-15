"""问题 2 三小时双轴演化：中心/表面温度与含水率 + 变物性响应。

a：T_center、T_surface 追随 T_a 升温，含水率同步下降，阴影为中心-表面温差带；
b：表面 D、ρ、Bi_m 随含水率下降的实时变化，说明变物性不可近似为常数。
数据来源 problem_2_results.json 的 dual_axis / property_evolution / radial_gap。
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _figcommon import *  # noqa: F401,F403

d2 = load('problem_2_results')
da = d2['dual_axis']
pe = d2['property_evolution']
rgap = d2['radial_gap']
t = np.asarray(da['times_s'], dtype=float) / 3600.0
Tc = np.asarray(da['T_center'], dtype=float)
Ts = np.asarray(da['T_surface'], dtype=float)
Ta = np.asarray(da['Ta'], dtype=float)
Cc = np.asarray(da['C_center'], dtype=float)
Cs = np.asarray(da['C_surface'], dtype=float)
Cb = np.asarray(da['C_bar'], dtype=float)
tp = np.asarray(pe['times_s'], dtype=float) / 3600.0
Ds = np.asarray(pe['D_surface'], dtype=float) * 1e9
Dc = np.asarray(pe['D_center'], dtype=float) * 1e9
Bim = np.asarray(pe['Bi_m_surface'], dtype=float)
peak_t = float(rgap['peak_t_h'])
peak_d = float(rgap['peak_dT'])
fin = d2['final_state']

fig = plt.figure(figsize=(6.4, 4.5), layout='constrained')
gs = fig.add_gridspec(2, 1, height_ratios=[1.20, 1.0], hspace=0.05)
ax = fig.add_subplot(gs[0, 0])
ax2 = ax.twinx()
axb = fig.add_subplot(gs[1, 0], sharex=ax)
axb2 = axb.twinx()

ax.fill_between(t, Tc, Ts, color=C_HEAT, alpha=0.16, lw=0, zorder=1)
la, = ax.plot(t, Ta, color=C_AMB, lw=1.6, ls=':', zorder=2, label='环境 $T_a$')
l1, = ax.plot(t, Ts, color=C_SURF, lw=2.0, zorder=4, label='表面 $T(R)$')
l2, = ax.plot(t, Tc, color=C_HEAT, lw=2.0, ls='--', zorder=4, label='中心 $T(0)$')
l3, = ax2.plot(t, Cs, color=C_MOIST, lw=2.0, zorder=4, label='表面 $C(R)$')
l4, = ax2.plot(t, Cc, color=C_ACC2, lw=2.0, ls='--', zorder=4, label='中心 $C(0)$')
l5, = ax2.plot(t, Cb, color=C_MEAN, lw=1.7, ls='-.', zorder=3, label='平均 $\\bar{C}$')
ax.set_ylabel('温度 $T$ (℃)')
ax2.set_ylabel('干基含水率 $C$ (kg/kg)')
ax.set_xlim(0, t[-1])
ax.set_ylim(27.0, 56.5)
ax2.set_ylim(0.86, 2.66)
ax2.grid(False)
ax.tick_params(labelbottom=False)
# 6 条曲线（3 温度 + 3 含水率）分占上下两半，双轴叠加后面板内无 2×3 大小的空白；
# 实测本 loc 下图例压住 18.1% 的曲线顶点，故保留不透明底衬。
ax.legend(handles=[la, l1, l2, l3, l4, l5], loc='center left', ncol=2,
          frameon=True, framealpha=0.92, fontsize=FS_S, borderpad=0.4,
          handlelength=1.6, columnspacing=0.9)
y_peak = 0.5 * (float(np.interp(peak_t, t, Tc)) + float(np.interp(peak_t, t, Ts)))
ax.annotate('峰值温差 %.4f ℃' % peak_d, xy=(peak_t, y_peak),
            xytext=(1.35, 54.4), fontsize=FS_S, color=C_TEXT, ha='center',
            arrowprops=dict(arrowstyle='-|>', color=C_LINE, lw=1.0,
                            shrinkA=2, shrinkB=3))
panel(ax, 'a')

axb.fill_between(tp, Ds, Dc, color=C_MOIST, alpha=0.14, lw=0, zorder=1)
lD, = axb.plot(tp, Ds, color=C_MOIST, lw=2.0, zorder=4, label='表面 $D$')
lE, = axb.plot(tp, Dc, color=C_BLUE2, lw=2.0, ls='--', zorder=4, label='中心 $D$')
lB, = axb2.plot(tp, Bim, color=C_ACC2, lw=1.9, ls='-.', zorder=4,
                label='表面 $Bi_m$')
axb.set_xlabel('时间 $t$ (h)')
axb.set_ylabel('$D$ ($10^{-9}$ m$^2$/s)', color=C_MOIST)
axb2.set_ylabel('传质毕奥数 $Bi_m$', color=C_ACC2)
axb.tick_params(axis='y', colors=C_MOIST)
axb2.tick_params(axis='y', colors=C_ACC2)
axb.set_xlim(0, t[-1])
axb.set_ylim(4.9, 13.4)
axb2.set_ylim(1.32, 3.08)
axb2.grid(False)
axb.legend(handles=[lD, lE, lB], loc='upper left', frameon=False,
           fontsize=FS_S, labelspacing=0.35, handlelength=1.6)
panel(axb, 'b')

out(fig, 'fig_q2_dual')
