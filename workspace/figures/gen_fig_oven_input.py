"""附件 1 烘房温湿度实测：双轴时间序列 + 尾段平台延拓依据。

上：T_a(t) 与 C_a(t) 双纵轴；下：逐点一阶差分，显示 241 点末段进入恒定平台。
数据来源 user_data/附件1.xlsx，经 code/utils.load_attachment1() 读入，不做平滑。
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _figcommon import *  # noqa: F401,F403
import utils as U

t_s, Ta, Ca = (np.asarray(x, dtype=float) for x in U.load_attachment1())
t_h = t_s / 3600.0

plateau_n = 61
p0 = t_h[len(t_h) - plateau_n]
Ta_bar = float(Ta[len(Ta) - plateau_n:].mean())
Ca_bar = float(Ca[len(Ca) - plateau_n:].mean())

fig = plt.figure(figsize=(6.4, 4.5), layout='constrained')
gs = fig.add_gridspec(2, 1, height_ratios=[2.15, 1.0], hspace=0.04)
ax = fig.add_subplot(gs[0, 0])
ax2 = ax.twinx()
axd = fig.add_subplot(gs[1, 0], sharex=ax)

ax.axvspan(p0, t_h[-1], color=C_GRN1, alpha=0.55, lw=0, zorder=0)
lt, = ax.plot(t_h, Ta, color=C_HEAT, lw=2.0, zorder=3, label='空气温度 $T_a$')
ax.scatter(t_h[::12], Ta[::12], s=13, color=C_HEAT, zorder=4,
           marker=MK[2], edgecolor='white', linewidth=0.5)
lc, = ax2.plot(t_h, Ca, color=C_MOIST, lw=2.0, ls='--', zorder=3,
               label='空气水分浓度 $C_a$')

ax.axhline(Ta_bar, color=C_REF, lw=1.2, ls=':', zorder=1)
ax.set_ylabel('空气温度 $T_a$ (℃)', color=C_HEAT)
ax2.set_ylabel('空气水分浓度 $C_a$ (kg/kg)', color=C_MOIST)
ax.tick_params(axis='y', colors=C_HEAT)
ax2.tick_params(axis='y', colors=C_MOIST)
ax.set_ylim(26.5, 52.0)
ax2.set_ylim(0.016, 0.055)
ax2.grid(False)
ax.set_xlim(0, t_h[-1])
ax.tick_params(labelbottom=False)
ax.legend(handles=[lt, lc], loc='lower right', frameon=False,
          fontsize=FS_S, labelspacing=0.35, handlelength=1.6)
ax.annotate('平台 $T_a$=%.4f ℃\n$C_a$=%.6f' % (Ta_bar, Ca_bar),
            xy=(p0, Ta_bar), xytext=(1.25, 41.0), fontsize=FS_S, color=C_TEXT,
            arrowprops=dict(arrowstyle='-|>', color=C_LINE, lw=1.1,
                            shrinkA=2, shrinkB=3))
panel(ax, 'a', dx=-0.105, dy=1.06)

dT = np.diff(Ta)
dC = np.diff(Ca)
tc = t_h[1:]
axd.axvspan(p0, t_h[-1], color=C_GRN1, alpha=0.55, lw=0, zorder=0)
axd.axhline(0, color=C_REF, lw=1.0, zorder=1)
axd.plot(tc, dT, color=C_HEAT, lw=1.3, zorder=3, label='$\\Delta T_a$ (℃)')
axd.plot(tc, dC * 1e3, color=C_MOIST, lw=1.3, ls='--', zorder=3,
         label='$\\Delta C_a\\times 10^{3}$')
axd.set_xlabel('时间 $t$ (h)')
axd.set_ylabel('相邻点增量')
axd.set_xlim(0, t_h[-1])
axd.set_ylim(-0.62, 1.42)
axd.legend(loc='upper right', ncol=2, frameon=False, fontsize=FS_S)
panel(axd, 'b', dx=-0.105, dy=1.16)

out(fig, 'fig_oven_input')
