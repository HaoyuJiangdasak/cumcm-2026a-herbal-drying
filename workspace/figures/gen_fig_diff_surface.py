"""附录 3 扩散系数 D(C,T) 三维曲面 + 实算状态轨迹俯视投影。

左：lg D 在 (C,T) 网格上的曲面，A/B/C 三点为问题 3 生产算例经过的真实状态；
右：同一函数的等高线俯视图，叠加问题 3 中心与表面状态点的实算轨迹。
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from mpl_toolkits.mplot3d import Axes3D  # noqa: F401
from _figcommon import *  # noqa: F401,F403
import utils as U

pr = U.Props(3)
Cg = np.linspace(0.05, C_INIT, 90)
Tg = np.linspace(T_INIT, 50.5, 70)
CC, TT = np.meshgrid(Cg, Tg)
DD = np.asarray(pr.D(CC.ravel(), TT.ravel()), dtype=float).reshape(CC.shape)
LD = np.log10(DD)

d3 = load('problem_3_results')
cs = d3['criterion_series']
C_surf = np.asarray(cs['C_surface'], dtype=float)
T_surf = np.asarray(cs['T_surface'], dtype=float)
C_cen = np.asarray(cs['C_max'], dtype=float)
T_cen = np.asarray(cs['T_center'], dtype=float)

# 3D 轴使整图跳过自动字号提升，故本图字号显式抬高，留足自动排版的内缩余量
FS3 = 12.5
plt.rcParams['xtick.labelsize'] = FS3
plt.rcParams['ytick.labelsize'] = FS3
plt.rcParams['axes.labelsize'] = FS3
fig = plt.figure(figsize=(5.85, 3.05), layout='constrained')
gs = fig.add_gridspec(1, 2, width_ratios=[1.12, 1.0], wspace=0.10)
ax3 = fig.add_subplot(gs[0, 0], projection='3d')
ax2 = fig.add_subplot(gs[0, 1])

cmap = seq_cmap('moist')
ax3.plot_surface(CC, TT, LD, cmap=cmap, linewidth=0, antialiased=True,
                 rstride=1, cstride=1, alpha=0.95, shade=True)
z_off = float(LD.min()) - 0.7
ax3.contour(CC, TT, LD, levels=8, zdir='z', offset=z_off,
            colors=[C_LINE], linewidths=0.55, alpha=0.6)
ax3.set_xlabel('C (kg/kg)', labelpad=-2, fontsize=FS3)
ax3.set_ylabel('T (℃)', labelpad=-2, fontsize=FS3)
ax3.set_zlabel('lg D', labelpad=1.5, fontsize=FS3)
ax3.set_zlim(z_off, float(LD.max()) + 0.3)
ax3.tick_params(labelsize=FS3, pad=-2.0)
ax3.set_xticks([0, 1, 2])
ax3.set_yticks([30, 40, 50])
ax3.view_init(elev=23, azim=-128)
ax3.set_box_aspect((1.05, 1.0, 0.60))
ax3.grid(False)

pts = [(C_INIT, T_INIT, 'A'), (C_TH, 49.9989, 'B'), (0.0526, 49.9989, 'C')]
for cv, tv, tag in pts:
    zv = float(np.log10(np.asarray(pr.D(np.array([cv]), np.array([tv])))[0]))
    ax3.scatter([cv], [tv], [zv], s=32, color=C_HEAT, marker=MK[1],
                depthshade=False, zorder=6, edgecolor='white', linewidth=0.6)
    ax3.text(cv, tv, zv + 0.75, tag, fontsize=FS3, color=C_HEAT,
             fontweight='bold', ha='center', zorder=8,
             bbox=dict(boxstyle='round,pad=0.16', facecolor='white',
                       edgecolor='none', alpha=1.0))
panel(ax3, 'a')

cf = ax2.contourf(CC, TT, LD, levels=14, cmap=cmap, alpha=0.9)
ax2.contour(CC, TT, LD, levels=7, colors=[C_LINE], linewidths=0.6, alpha=0.75)
ax2.plot(C_surf, T_surf, color=C_SURF, lw=2.0, zorder=5, label='表面状态轨迹')
ax2.plot(C_cen, T_cen, color=C_HEAT, lw=2.0, ls='--', zorder=5, label='中心状态轨迹')
lab_xy = {'A': (1.92, 31.8), 'B': (0.72, 45.4), 'C': (0.26, 38.2)}
for cv, tv, tag in pts:
    ax2.scatter([cv], [tv], s=30, color=C_HEAT, marker=MK[1], zorder=6,
                edgecolor='white', linewidth=0.6)
    ax2.annotate(tag, xy=(cv, tv), xytext=lab_xy[tag], fontsize=FS_M,
                 fontweight='bold', color=C_HEAT, ha='center', va='center',
                 zorder=8,
                 bbox=dict(boxstyle='round,pad=0.16', facecolor='white',
                           edgecolor='none', alpha=1.0),
                 arrowprops=dict(arrowstyle='-', color=C_HEAT, lw=0.9,
                                 shrinkA=4, shrinkB=3, alpha=0.8))
ax2.set_xlabel('干基含水率 C (kg/kg)')
ax2.set_ylabel('药材温度 T (℃)')
ax2.set_xlim(0, C_INIT)
ax2.set_ylim(T_INIT - 0.6, 50.8)
ax2.grid(False)
# 本图 b 是满幅 contourf（14 级填充铺满整个面板），没有任何空白角落可放图例；
# 去框后文字直接压在中调色块上不可读，故保留不透明底衬——可读性优先于去框风格。
ax2.legend(loc='lower right', frameon=True, framealpha=0.92, fontsize=FS_S)
cb = fig.colorbar(cf, ax=ax2, fraction=0.055, pad=0.03)
cb.set_label('lg D  (D 单位 m²/s)', fontsize=FS_M)
cb.ax.tick_params(labelsize=FS_S)
panel(ax2, 'b')

out(fig, 'fig_diff_surface')
