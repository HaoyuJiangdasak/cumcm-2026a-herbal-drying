"""问题 4 收缩域水分场：(r,t) 等高线 + 移动边界 R(t) 轨迹。

实心区为药材本体，R(t) 右侧留白即已收缩掉的区域（不填 0，遵循 K14）；
虚线为干燥前沿等值线 C=0.15，与 R(t) 的交汇决定 t_dry。
数据来源 problem_4_results.json 的 result_grid / main.t_dry。
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _figcommon import *  # noqa: F401,F403

d4 = load('problem_4_results')
rg = d4['result_grid']
mn = d4['main']
td = mn['t_dry']
t_h = np.asarray(rg['times_s'], dtype=float) / 3600.0
r = np.asarray(rg['radii_cm'], dtype=float)
R_cm = np.asarray(rg['R_cm'], dtype=float)
C_raw = rg['C']
Cg = np.array([[np.nan if v is None else float(v) for v in row] for row in C_raw])
Csur = np.asarray(rg['C_surface'], dtype=float)
t_dry_h = float(td['hours'])
R_dry = float(td['R_at_t_dry_cm'])

Cm = np.ma.masked_invalid(Cg)
cmap = seq_cmap('moist')

fig = plt.figure(figsize=(6.4, 3.35), layout='constrained')
gs = fig.add_gridspec(1, 2, width_ratios=[1.34, 1.0], wspace=0.10)
ax = fig.add_subplot(gs[0, 0])
axr = fig.add_subplot(gs[0, 1])

lv = np.linspace(0, C_INIT, 15)
cf = ax.contourf(r, t_h, Cm, levels=lv, cmap=cmap, extend='neither')
ax.contour(r, t_h, Cm, levels=[0.5, 1.0, 1.5, 2.0], colors=[C_LINE],
           linewidths=0.6, alpha=0.7)
cfr = ax.contour(r, t_h, Cm, levels=[C_TH], colors=[C_HEAT], linewidths=1.8,
                 linestyles='--')
ax.plot(R_cm, t_h, color=C_SURF, lw=2.2, zorder=5, label='移动边界 $R(t)$')
ax.plot([], [], color=C_HEAT, lw=1.8, ls='--', label='干燥前沿 $C=0.15$')
ax.axhline(t_dry_h, color=C_ACC2, lw=1.4, ls='-', zorder=4,
           label='$t_{dry}$=%.4f h' % t_dry_h)
ax.set_xlabel('径向位置 $r$ (cm)')
ax.set_ylabel('时间 $t$ (h)')
ax.set_xlim(0, R0_CM)
ax.set_ylim(0, t_h[-1])
ax.grid(False)
# R(t) 右侧留白随时间才变宽，中段仅剩 27% 图宽，放不下三条中文图例；实测本 loc
# 压住 26.3% 的曲线顶点（含 R(t) 轨迹本身），故保留不透明底衬。
ax.legend(loc='center right', frameon=True, framealpha=0.92, fontsize=FS_S,
          borderpad=0.35, handlelength=1.6)
cb = fig.colorbar(cf, ax=ax, fraction=0.052, pad=0.03)
cb.set_label('$C$ (kg/kg)', fontsize=FS_M)
cb.ax.tick_params(labelsize=FS_S)
panel(ax, 'a')

cser = mn['criterion_series']
tc = np.asarray(cser['times_h'], dtype=float)
Rc = np.asarray(cser['R_cm'], dtype=float)
Cmx = np.asarray(cser['C_max'], dtype=float)
axr2 = axr.twinx()
axr.fill_between(tc, 0, Rc, color=C_SURF, alpha=0.20, lw=0, zorder=1)
lr, = axr.plot(tc, Rc, color=C_SURF, lw=2.1, zorder=4,
               label='半径 $R(t)$ → %.2f cm' % R_dry)
lc, = axr2.plot(tc, Cmx, color=C_HEAT, lw=2.1, ls='--', zorder=4,
                label='判据量 $C_{\\max}$')
axr2.axhline(C_TH, color=C_LINE, lw=1.2, ls=':', zorder=2)
axr.axvline(t_dry_h, color=C_ACC2, lw=1.4, zorder=3)
axr.set_xlabel('时间 $t$ (h)')
axr.set_ylabel('半径 $R$ (cm)', color=C_SURF_T)
axr2.set_ylabel('$C_{\\max}$ (kg/kg)', color=C_HEAT)
axr.tick_params(axis='y', colors=C_SURF_T)
axr2.tick_params(axis='y', colors=C_HEAT)
axr.set_xlim(0, tc[-1])
axr.set_ylim(0, 3.05)
axr2.set_ylim(0, 3.76)
axr.set_yticks([0.0, 0.5, 1.0, 1.5, 2.0])
axr2.set_yticks([0.0, 0.5, 1.0, 1.5, 2.0, 2.5])
axr2.grid(False)
# t_dry 落在 x 轴 99.98% 处，upper right 会被该竖线穿过；两曲线峰值仅占 68% 高度。
axr.legend(handles=[lr, lc], loc='upper left', frameon=False,
           fontsize=FS_S, labelspacing=0.35, handlelength=1.6)
panel(axr, 'b')

out(fig, 'fig_q4_contour')
