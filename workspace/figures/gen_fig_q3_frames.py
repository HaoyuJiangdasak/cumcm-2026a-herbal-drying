"""问题 3 圆截面水分场快照：4 个时刻的极坐标俯视盘面。

把一维径向解按轴对称铺成圆盘，直观显示干燥前沿由外向内推进；末帧对应 t_dry。
色标统一，便于跨帧比较。数据来源 problem_3_results.json 的 result_grid。
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _figcommon import *  # noqa: F401,F403

d3 = load('problem_3_results')
rg = d3['result_grid']
td = d3['t_dry']
t_s = np.asarray(rg['times_s'], dtype=float)
r = np.asarray(rg['radii_cm'], dtype=float)
Cg = np.asarray(rg['C'], dtype=float)
t_dry_h = float(td['hours'])

want_h = [1.0, 10.0, 30.0, t_dry_h]
idx = [int(np.argmin(np.abs(t_s / 3600.0 - w))) for w in want_h]

th = np.linspace(0, 2 * np.pi, 240)
TH, RR = np.meshgrid(th, r, indexing='ij')
cmap = seq_cmap('moist')
vmin, vmax = 0.0, C_INIT

fig = plt.figure(figsize=(6.4, 2.15), layout='constrained')
gs = fig.add_gridspec(1, 4, wspace=0.06)
axes = [fig.add_subplot(gs[0, k], projection='polar') for k in range(4)]

for ax, k, tag in zip(axes, idx, 'abcd'):
    field = np.tile(Cg[k], (th.size, 1))
    pcm = ax.pcolormesh(TH, RR, field, cmap=cmap, vmin=vmin, vmax=vmax,
                        shading='auto', rasterized=True)
    ax.contour(TH, RR, field, levels=[C_TH, 0.5, 1.0, 1.8], colors=[C_LINE],
               linewidths=0.6, alpha=0.7)
    ax.set_rmax(R0_CM)
    ax.set_rticks([1.0, 2.0])
    ax.set_yticklabels([])
    ax.set_xticks([])
    ax.grid(False)
    ax.set_title('%s   %.2f h' % (tag, t_s[k] / 3600.0), loc='left',
                 fontsize=FS_M, color=C_TEXT, pad=2.5)

cb = fig.colorbar(pcm, ax=axes, fraction=0.030, pad=0.02, aspect=18)
cb.set_label('$C$ (kg/kg)', fontsize=FS_M)
cb.ax.tick_params(labelsize=FS_S)

out(fig, 'fig_q3_frames')
