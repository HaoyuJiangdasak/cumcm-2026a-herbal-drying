"""问题 1 时空 Hovmöller 图：T(r,t) 与 C(r,t) 在 (r,t) 平面的演化。

上下两条色带分别为温度与水分，横轴半径、纵轴时间；叠加等值线显示热波与
干燥前沿的推进方向。数据为 problem_1_results.json 的 hovmoller（200×21）。
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _figcommon import *  # noqa: F401,F403
import matplotlib.patheffects as pe

# clabel(inline=True) 把标签直接压在色带上，深灰字落到色图暗端时对比度不够
# （实测 C 图 '2.40' 处只有 2.57、T 图 '35' 处 3.71，WCAG AA 正文要 4.5）。
# 用白色描边把字紧邻的一圈提亮：只勾字形轮廓，不加不透明底框，字缝里的场还看得见。
HALO = [pe.withStroke(linewidth=1.8, foreground='white')]

d1 = load('problem_1_results')
hv = d1['hovmoller']
r = np.asarray(hv['radii_cm'], dtype=float)
t = np.asarray(hv['times_s'], dtype=float) / 60.0
TT = np.asarray(hv['T'], dtype=float)
CC = np.asarray(hv['C'], dtype=float)

fig = plt.figure(figsize=(6.4, 3.15), layout='constrained')
gs = fig.add_gridspec(1, 2, wspace=0.08)
axT = fig.add_subplot(gs[0, 0])
axC = fig.add_subplot(gs[0, 1])

ext = [r[0], r[-1], t[0], t[-1]]
imT = axT.imshow(TT, origin='lower', aspect='auto', extent=ext,
                 cmap=seq_cmap('heat'), interpolation='bilinear')
csT = axT.contour(r, t, TT, levels=[29, 31, 33, 35], colors=[C_LINE],
                  linewidths=0.7, alpha=0.75)
for _tx in axT.clabel(csT, inline=True, fontsize=FS_S, fmt='%.0f'):
    _tx.set_path_effects(HALO)
axT.set_xlabel('径向位置 $r$ (cm)')
axT.set_ylabel('时间 $t$ (min)')
cbT = fig.colorbar(imT, ax=axT, fraction=0.055, pad=0.03)
cbT.set_label('$T$ (℃)', fontsize=FS_M)
cbT.ax.tick_params(labelsize=FS_S)
axT.grid(False)
panel(axT, 'a')

imC = axC.imshow(CC, origin='lower', aspect='auto', extent=ext,
                 cmap=seq_cmap('moist'), interpolation='bilinear')
csC = axC.contour(r, t, CC, levels=[1.8, 2.1, 2.4, 2.52], colors=[C_LINE],
                  linewidths=0.7, alpha=0.75)
for _tx in axC.clabel(csC, inline=True, fontsize=FS_S, fmt='%.2f'):
    _tx.set_path_effects(HALO)
axC.set_xlabel('径向位置 $r$ (cm)')
axC.set_ylabel('时间 $t$ (min)')
cbC = fig.colorbar(imC, ax=axC, fraction=0.055, pad=0.03)
cbC.set_label('$C$ (kg/kg)', fontsize=FS_M)
cbC.ax.tick_params(labelsize=FS_S)
axC.grid(False)
panel(axC, 'b')

out(fig, 'fig_q1_hovmoller')
