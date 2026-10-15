"""流程图六：停时判据的三候选决策树与判定时刻的先后次序。

上层列出三个候选判据，中层是各自的后果，三者比较后取逐点最大判据（实线为采
纳路径，虚线为对照路径）。底部色带只表示三种判据给出的停时由早到晚的次序，
不含数值；具体时长见问题 3 的结果图。
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _figcommon import *  # noqa: F401,F403

FSB = 13.0
FSL = 12.0

fig = plt.figure(figsize=(6.4, 3.7), layout='constrained')
ax = fig.add_subplot(111)
bare(ax)
ax.set_xlim(0.0, 10.0)
ax.set_ylim(0.0, 10.0)

HB = 1.05
XM = 5.00
XB = [1.75, 5.00, 8.25]
WC = 2.94
Y_TOP, Y_CAND, Y_VERD, Y_END = 8.55, 6.20, 4.35, 2.20
Y_SPL, Y_JOIN = 7.75, 3.75
Y_BAR = 0.95

CAND = [('表面含水率', C_RED1, C_HEAT, '--'),
        ('截面平均含水率', NP['neutral_light'], C_LINE, '--'),
        ('逐点最大含水率', C_GRN2, C_MEAN_T, '-')]
VERD = ['偏离内部状态', '低估达标时间', '全域达标可靠']

box(ax, (XM - 2.10, Y_TOP), 4.20, HB, '停时判据候选', C_RED2, ec=C_HEAT,
    fs=FSB, z=3)
for xb, (name, fc, ec, ls), vd in zip(XB, CAND, VERD):
    box(ax, (xb - WC / 2, Y_CAND), WC, HB, name, fc, ec=ec, fs=FSB, z=3)
    box(ax, (xb - WC / 2, Y_VERD), WC, HB, vd, fc, ec=ec, fs=FSB, z=3)
    flow(ax, (xb, Y_SPL), (xb, Y_CAND + HB), color=ec, lw=1.5, ms=7.0, z=5)
    flow(ax, (xb, Y_CAND), (xb, Y_VERD + HB), color=ec, lw=1.5, ms=7.0, z=5)
    ax.plot([xb, xb], [Y_VERD, Y_JOIN], color=ec, lw=1.5, ls=ls, zorder=5,
            solid_capstyle='butt')
box(ax, (XM - 2.30, Y_END), 4.60, HB, '采用逐点最大判据', C_GRN1, ec=C_MEAN_T,
    fs=FSB, z=3)

ax.plot([XM, XM], [Y_TOP, Y_SPL], color=C_LINE, lw=1.6, zorder=5,
        solid_capstyle='butt')
ax.plot([XB[0], XB[-1]], [Y_SPL, Y_SPL], color=C_LINE, lw=1.6, zorder=5,
        solid_capstyle='butt')
ax.plot([XB[0], XB[-1]], [Y_JOIN, Y_JOIN], color=C_LINE, lw=1.6, zorder=5,
        solid_capstyle='butt')
flow(ax, (XM, Y_JOIN), (XM, Y_END + HB), color=C_MEAN_T, lw=1.7, ms=7.5, z=5)

# 底部次序带：只表示先后，不带刻度数值
flow(ax, (1.30, Y_BAR), (8.90, Y_BAR), color=C_LINE, lw=1.5, ms=8.0, z=4)
for xb, (name, _f, ec, _l) in zip((2.20, 5.00, 7.80), CAND):
    ax.plot([xb], [Y_BAR], marker=MK[0], ms=6.0, ls='none', color=ec,
            markeredgecolor='white', markeredgewidth=1.0, zorder=6)
    ax.text(xb, Y_BAR + 0.52, name[:2], fontsize=FSL, color=ec, ha='center',
            va='bottom', zorder=8)
ax.text(1.30, Y_BAR - 0.42, '早', fontsize=FSL, color=C_TEXT, ha='center',
        va='top', zorder=8)
ax.text(8.90, Y_BAR - 0.42, '晚', fontsize=FSL, color=C_TEXT, ha='center',
        va='top', zorder=8)

out(fig, 'fig_criterion_logic')
