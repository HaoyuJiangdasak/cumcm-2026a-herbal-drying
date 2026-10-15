"""问题 3 与问题 4 的效应分解（哑铃图）+ 两种坐标系的干物质守恒残差。

a：每行一个分支，哑铃左端固定为问题 3 的 57.5316 h，右端为该分支实算值，
连线长度即该处理带来的净效应；可见"换附录 4"单独使干燥变慢，而"半径收缩"
把时间又拉回并超过基准。b：物质坐标残差在机器精度量级，实验室坐标达 30.3%，
远超 1% 容差，故后者只作对照不作主模型。
数据来源 problem_4_results.json 的 cross_problem / shrink_isolation /
frame_comparison / dry_mass_audit。
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _figcommon import *  # noqa: F401,F403

d4 = load('problem_4_results')
cp = d4['cross_problem']
si = d4['shrink_isolation']
fc = d4['frame_comparison']
dma = d4['dry_mass_audit']

base = float(cp['t_dry_problem3_h'])
ITEMS = [('附录 4 + 收缩域\n（问题 4 主模型）', float(cp['t_dry_problem4_h'])),
         ('附录 4 + 固定域\n（隔离几何效应）', float(si['fixed_domain_h'])),
         ('附录 4 + 收缩域\n实验室坐标对照', float(fc['lab_h']))]

fig = plt.figure(figsize=(6.4, 3.3), layout='constrained')
gs = fig.add_gridspec(1, 2, width_ratios=[1.46, 1.0], wspace=0.10)
ax = fig.add_subplot(gs[0, 0])
axb = fig.add_subplot(gs[0, 1])

y = np.arange(len(ITEMS))
vals = np.array([it[1] for it in ITEMS], dtype=float)
cols = [C_HEAT if v > base else C_MOIST for v in vals]
# 哑铃杠用 hlines、两端用 scatter：与 plot 画法等效，但保持与规划配方同一套图元
lc = ax.hlines(y, base, vals, colors=cols, lw=3.0, alpha=0.5, zorder=2)
lc.set_capstyle('round')
ax.scatter(vals, y, marker=MK[0], s=72.0, c=cols, edgecolors='white',
           linewidths=1.2, zorder=4)
ax.scatter(np.full(len(y), base), y, marker=MK[1], s=56.0, c=C_AMB,
           edgecolors='white', linewidths=1.0, zorder=4)
ax.axvline(base, color=C_AMB, lw=1.2, ls='--', zorder=1)
ax.set_yticks(y)
ax.set_yticklabels([it[0] for it in ITEMS])
ax.set_xlabel('干燥终点 $t_{dry}$ (h)')
ax.set_xlim(38.0, 143.0)
ax.set_ylim(-0.62, len(ITEMS) - 0.38)
ax.grid(axis='y', visible=False)
hs = [Line2D([], [], color=C_AMB, marker=MK[1], ls='--', ms=7.5, lw=1.2),
      Line2D([], [], color=C_MOIST, marker=MK[0], ls='none', ms=8.5),
      Line2D([], [], color=C_HEAT, marker=MK[0], ls='none', ms=8.5)]
ax.legend(hs, ['问题 3 基准 %.2f h' % base, '快于基准', '慢于基准'],
          loc='lower right', frameon=False, fontsize=FS_S,
          labelspacing=0.35, handlelength=1.5)
panel(ax, 'a')

# b：干物质守恒残差（对数轴），容差 1%。
# 不用柱状：对数轴上柱长由人为设定的地板决定（这里两值真跨 12 个数量级，
# 地板取 1e-13 还是 1e-16 会让"主模型那根柱子"长度差一截）。改为以容差线
# 为基线的杆状——杆长即"离容差还有几个数量级"，是有物理含义的量。
TOL = 0.01
eps = [float(fc['eps_mass_material']), float(fc['eps_mass_lab'])]
names = ['物质坐标\n（主模型）', '实验室坐标\n（对照）']
xb = np.arange(2)
cols = [C_MEAN, C_RED2]
edges = [C_MEAN_T, C_HEAT]


def sci(v):
    e = int(np.floor(np.log10(abs(v))))
    return '$%.2f{\\times}10^{%d}$' % (v / 10.0 ** e, e)


axb.set_yscale('log')
axb.set_ylim(1e-13, 2.0)
axb.axhspan(1e-13, TOL, color=C_GRN1, alpha=0.28, lw=0, zorder=0)
axb.vlines(xb, TOL, eps, colors=edges, lw=2.6, alpha=0.7, zorder=3)
axb.scatter(xb, eps, marker=MK[0], s=88.0, c=cols, edgecolors=edges,
            linewidths=1.2, zorder=5)
axb.axhline(TOL, color=C_LINE, lw=1.3, ls=':', zorder=4,
            label='容差 $10^{-2}$')
axb.set_yticks([1e-12, 1e-9, 1e-6, 1e-3, 1e0])
axb.set_xticks(xb)
axb.set_xticklabels(names)
axb.set_xlim(-0.62, 1.62)
axb.set_ylabel('干物质守恒残差 $\\varepsilon_{mass}$')
axb.grid(axis='x', visible=False)
# 两个标注都侧放：低点贴着地板、高点贴着天花板，上下都没有放字的余量。
# 一左一右朝向面板内侧，两者纵向相距 91% 面板高，不会互相压。
axb.annotate(sci(eps[0]), xy=(0, eps[0]), xytext=(11, 0),
             textcoords='offset points', fontsize=FS_S, color=C_TEXT,
             ha='left', va='center')
axb.annotate(sci(eps[1]), xy=(1, eps[1]), xytext=(-11, 0),
             textcoords='offset points', fontsize=FS_S, color=C_TEXT,
             ha='right', va='center')
axb.legend(loc='lower right', frameon=False, fontsize=FS_S,
           labelspacing=0.35, handlelength=1.6)
panel(axb, 'b')

out(fig, 'fig_q34_compare')
