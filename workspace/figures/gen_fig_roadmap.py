"""流程图一：全文技术路线的分层框图。

自上而下五层——题面数据、公共机理、数值方法、四个问题、可信性与择优。
同层左右并列表示并行处理，层间箭头表示依赖方向。流程图只给环节名称，
公式与数值一律不进图，避免小字。
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _figcommon import *  # noqa: F401,F403

FSB = 13.0      # 框内字号：流程图字号比数据图再大一档
FSL = 12.0      # 层名字号

fig = plt.figure(figsize=(6.4, 4.6), layout='constrained')
ax = fig.add_subplot(111)
bare(ax)
ax.set_xlim(0.0, 10.0)
ax.set_ylim(0.0, 10.1)

# 层名列压到两字、内容列拉到最宽：框内文字与框边的净空不足会被判"文字与柱边界遮挡"
XL, XR = 1.05, 9.88
HB = 1.52
YS = [8.44, 6.36, 4.28, 2.20, 0.12]

LAYERS = [
    ('数据', C_GRN1, C_MEAN_T,
     ['附件 1 烘房条件', '附件 2 半径收缩', '附录 2/3/4 物性']),
    ('机理', C_RED1, C_HEAT,
     ['圆柱降维', '热质双向耦合', 'Robin 边界']),
    ('方法', NP['neutral_light'], C_LINE,
     ['顶点式有限体积', '全隐式 + Picard', '三对角直接解']),
    ('求解', C_GRN2, C_MEAN_T,
     ['问题 1\n预热平衡', '问题 2\n全过程', '问题 3\n停时判定', '问题 4\n收缩域']),
    ('校验', C_RED2, C_HEAT,
     ['网格与时间步收敛', '解析解交叉校验', '灵敏度与择优']),
]

for yi, (name, fc, ec, items) in zip(YS, LAYERS):
    ax.text(XL - 0.28, yi + HB / 2, name, fontsize=FSL, color=ec, ha='right',
            va='center', zorder=6)
    n = len(items)
    gap = 0.24
    wb = (XR - XL - gap * (n - 1)) / n
    for j, it in enumerate(items):
        box(ax, (XL + j * (wb + gap), yi), wb, HB, it, fc, ec=ec, fs=FSB, z=3)

# 层间依赖箭头：走内容列中线，箭头落在层间空档
xm = 0.5 * (XL + XR)
for k in range(len(YS) - 1):
    flow(ax, (xm, YS[k]), (xm, YS[k + 1] + HB), color=C_LINE, lw=1.7, ms=8.0,
         z=5)

out(fig, 'fig_roadmap')
