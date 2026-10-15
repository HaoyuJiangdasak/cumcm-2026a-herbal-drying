"""流程图二：问题 1 预热平衡的求解流程，按三条纵向泳道组织。

左中右为准备、离散、判定三个阶段，泳道内自上而下为执行顺序，泳道之间的横向
箭头表示阶段交接。只给环节名称，公式与数值不进图。
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _figcommon import *  # noqa: F401,F403

FSB = 13.0
FSL = 12.0

fig = plt.figure(figsize=(6.4, 3.4), layout='constrained')
ax = fig.add_subplot(111)
bare(ax)
ax.set_xlim(0.0, 10.0)
# 纵向范围贴着泳道底带裁：留到 0 会在图底留下近 20% 的空白
ax.set_ylim(1.80, 9.98)

XC = [1.72, 5.00, 8.28]           # 三条泳道的中线
WB = 3.02                         # 框宽：留足与泳道底带的净空
Y_HEAD, H_HEAD = 8.62, 0.94
YS = [6.62, 4.42, 2.22]
HB = 1.56

LANES = [
    ('准备', C_GRN1, C_MEAN_T, NP['green_1'],
     ['读附件 1 烘房条件', '取附录 2 定常物性', '设定初始均匀温度']),
    ('离散', C_RED1, C_HEAT, NP['red_1'],
     ['径向控制体划分', '表面对流边界', '全隐式时间推进']),
    ('判定', NP['neutral_light'], C_LINE, NP['neutral_light'],
     ['中心温度稳定判定', '给出预热时长', '输出温度径向剖面']),
]

for xc, (name, fc, ec, band, items) in zip(XC, LANES):
    # 泳道底带：把同一阶段的环节在视觉上圈到一起
    ax.add_patch(Rectangle((xc - WB / 2 - 0.16, YS[-1] - 0.34),
                           WB + 0.32, Y_HEAD + H_HEAD - YS[-1] + 0.68,
                           facecolor=band, alpha=0.30, edgecolor='none',
                           zorder=1))
    box(ax, (xc - WB / 2, Y_HEAD), WB, H_HEAD, name, fc, ec=ec, fs=FSL, z=3)
    for yi, it in zip(YS, items):
        box(ax, (xc - WB / 2, yi), WB, HB, it, fc, ec=ec, fs=FSB, z=3)
    for k in range(len(YS) - 1):
        flow(ax, (xc, YS[k]), (xc, YS[k + 1] + HB), color=ec, lw=1.5, ms=7.0,
             z=5)

# 阶段交接：横向箭头落在泳道之间的空档，取中间一行的高度
y_mid = YS[1] + HB / 2
for a, b in zip(XC[:-1], XC[1:]):
    flow(ax, (a + WB / 2 + 0.20, y_mid), (b - WB / 2 - 0.20, y_mid),
         color=C_LINE, lw=1.6, ms=7.5, z=5)

out(fig, 'fig_flow_q1')
