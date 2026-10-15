"""流程图三：问题 2 全过程求解的双重迭代结构。

主链自上而下为一次时间步内的动作；右侧两条回路分别是步内的 Picard 非线性
迭代与步间的时间推进，回路样式在图例中区分。只给环节名称，不进公式与数值。
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from matplotlib.lines import Line2D
from _figcommon import *  # noqa: F401,F403

FSB = 13.0
FSL = 12.0

fig = plt.figure(figsize=(6.4, 4.6), layout='constrained')
ax = fig.add_subplot(111)
bare(ax)
# 坐标范围贴着墨迹裁：主链左边界在 1.30、末框底边在 0.96，留到 0 会在左侧和
# 底部各浪费一成图幅
ax.set_xlim(1.15, 10.0)
ax.set_ylim(0.82, 9.86)

XA, WB = 1.30, 4.60             # 主链左边界与框宽
HB, GAP = 1.24, 0.26
YS = [8.46 - i * (HB + GAP) for i in range(6)]
xc = XA + WB / 2

STEPS = [
    ('初始温度与含水率场', C_GRN1, C_MEAN_T),
    ('推进到下一时刻', NP['neutral_light'], C_LINE),
    ('按现值组装隐式系数', C_RED1, C_HEAT),
    ('三对角直接求解两场', C_RED1, C_HEAT),
    ('迭代收敛判定', C_GRN2, C_MEAN_T),
    ('输出全过程时序', C_GRN1, C_MEAN_T),
]
for yi, (txt, fc, ec) in zip(YS, STEPS):
    box(ax, (XA, yi), WB, HB, txt, fc, ec=ec, fs=FSB, z=3)
for k in range(len(YS) - 1):
    flow(ax, (xc, YS[k]), (xc, YS[k + 1] + HB), color=C_LINE, lw=1.6, ms=7.0,
         z=5)

# 两条回路走主链右侧的独立通道，回到各自的起点框
XI, XO = 6.55, 9.10
y_judge = YS[4] + HB / 2
elbow(ax, [(XA + WB, y_judge - 0.28), (XI, y_judge - 0.28),
           (XI, YS[2] + HB / 2), (XA + WB, YS[2] + HB / 2)],
      color=C_HEAT, lw=1.5, ms=7.0)
elbow(ax, [(XA + WB, y_judge + 0.28), (XO, y_judge + 0.28),
           (XO, YS[1] + HB / 2), (XA + WB, YS[1] + HB / 2)],
      color=C_MEAN_T, lw=1.5, ms=7.0, ls='--')

hs = [Line2D([], [], color=C_HEAT, lw=1.5),
      Line2D([], [], color=C_MEAN_T, lw=1.5, ls='--')]
ax.legend(hs, ['步内 Picard 迭代', '步间时间推进'], loc='upper right',
          frameon=False, fontsize=FSL, handlelength=1.4, handletextpad=0.5,
          borderpad=0.2, labelspacing=0.5)

out(fig, 'fig_flow_q2')
