"""流程图四：问题 3 停止时刻判定的分支与插值定位。

主链逐时刻推进，判据框有两个出口：未达阈值则回到下一时刻，已达阈值则在跨越
区间内线性插值定位停时。分支出口用两字锚标注，其余说明留给图题。
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _figcommon import *  # noqa: F401,F403

FSB = 13.0
FSL = 12.0

fig = plt.figure(figsize=(6.4, 4.5), layout='constrained')
ax = fig.add_subplot(111)
bare(ax)
# 坐标范围只包住主链与回路通道：右侧留到 10 会空出近两成图幅，而回路最远只到
# 7.30。判据框下方的间隔特意留宽，用来放"已达"的分支锚。
ax.set_xlim(0.45, 7.45)
ax.set_ylim(0.30, 9.70)

XA, WB, HB = 0.60, 4.90, 1.40
xc, XR = XA + WB / 2, XA + WB
YS = [8.20, 6.38, 4.56, 2.30, 0.48]

STEPS = [
    ('载入全过程含水率场', C_GRN1, C_MEAN_T),
    ('逐时刻取径向最大值', C_GRN2, C_MEAN_T),
    ('与停止阈值比较', C_RED2, C_HEAT),
    ('跨越区间线性插值', C_RED1, C_HEAT),
    ('输出停止时刻与位置', C_GRN1, C_MEAN_T),
]
for yi, (txt, fc, ec) in zip(YS, STEPS):
    box(ax, (XA, yi), WB, HB, txt, fc, ec=ec, fs=FSB, z=3)
for k in range(len(YS) - 1):
    flow(ax, (xc, YS[k]), (xc, YS[k + 1] + HB), color=C_LINE, lw=1.6, ms=7.0,
         z=5)

# 未达阈值：回到上一环节取下一时刻
XL = 7.30
y_dec, y_scan = YS[2] + HB / 2, YS[1] + HB / 2
elbow(ax, [(XR, y_dec), (XL, y_dec), (XL, y_scan), (XR, y_scan)],
      color=C_MEAN_T, lw=1.5, ms=7.0, ls='--')
# 锚标放在主链与回路通道之间的空档，两侧各留半格净空
ax.text(0.5 * (XR + XL), 0.5 * (y_dec + y_scan), '未达', fontsize=FSL,
        color=C_MEAN_T, ha='center', va='center', zorder=8)
ax.text(xc + 0.22, 0.5 * (YS[2] + YS[3] + HB), '已达', fontsize=FSL,
        color=C_HEAT, ha='left', va='center', zorder=8)

out(fig, 'fig_flow_q3')
