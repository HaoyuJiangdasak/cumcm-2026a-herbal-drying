"""流程图五：问题 4 收缩域的两条坐标路线与合流择优。

共同输入后分成左右两支：左支把物理域映射到固定域，右支保留物理域而让网格随
时间移动。两支各自算完后由干物质守恒审计合流，据此选定主模型。
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _figcommon import *  # noqa: F401,F403

FSB = 13.0

fig = plt.figure(figsize=(6.4, 4.2), layout='constrained')
ax = fig.add_subplot(111)
bare(ax)
# 坐标范围贴着墨迹裁：两支外缘在 0.60 与 9.40、末框底边在 0.90
ax.set_xlim(0.46, 9.54)
ax.set_ylim(0.74, 9.76)

HB = 1.15
XM = 5.00
XL, XR = 2.45, 7.55              # 两支的中线
WC = 3.70                        # 两支的框宽
Y_TOP, Y_R1, Y_R2, Y_MRG, Y_END = 8.45, 6.30, 4.70, 2.60, 0.90
Y_SPL, Y_JOIN = 7.95, 4.20       # 分流与合流的横向通道

box(ax, (XM - 2.30, Y_TOP), 4.60, HB, '读附件 2 半径收缩序列',
    NP['neutral_light'], ec=C_LINE, fs=FSB, z=3)

LEFT = ['物质坐标固定网格', '收缩进入扩散项']
RIGHT = ['实验室坐标动网格', '表观对流项显式']
for yi, tl, tr in zip((Y_R1, Y_R2), LEFT, RIGHT):
    box(ax, (XL - WC / 2, yi), WC, HB, tl, C_GRN2, ec=C_MEAN_T, fs=FSB, z=3)
    box(ax, (XR - WC / 2, yi), WC, HB, tr, C_RED1, ec=C_HEAT, fs=FSB, z=3)

box(ax, (XM - 1.90, Y_MRG), 3.80, HB, '干物质守恒审计', C_RED2, ec=C_HEAT,
    fs=FSB, z=3)
box(ax, (XM - 2.10, Y_END), 4.20, HB, '选定物质坐标主模型', C_GRN1,
    ec=C_MEAN_T, fs=FSB, z=3)

# 分流：先落到横向通道，再分别下扎进两支
ax.plot([XM, XM], [Y_TOP, Y_SPL], color=C_LINE, lw=1.6, zorder=5,
        solid_capstyle='butt')
ax.plot([XL, XR], [Y_SPL, Y_SPL], color=C_LINE, lw=1.6, zorder=5,
        solid_capstyle='butt')
for xb, col in ((XL, C_MEAN_T), (XR, C_HEAT)):
    flow(ax, (xb, Y_SPL), (xb, Y_R1 + HB), color=col, lw=1.6, ms=7.0, z=5)
    flow(ax, (xb, Y_R1), (xb, Y_R2 + HB), color=col, lw=1.6, ms=7.0, z=5)
    ax.plot([xb, xb], [Y_R2, Y_JOIN], color=col, lw=1.6, zorder=5,
            solid_capstyle='butt')

# 合流：两支汇到横向通道后单箭头进审计框
ax.plot([XL, XR], [Y_JOIN, Y_JOIN], color=C_LINE, lw=1.6, zorder=5,
        solid_capstyle='butt')
flow(ax, (XM, Y_JOIN), (XM, Y_MRG + HB), color=C_LINE, lw=1.7, ms=7.5, z=5)
flow(ax, (XM, Y_MRG), (XM, Y_END + HB), color=C_LINE, lw=1.7, ms=7.5, z=5)

out(fig, 'fig_flow_q4')
