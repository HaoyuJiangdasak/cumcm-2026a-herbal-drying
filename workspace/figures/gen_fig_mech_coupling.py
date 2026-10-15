"""机理图二：控制体守恒的面积加权，与温度—含水率的双向物性耦合。

a：取 r 到 r+dr 的环形控制体，内外两个柱面的通过面积分别为 2πr 与 2π(r+dr)，
两者不等正是圆柱几何区别于平板的根源，整理后在方程里留下 (1/r)∂(r·)/∂r 的加权。
b：两条守恒律并不独立——温度经 Arrhenius 因子强烈改变扩散系数，含水率经
ρ、c_p、λ 反过来改变温度场的时间常数，构成单向强、反向弱的耦合。
本图只表达守恒律与耦合的结构，不含实算数据；耦合强弱的定量证据见灵敏度图。
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _figcommon import *  # noqa: F401,F403

FSG = 12.0
FSS = 11.0

fig = plt.figure(figsize=(6.4, 2.9), layout='constrained')
gs = fig.add_gridspec(1, 2, width_ratios=[1.0, 1.30], wspace=0.06)
ax = fig.add_subplot(gs[0, 0])
axb = fig.add_subplot(gs[0, 1])

# ── a：环形控制体与两个柱面的面积差 ──────────────────────────────────
bare(ax)
# 圆必须锁 equal aspect；数据范围按格子长宽比配平，避免锁纵横后左右留大白
ax.set_xlim(0, 10)
ax.set_ylim(0, 9.6)
ax.set_aspect('equal', adjustable='box')

cx, cy = 4.90, 4.72
R_out, r_i, r_o = 4.30, 2.16, 2.97
ax.add_patch(Circle((cx, cy), R_out, facecolor=C_GRN1, edgecolor=C_LINE,
                    lw=1.6, zorder=2))
ax.add_patch(Wedge((cx, cy), r_o, 0, 360, width=r_o - r_i, facecolor=C_RED1,
                   alpha=0.85, edgecolor=C_LINE, lw=1.3, zorder=3))
ax.plot([cx], [cy], marker=MK[1], ms=6.0, color=C_CENTER, zorder=6,
        markeredgecolor='white', markeredgewidth=1.0)

# 半径标尺：沿 -58° 引出，两个刻度落在控制体的内外面上
a_r = np.radians(-58.0)
ux, uy = np.cos(a_r), np.sin(a_r)
arrow(ax, (cx, cy), (cx + (R_out + 0.55) * ux, cy + (R_out + 0.55) * uy),
      color=C_AMB, lw=1.1, ms=10, z=5)
for rr, tag in ((r_i, '$r$'), (r_o, '$r+\\mathrm{d}r$')):
    px, py = cx + rr * ux, cy + rr * uy
    ax.plot([px - 0.30 * uy, px + 0.30 * uy], [py + 0.30 * ux, py - 0.30 * ux],
            color=C_LINE, lw=1.5, zorder=6)
# 半径标签压在标尺刻度线上：字色 #4a4a4a 与结构线 #4b4b4b 的对比度只有 1.0:1，
# 故垫不透明白底，把字下的背景固定成白色（与三维图 A/B/C 点标签同一做法）。
_wbg = dict(boxstyle='round,pad=0.16', facecolor='white', edgecolor='none',
            alpha=1.0)
ax.text(cx + r_i * ux - 0.62, cy + r_i * uy - 0.30, '$r$', fontsize=FSG,
        color=C_TEXT, ha='right', va='top', zorder=7, bbox=_wbg)
ax.text(cx + r_o * ux + 0.34, cy + r_o * uy - 0.24, '$r+\\mathrm{d}r$',
        fontsize=FSG, color=C_TEXT, ha='left', va='top', zorder=7, bbox=_wbg)

# 两个柱面的通过面积：短弧长 vs 长弧长，用弧线粗细呼应面积差
for rr, a0, a1, col, lw in ((r_i, 108, 152, C_MOIST, 2.6),
                            (r_o, 108, 152, C_HEAT, 4.0)):
    ax.add_patch(Arc((cx, cy), 2 * rr, 2 * rr, theta1=a0, theta2=a1,
                     color=col, lw=lw, zorder=7))
ax.text(0.10, 9.42, '$2\\pi (r+\\mathrm{d}r)$', fontsize=FSS, color=C_HEAT,
        ha='left', va='center', zorder=8)
ax.text(0.10, 8.30, '$2\\pi r$', fontsize=FSS, color=C_MOIST, ha='left',
        va='center', zorder=8)
ax.text(cx + 0.15, cy - r_o + 0.36, '控制体', fontsize=FSS, color=C_TEXT,
        ha='center', va='center', zorder=8)
panel(ax, 'a')

# ── b：双向耦合 ───────────────────────────────────────────────────────
bare(axb)
axb.set_xlim(0, 10)
axb.set_ylim(0, 9.6)

box(axb, (0.15, 3.90), 3.25, 2.15,
    '温度场 $T$\n$\\rho c_{p}\\,\\partial T/\\partial t$',
    C_RED1, ec=C_HEAT, fs=FSS, z=3)
box(axb, (6.60, 3.90), 3.25, 2.15,
    '含水率场 $C$\n$\\partial C/\\partial t$',
    C_GRN1, ec=C_MOIST, fs=FSS, z=3)

arrow(axb, (3.52, 5.55), (6.48, 5.55), color=C_HEAT, lw=2.6, ms=15, z=4)
arrow(axb, (6.48, 4.40), (3.52, 4.40), color=C_MOIST, lw=1.3, ms=12, z=4,
      ls='--')
# 指数写成 \exp(...) 而非 e^{...}：嵌套上下标是二级 mathtext（0.49 倍字号），
# 印刷字号补偿只除一次 0.70，补不回来，会在终检里留下 6.7pt 的字形。
axb.text(5.0, 8.45, '强耦合\n$D \\propto \\exp(-3850/T_{\\mathrm{K}})$',
         fontsize=FSS, color=C_HEAT, ha='center', va='center', zorder=6,
         linespacing=1.45)
axb.text(5.0, 1.45, '弱反馈\n$\\lambda(C),\\ \\rho(C),\\ c_{p}(C)$',
         fontsize=FSS, color=C_MOIST, ha='center', va='center', zorder=6,
         linespacing=1.45)
panel(axb, 'b')

out(fig, 'fig_mech_coupling')
