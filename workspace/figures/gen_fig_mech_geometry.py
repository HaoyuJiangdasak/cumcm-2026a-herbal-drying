"""机理图一：几何降维与边界条件。

a：细长圆柱药材，轴向长度远大于半径且两端封闭，故轴向梯度可忽略、周向轴对称，
三维问题降为沿半径 r 的一维问题。b：任一横截面上的 Robin 边界——热量由热空气
流入表面，水分由表面流出至空气，两者共用同一表面温度与表面含水率。
本图不含实算数据，仅表达题面几何与守恒律的结构。
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _figcommon import *  # noqa: F401,F403

FSG = 12.0     # 机理图用字号（比数据图大一档，缩放后仍清晰）
FSS = 11.0

fig = plt.figure(figsize=(6.4, 2.85), layout='constrained')
gs = fig.add_gridspec(1, 2, width_ratios=[1.16, 1.0], wspace=0.06)
ax = fig.add_subplot(gs[0, 0])
axb = fig.add_subplot(gs[0, 1])

# ── a：圆柱与降维 ──────────────────────────────────────────────────────
bare(ax)
# 本panel 不锁 equal aspect：锁了会让方形绘图区在宽格子里左右留大白，
# 而这里只是示意圆柱，轻微纵横拉伸不影响表达（panel b 的圆必须锁）。
ax.set_xlim(0.95, 9.95)
ax.set_ylim(0.85, 8.55)

cx0, cx1, cyc, ry = 1.7, 8.0, 6.4, 1.55
ax.add_patch(Rectangle((cx0, cyc - ry), cx1 - cx0, 2 * ry, facecolor=C_GRN1,
                       edgecolor='none', zorder=2))
ax.add_patch(Ellipse((cx1, cyc), 0.86, 2 * ry, facecolor=C_GRN2,
                     edgecolor=C_LINE, lw=1.5, zorder=4))
ax.add_patch(Ellipse((cx0, cyc), 0.86, 2 * ry, facecolor=C_GRN1,
                     edgecolor=C_LINE, lw=1.5, ls='--', zorder=1))
for sgn in (-1, 1):
    ax.plot([cx0, cx1], [cyc + sgn * ry, cyc + sgn * ry], color=C_LINE,
            lw=1.5, zorder=3)
ax.plot([cx0 - 0.9, cx1 + 0.9], [cyc, cyc], color=C_AMB, lw=1.1, ls='-.',
        zorder=5)
arrow(ax, (cx1, cyc), (cx1, cyc + ry), color=C_MOIST, lw=1.7, z=6)
ax.text(cx1 + 0.30, cyc + 0.80 * ry, '$R_{0}$', fontsize=FSG, color=C_MOIST,
        ha='left', va='center', zorder=7)
ax.text(cx0 + 1.5, cyc + 0.72 * ry, '药材本体', fontsize=FSS, color=C_TEXT,
        ha='center', va='center', zorder=7)
arrow(ax, (cx0 + 0.2, 4.42), (cx1, 4.42), color=C_TEXT, lw=1.3,
      style='<|-|>', z=6)
ax.text(cx1 + 0.34, 4.42, '$L \\gg R_{0}$', fontsize=FSG, color=C_TEXT,
        ha='left', va='center', zorder=7)

# 降维箭头 + 结果：一维半径线
arrow(ax, (5.0, 3.60), (5.0, 2.60), color=C_ACC2, lw=2.0, ms=13, z=6)
ax.plot([2.7, 7.7], [1.95, 1.95], color=C_LINE, lw=2.2, zorder=5)
for xv, tag, col in ((2.7, '$r=0$', C_CENTER), (7.7, '$r=R_{0}$', C_SURF_T)):
    ax.plot([xv], [1.95], marker=MK[0], ms=8.0, color=col, zorder=6,
            markeredgecolor='white', markeredgewidth=1.2)
    ax.text(xv, 1.28, tag, fontsize=FSG, color=col, ha='center', va='top',
            zorder=7)
panel(ax, 'a')

# ── b：横截面与 Robin 边界 ────────────────────────────────────────────
bare(axb)
axb.set_xlim(0, 10)
axb.set_ylim(0, 10)
axb.set_aspect('equal', adjustable='box')

ox, oy, Rr = 5.0, 5.35, 2.55
axb.add_patch(Circle((ox, oy), Rr * 1.42, facecolor=C_RED1, alpha=0.30,
                     edgecolor='none', zorder=1))
for k, a in enumerate(np.linspace(1.0, 0.28, 5)):
    axb.add_patch(Circle((ox, oy), Rr * a, facecolor='none',
                         edgecolor=C_REF, lw=0.9, ls=':', zorder=2))
axb.add_patch(Circle((ox, oy), Rr, facecolor=C_GRN1, edgecolor=C_LINE,
                     lw=1.8, zorder=3))
axb.plot([ox], [oy], marker=MK[1], ms=7.0, color=C_CENTER, zorder=6,
         markeredgecolor='white', markeredgewidth=1.0)

# 热流入（红，指向内）与湿流出（蓝，指向外），角度错开避免压线
for a_deg in (128.0, 168.0):
    a = np.radians(a_deg)
    p1 = (ox + Rr * np.cos(a), oy + Rr * np.sin(a))
    p0 = (ox + Rr * 1.44 * np.cos(a), oy + Rr * 1.44 * np.sin(a))
    arrow(axb, p0, p1, color=C_HEAT, lw=1.9, ms=12, z=7)
for a_deg in (-52.0, -12.0):
    a = np.radians(a_deg)
    p0 = (ox + Rr * np.cos(a), oy + Rr * np.sin(a))
    p1 = (ox + Rr * 1.44 * np.cos(a), oy + Rr * 1.44 * np.sin(a))
    arrow(axb, p0, p1, color=C_MOIST, lw=1.9, ms=12, z=7)

# 下标一律写成 _{x} 而非 _x：落盘期的印刷字号补偿只认带花括号的 mathtext
# token，写 _x 会让 0.7 倍的下标字形绕过补偿、在终检里跌破 8pt。
axb.text(0.30, 9.35, '$-\\lambda\\,\\partial T/\\partial r = h\\,(T_{R}-T_{a})$',
         fontsize=FSS, color=C_HEAT, ha='left', va='center', zorder=8)
axb.text(9.70, 1.05, '$-D\\,\\partial C/\\partial r = h_{m} (C_{R}-C_{a})$',
         fontsize=FSS, color=C_MOIST, ha='right', va='center', zorder=8)
axb.text(ox, oy - 0.62, '$\\partial T/\\partial r=0$', fontsize=FSS,
         color=C_CENTER, ha='center', va='top', zorder=8)
axb.text(8.55, 8.30, '$T_{a},\\ C_{a}$', fontsize=FSG, color=C_AMB,
         ha='center', va='center', zorder=8)
panel(axb, 'b')

out(fig, 'fig_mech_geometry')
