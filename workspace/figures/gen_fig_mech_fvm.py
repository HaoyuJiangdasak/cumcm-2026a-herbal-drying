"""机理图三：顶点式有限体积网格与边界的自然嵌入，及其三对角代数结构。

a：节点 $s_i$ 落在输出要求的位置上，控制体取相邻节点的中点连线，两端为半控制体。
内侧界面的面积因子 s=0 使轴线零通量精确满足，不需虚节点；外侧界面即物理表面，
Robin 条件以通量形式直接进入方程。b：由此得到的系数矩阵严格三对角，首末行由
两个边界条件闭合，用 Thomas 算法按 O(N) 直接求解。
本图为离散格式的结构示意，网格数取 8 仅为作图清晰，不代表实算网格。
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from matplotlib.colors import ListedColormap
from _figcommon import *  # noqa: F401,F403

FSG = 12.0
FSS = 11.0
NG = 8                      # 示意网格数，非实算网格

fig = plt.figure(figsize=(6.4, 2.75), layout='constrained')
gs = fig.add_gridspec(1, 2, width_ratios=[1.62, 1.0], wspace=0.06)
ax = fig.add_subplot(gs[0, 0])
axb = fig.add_subplot(gs[0, 1])

# ── a：网格与控制体 ───────────────────────────────────────────────────
bare(ax)
ax.set_xlim(0.0, 10.0)
ax.set_ylim(0.0, 9.0)

xs = np.linspace(0.85, 9.15, NG + 1)
dx = xs[1] - xs[0]
y0, y1 = 3.55, 6.05
face = np.concatenate(([xs[0]], 0.5 * (xs[:-1] + xs[1:]), [xs[-1]]))
IG = 4                      # 被高亮的一般控制体

for i in range(NG + 1):
    if i == IG:
        fc, ec, lw = C_GRN2, C_MEAN_T, 1.6
    elif i in (0, NG):
        fc, ec, lw = C_RED1, C_HEAT, 1.6
    else:
        fc, ec, lw = 'white', C_LINE, 0.9
    ax.add_patch(Rectangle((face[i], y0), face[i + 1] - face[i], y1 - y0,
                           facecolor=fc, edgecolor=ec, lw=lw, zorder=3))
ax.plot([0.35, 9.65], [y0, y0], color=C_LINE, lw=1.5, zorder=4)
ax.plot(xs, np.full_like(xs, y0), marker=MK[0], ms=5.5, ls='none',
        color=C_CENTER, zorder=6, markeredgecolor='white', markeredgewidth=1.0)

# 一般控制体的两个界面通量
for xf, sgn in ((face[IG], 1.0), (face[IG + 1], 1.0)):
    arrow(ax, (xf - 0.62 * sgn, 4.80), (xf + 0.30 * sgn, 4.80), color=C_MOIST,
          lw=1.7, ms=12, z=7)
ax.text(face[IG] - 0.10, 7.20, '$j_{i-1/2}$', fontsize=FSS, color=C_MOIST,
        ha='center', va='center', zorder=8)
ax.text(face[IG + 1] + 0.72, 7.20, '$j_{i+1/2}$', fontsize=FSS, color=C_MOIST,
        ha='center', va='center', zorder=8)

# 两端边界：轴线处面积因子为零；表面处 Robin 进入通量项
ax.text(xs[0], 7.20, '$s=0$', fontsize=FSS, color=C_HEAT, ha='left',
        va='center', zorder=8)
ax.text(xs[-1] + 0.42, 7.20, '$h,\\ h_{m}$', fontsize=FSS, color=C_HEAT,
        ha='right', va='center', zorder=8)

# 节点标注与步长
for i, tag in ((0, '$s_{0}$'), (IG, '$s_{i}$'), (NG, '$s_{N}$')):
    ax.text(xs[i], y0 - 0.62, tag, fontsize=FSG, color=C_TEXT, ha='center',
            va='top', zorder=8)
arrow(ax, (xs[1], 1.42), (xs[2], 1.42), color=C_AMB, lw=1.2, style='<|-|>',
      ms=10, z=6)
ax.text(0.5 * (xs[1] + xs[2]), 0.72, '$\\Delta s$', fontsize=FSG, color=C_TEXT,
        ha='center', va='top', zorder=8)
panel(ax, 'a')

# ── b：三对角结构 ─────────────────────────────────────────────────────
n = NG + 1
M = np.zeros((n, n), dtype=float)
for i in range(n):
    M[i, i] = 2.0
    if i > 0:
        M[i, i - 1] = 1.0
    if i < n - 1:
        M[i, i + 1] = 1.0
for i in (0, n - 1):                      # 首末行由两个边界条件闭合
    M[i, M[i] > 0] = 3.0

cmap = ListedColormap(['white', C_GRN1, C_GRN2, C_RED2])
axb.imshow(M, cmap=cmap, vmin=-0.5, vmax=3.5, origin='upper',
           aspect='equal', interpolation='nearest')
axb.set_xticks([0, IG, n - 1])
axb.set_xticklabels(['$0$', '$i$', '$N$'], fontsize=FSS)
axb.set_yticks([0, IG, n - 1])
axb.set_yticklabels(['$0$', '$i$', '$N$'], fontsize=FSS)
axb.set_xticks(np.arange(-0.5, n, 1), minor=True)
axb.set_yticks(np.arange(-0.5, n, 1), minor=True)
axb.grid(which='minor', color=C_REF, linewidth=0.7)
axb.grid(which='major', visible=False)
axb.tick_params(which='minor', length=0)
axb.tick_params(which='major', length=2.5, labelsize=FSS)
hs = [Rectangle((0, 0), 1, 1, facecolor=C_GRN2, edgecolor=C_MEAN_T, lw=0.9),
      Rectangle((0, 0), 1, 1, facecolor=C_GRN1, edgecolor=C_MEAN_T, lw=0.9),
      Rectangle((0, 0), 1, 1, facecolor=C_RED2, edgecolor=C_HEAT, lw=0.9)]
axb.legend(hs, ['主对角', '次对角', '边界行'], loc='upper center',
           bbox_to_anchor=(0.5, -0.13), ncol=3, frameon=False, fontsize=FSS,
           handlelength=1.1, columnspacing=0.9, handletextpad=0.4,
           borderpad=0.1)
panel(axb, 'b')

out(fig, 'fig_mech_fvm')
