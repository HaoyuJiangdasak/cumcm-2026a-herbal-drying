"""机理图四：物质坐标 s=r/R(t) 的映射，与收缩带来的 1/R² 放大。

a：物质坐标把随时间收缩的物理域 [0,R(t)] 映射到固定域 [0,1]。节点固定在 s 上，
随物质一起向轴线移动，故方程中不出现 Ṙ 对流项——表观对流项与坐标变换项精确相消。
b：收缩以 1/R²(t) 的形式进入扩散项，是"收缩加速干燥"的数学来源；曲线由附件 2 的
R(t) 实测插值算出，右端标出全程放大倍数。
数据来源 problem_4_results.json 的 radius_profile；panel a 为结构示意。
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _figcommon import *  # noqa: F401,F403

FSG = 12.0
FSS = 11.0

d4 = load('problem_4_results')
rp = d4['radius_profile']
t_h = np.asarray(rp['times_s'], dtype=float) / 3600.0
R_cm = np.asarray(rp['R_cm'], dtype=float)
amp = (R_cm[0] / R_cm) ** 2

fig = plt.figure(figsize=(6.4, 3.0), layout='constrained')
gs = fig.add_gridspec(1, 2, width_ratios=[1.18, 1.0], wspace=0.10)
ax = fig.add_subplot(gs[0, 0])
axb = fig.add_subplot(gs[0, 1])

# ── a：固定域与收缩物理域的映射 ───────────────────────────────────────
bare(ax)
ax.set_xlim(0.0, 10.0)
ax.set_ylim(0.0, 9.4)

NN = 5
sv = np.linspace(0.0, 1.0, NN)
x0, span = 1.15, 7.60
y_fix = 8.05
ax.add_patch(Rectangle((x0, y_fix - 0.42), span, 0.84, facecolor=C_GRN1,
                       edgecolor=C_MEAN_T, lw=1.5, zorder=3))
ax.plot(x0 + span * sv, np.full(NN, y_fix), marker=MK[0], ms=5.5, ls='none',
        color=C_CENTER, zorder=6, markeredgecolor='white', markeredgewidth=1.0)
ax.text(x0 - 0.22, y_fix, '$s$', fontsize=FSG, color=C_TEXT, ha='right',
        va='center', zorder=7)
ax.text(x0 + span + 0.34, y_fix, '$1$', fontsize=FSS, color=C_TEXT, ha='left',
        va='center', zorder=7)

# 三个时刻的物理域：长度按实测 R(t) 的三个分位取，收缩比真实
frac = [1.0, float(R_cm[len(R_cm) // 3] / R_cm[0]),
        float(R_cm[-1] / R_cm[0])]
y_ph = [5.55, 3.55, 1.55]
cols = [C_RED1, C_RED2, C_HEAT]
for k, (fr, yy) in enumerate(zip(frac, y_ph)):
    w = span * fr
    ax.add_patch(Rectangle((x0, yy - 0.34), w, 0.68, facecolor=cols[k],
                           alpha=0.85 if k < 2 else 0.55, edgecolor=C_HEAT,
                           lw=1.3, zorder=3))
    ax.plot(x0 + w * sv, np.full(NN, yy), marker=MK[0], ms=4.6, ls='none',
            color=C_CENTER, zorder=6, markeredgecolor='white',
            markeredgewidth=0.8)
    # 条带长度本身就是 R(t)，标签只给时刻，短标签才留得下与条带边缘的净空
    ax.text(x0 + w + 0.34, yy, '$t_{%d}$' % k, fontsize=FSS, color=C_HEAT,
            ha='left', va='center', zorder=7)

# 物质点的轨迹：同一个 s 在三个时刻的物理位置连成线，直观显示节点随物质内移
# s=1 的轨迹线不画：它就是各条带的右边缘（已可见），画出来反而压住 R(t) 标签
for j in range(1, NN - 1):
    xs_tr = [x0 + span * sv[j]] + [x0 + span * fr * sv[j] for fr in frac]
    ys_tr = [y_fix] + y_ph
    ax.plot(xs_tr, ys_tr, color=C_AMB, lw=0.9, ls=':', zorder=4)
ax.text(x0 - 0.22, y_ph[1], '$r$', fontsize=FSG, color=C_TEXT, ha='right',
        va='center', zorder=7)
ax.text(x0 + 0.08, 0.52, '$r = s\\,R(t)$', fontsize=FSG, color=C_TEXT,
        ha='left', va='center', zorder=7)
panel(ax, 'a')

# ── b：1/R² 放大（实测 R(t) 算出） ────────────────────────────────────
axb.fill_between(t_h, 1.0, amp, color=C_ACC2, alpha=0.16, lw=0, zorder=1)
axb.plot(t_h, amp, color=C_ACC2, lw=2.1, zorder=4,
         label='$[R_{0}/R(t)]^{2}$ 实测插值')
axb.axhline(1.0, color=C_LINE, lw=1.2, ls=':', zorder=3)
axb.plot([], [], ls='none', label='终端放大 %.2f 倍' % amp[-1])
axb.set_xlabel('时间 $t$ (h)')
axb.set_ylabel('扩散项放大倍数')
axb.set_xlim(0.0, float(t_h[-1]))
axb.set_ylim(0.90, float(amp.max()) * 1.06)
# R(t) 前期收缩快 → 放大倍数在左上角就已抬起来（实测该处压住 45% 的曲线顶点），
# 图例改放右下：那里在曲线下方、只有 alpha=0.16 的淡填充，去框也不影响可读
axb.legend(loc='lower right', frameon=False, fontsize=FS_S,
           labelspacing=0.35, handlelength=1.5)
panel(axb, 'b')

out(fig, 'fig_mech_moving')
