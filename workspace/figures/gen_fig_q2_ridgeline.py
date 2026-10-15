"""问题 2 山脊图：6 个整时刻的径向剖面沿时间轴堆叠。

每条山脊为一个时刻的径向分布，纵向偏移即时间；左温度、右含水率。
山脊图把 6 条剖面在一个纵深方向排开，避免折线族互相压线。
数据来源 problem_2_results.json 的 ridgeline（0.5 h 间隔，21 个径向位置）。
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _figcommon import *  # noqa: F401,F403

d2 = load('problem_2_results')
rd = d2['ridgeline']
r = np.asarray(rd['radii_cm'], dtype=float)
th = np.asarray(rd['times_h'], dtype=float)
Tm = np.asarray(rd['T'], dtype=float)
Cm = np.asarray(rd['C'], dtype=float)
n = th.size

fig = plt.figure(figsize=(6.4, 3.4), layout='constrained')
gs = fig.add_gridspec(1, 2, wspace=0.10)
axT = fig.add_subplot(gs[0, 0])
axC = fig.add_subplot(gs[0, 1])


def ridge(ax, M, cols, step, base_lo):
    """按时刻自下而上堆叠：每条山脊做基线填充 + 描边。"""
    for k in range(n):
        off = k * step
        y = M[k] - base_lo + off
        ax.fill_between(r, off, y, color=cols[k], alpha=0.55, lw=0, zorder=k + 1)
        ax.plot(r, y, color=cols[k], lw=1.7, zorder=n + k + 1)
        ax.axhline(off, color=C_REF, lw=0.8, ls=':', zorder=1)
    ax.set_yticks([k * step for k in range(n)])
    ax.set_yticklabels(['%.1f h' % v for v in th])
    ax.set_xlim(0, R0_CM)
    ax.set_xlabel('径向位置 $r$ (cm)')


shT = fade(C_HEAT, n, lo=0.42, hi=1.0)
shC = fade(C_MOIST, n, lo=0.42, hi=1.0)
T_lo = float(Tm.min())
C_lo = float(Cm.min())
stepT, stepC = 7.0, 0.44
ridge(axT, Tm, shT, stepT, T_lo)
ridge(axC, Cm, shC, stepC, C_lo)

axT.set_ylabel('时刻（山脊基线，高度为 $T-%.1f$ ℃）' % T_lo)
axT.set_ylim(-1.6, (n - 1) * stepT + float(Tm.max()) - T_lo + 1.8)
panel(axT, 'a')

axC.set_ylabel('时刻（山脊基线，高度为 $C-%.3f$）' % C_lo)
axC.set_ylim(-0.11, (n - 1) * stepC + float(Cm.max()) - C_lo + 0.12)
panel(axC, 'b')

out(fig, 'fig_q2_ridgeline')
