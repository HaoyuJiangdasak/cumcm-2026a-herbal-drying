"""问题 1 径向剖面族：7 个输出时刻的 T(r) 与 C(r)，同族深浅渐变区分时刻。

左：温度剖面由外向内推进；右：水分剖面仅表层下降，中心 30 min 内几乎未动。
数据取 problem_1_results.json 的 profiles（21 个径向位置，0.1 cm 间隔）。
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _figcommon import *  # noqa: F401,F403

d1 = load('problem_1_results')
pf = d1['profiles']
r = np.asarray(pf['radii_cm'], dtype=float)
ts = np.asarray(pf['times_s'], dtype=float)
Tm = np.asarray(pf['T'], dtype=float)
Cm = np.asarray(pf['C'], dtype=float)
fin = d1['final_state']

fig = plt.figure(figsize=(6.4, 3.0), layout='constrained')
gs = fig.add_gridspec(1, 2, wspace=0.06)
axT = fig.add_subplot(gs[0, 0])
axC = fig.add_subplot(gs[0, 1])

shT = fade(C_HEAT, ts.size, lo=0.32, hi=1.0)
shC = fade(C_MOIST, ts.size, lo=0.32, hi=1.0)
for k in range(ts.size):
    axT.plot(r, Tm[k], color=shT[k], lw=1.7, marker=MK[k % len(MK)], ms=3.4,
             markevery=4, label='%d min' % round(ts[k] / 60.0))
    axC.plot(r, Cm[k], color=shC[k], lw=1.7, marker=MK[k % len(MK)], ms=3.4,
             markevery=4)

axT.axhline(T_INIT, color=C_REF, lw=1.1, ls=':', zorder=1)
axT.set_xlabel('径向位置 $r$ (cm)')
axT.set_ylabel('温度 $T$ (℃)')
axT.set_xlim(0, R0_CM)
axT.set_ylim(27.4, 39.2)
axT.legend(loc='upper left', ncol=2, frameon=False, fontsize=FS_S,
           labelspacing=0.3, handlelength=1.5, columnspacing=0.9)
axT.annotate('%.3f' % fin['T_surface'], xy=(R0_CM, fin['T_surface']),
             xytext=(1.66, 38.5), fontsize=FS_S, color=C_TEXT, ha='center',
             arrowprops=dict(arrowstyle='-|>', color=C_LINE, lw=1.0,
                             shrinkA=2, shrinkB=2))
panel(axT, 'a')

axC.axhline(C_INIT, color=C_REF, lw=1.1, ls=':', zorder=1)
axC.set_xlabel('径向位置 $r$ (cm)')
axC.set_ylabel('干基含水率 $C$ (kg/kg)')
axC.set_xlim(0, R0_CM)
axC.set_ylim(1.42, 2.62)
axC.annotate('%.4f' % fin['C_surface'], xy=(R0_CM, fin['C_surface']),
             xytext=(1.30, 1.78), fontsize=FS_S, color=C_TEXT, ha='center',
             arrowprops=dict(arrowstyle='-|>', color=C_LINE, lw=1.0,
                             shrinkA=2, shrinkB=2))
axC.text(0.62, 2.53, '中心 $C$ 未动', fontsize=FS_S, color=C_MOIST, va='center')
panel(axC, 'b')

out(fig, 'fig_q1_profiles')
