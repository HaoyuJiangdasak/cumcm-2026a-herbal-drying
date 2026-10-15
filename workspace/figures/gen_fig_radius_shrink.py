"""附件 2 半径收缩：R(t) 渐变填充折线 + 收缩速率与截面积占比。

主图为 145 点实测半径与 PCHIP 保形插值，填充带表示药材实体截面；插图给出
收缩速率 |dR/dt| 与面积比 (R/R0)^2。数据来源 user_data/附件2.xlsx。
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _figcommon import *  # noqa: F401,F403
import utils as U

t_s, R_cm = (np.asarray(x, dtype=float) for x in U.load_attachment2())
t_h = t_s / 3600.0
rad, rad_dot, _ = U.make_radius()
tf = np.linspace(t_s[0], t_s[-1], 2000)
Rf = np.asarray([rad(tt) for tt in tf]) * 100.0
tfh = tf / 3600.0

R_end = float(R_cm[-1])
ratio_end = R_end / R0_CM
i_last = int(np.argmax(R_cm <= R_end + 1e-12))
t_flat_h = t_h[i_last]

fig = plt.figure(figsize=(6.4, 4.3), layout='constrained')
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 0.72], width_ratios=[1.0, 1.0],
                      hspace=0.05, wspace=0.22)
ax = fig.add_subplot(gs[0, :])
axr = fig.add_subplot(gs[1, 0])
axa = fig.add_subplot(gs[1, 1])

ax.fill_between(tfh, 0, Rf, color=C_MOIST, alpha=0.16, lw=0, zorder=1)
ax.fill_between(tfh, Rf, R0_CM, color=C_REF, alpha=0.5, lw=0, zorder=1)
ax.plot(tfh, Rf, color=C_MOIST, lw=2.2, zorder=3, label='PCHIP 保形插值')
ax.scatter(t_h, R_cm, s=16, color=C_MOIST, marker=MK[2], zorder=4,
           edgecolor='white', linewidth=0.5, label='附件 2 实测 145 点')
ax.axhline(R0_CM, color=C_LINE, lw=1.1, ls=':', zorder=2)
ax.axhline(R_end, color=C_SURF, lw=1.4, ls='--', zorder=2)
ax.set_ylabel('药材半径 $R$ (cm)')
ax.set_xlim(0, t_h[-1])
ax.set_ylim(0.0, 2.35)
ax.set_yticks([0.0, 0.5, 1.0, 1.5, 2.0])
ax.tick_params(labelbottom=False)
ax.legend(loc='upper right', frameon=False, fontsize=FS_S,
          labelspacing=0.35, handlelength=1.6)
ax.annotate('$R$=%.3f cm\n$R/R_0$=%.3f' % (R_end, ratio_end),
            xy=(t_flat_h, R_end), xytext=(t_flat_h - 21.0, 0.72), fontsize=FS_S,
            color=C_TEXT, arrowprops=dict(arrowstyle='-|>', color=C_LINE, lw=1.1,
                                          shrinkA=2, shrinkB=3))
ax.text(3.0, 1.55, '实体截面', fontsize=FS_S, color=C_MOIST, va='center')
panel(ax, 'a', dx=-0.075, dy=1.07)

rdot = np.asarray([-rad_dot(tt) for tt in tf]) * 100.0 * 1e4  # cm/s → 10^-4 cm/s
axr.fill_between(tfh, 0, rdot, color=C_ACC1, alpha=0.28, lw=0)
axr.plot(tfh, rdot, color=C_ACC1, lw=1.8)
axr.axhline(0, color=C_REF, lw=1.0)
axr.set_xlabel('时间 $t$ (h)')
axr.set_ylabel('$|\\mathrm{d}R/\\mathrm{d}t|$\n($10^{-4}$ cm/s)')
axr.set_xlim(0, t_h[-1])
panel(axr, 'b', dx=-0.30, dy=1.20)

area = (Rf / R0_CM) ** 2
axa.fill_between(tfh, 0, area, color=C_ACC2, alpha=0.24, lw=0)
axa.plot(tfh, area, color=C_ACC2, lw=1.8)
axa.axhline(ratio_end ** 2, color=C_SURF, lw=1.3, ls='--')
axa.set_xlabel('时间 $t$ (h)')
axa.set_ylabel('截面积比 $(R/R_0)^2$')
axa.set_xlim(0, t_h[-1])
axa.set_ylim(0, 1.14)
axa.annotate('%.3f' % (ratio_end ** 2), xy=(t_h[-1] * 0.62, ratio_end ** 2),
             xytext=(t_h[-1] * 0.40, 0.14), fontsize=FS_S, color=C_TEXT,
             arrowprops=dict(arrowstyle='-|>', color=C_LINE, lw=1.0,
                             shrinkA=2, shrinkB=2))
panel(axa, 'c', dx=-0.26, dy=1.20)

out(fig, 'fig_radius_shrink')
