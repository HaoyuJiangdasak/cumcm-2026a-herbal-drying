"""数值验证三联：网格加倍相对变化、时间步一阶收敛、与 Bessel 解析解的偏差。

a：三个问题各自的关键输出在网格加倍下的相对变化，全部落在 0.1% 容差线以下；
b：固定网格下 T 中心随 Δt 收敛，斜率拟合出的观测阶约为 1，与后向 Euler 的
理论一阶一致，虚线为 Richardson 外推极限；c：常物性算例与十项 Bessel 级数
解析解的中心温度绝对偏差随 Δt 单调下降。
数据来源 verification_results.json 的 grid_convergence / time_convergence /
analytic_cross_check / mass_conservation。
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _figcommon import *  # noqa: F401,F403

dv = load('verification_results')
gc = dv['grid_convergence']
tc = dv['time_convergence']['problem1_first_order']
ac = dv['analytic_cross_check']
mc = dv['mass_conservation']
tol = float(gc['tol_pct'])

fig = plt.figure(figsize=(6.4, 3.15), layout='constrained')
gs = fig.add_gridspec(2, 3, height_ratios=[1.0, 0.30], wspace=0.10, hspace=0.04)
ax = fig.add_subplot(gs[0, 0])
axb = fig.add_subplot(gs[0, 1])
axc = fig.add_subplot(gs[0, 2])
legax = [fig.add_subplot(gs[1, k]) for k in range(3)]
for _a in legax:
    _a.set_axis_off()

# a：网格加倍的相对变化
SETS = [('problem1', '问题 1', C_HEAT, MK[0]),
        ('problem3', '问题 3', C_MOIST, MK[1]),
        ('problem4', '问题 4', C_ACC2, MK[2])]
for key, lab, col, mk in SETS:
    ref = gc[key]['refinement']
    xs = np.array([r['to'] for r in ref], dtype=float)
    ys = np.array([abs(r['delta_pct']) for r in ref], dtype=float)
    ax.plot(xs, ys, color=col, lw=1.9, marker=mk, ms=6.0, zorder=4,
            markerfacecolor='white', markeredgewidth=1.3, label=lab)
ax.axhline(tol, color=C_LINE, lw=1.3, ls=':', zorder=3,
           label='容差 %.1f%%' % tol)
ax.set_xscale('log', base=2)
ax.set_yscale('log')
ax.set_xlabel('加倍后的网格数 $N$')
ax.set_ylabel('关键输出相对变化 (%)')
ax.set_xticks([40, 80, 160])
ax.xaxis.set_major_formatter(FuncFormatter(lambda v, _: '%d' % v))
ax.set_ylim(4e-3, 0.6)
hA, lA = ax.get_legend_handles_labels()
legax[0].legend(hA, lA, loc='center', ncol=2, frameon=False, fontsize=FS_S,
                handlelength=1.3, columnspacing=0.9, labelspacing=0.3,
                handletextpad=0.45, borderpad=0.1)
panel(ax, 'a')

# b：时间步一阶收敛
dts = np.asarray(tc['dt_s'], dtype=float)
Tcs = np.asarray(tc['T_center'], dtype=float)
rich = float(tc['richardson_limit'])
order = float(tc['observed_order_halving'])
o = np.argsort(dts)
axb.plot(dts[o], Tcs[o], color=C_HEAT, lw=2.0, marker=MK[0], ms=6.0,
         zorder=4, markerfacecolor='white', markeredgewidth=1.3,
         label='实算 $T$ 中心')
axb.axhline(rich, color=C_ACC1, lw=1.4, ls='--', zorder=3,
            label='Richardson %.4f ℃' % rich)
axb.plot([], [], ls='none', label='观测阶 %.3f' % order)
axb.set_xscale('log')
axb.set_xlabel('时间步 $\\Delta t$ (s)')
axb.set_ylabel('中心温度 $T(0)$ (℃)')
axb.set_ylim(33.5748, 33.5800)
axb.yaxis.set_major_formatter(FuncFormatter(lambda v, _: '%.4f' % v))
hB, lB = axb.get_legend_handles_labels()
legax[1].legend(hB, lB, loc='center', ncol=1, frameon=False, fontsize=FS_S,
                handlelength=1.3, labelspacing=0.25, handletextpad=0.45,
                borderpad=0.1)
panel(axb, 'b')

# c：与解析解的偏差
cs = ac['cases']
xd = np.array([float(c['dt_s']) for c in cs])
ec = np.array([float(c['err_center']) for c in cs])
o2 = np.argsort(xd)
axc.fill_between(xd[o2], 1e-5, ec[o2], color=C_MOIST, alpha=0.16, lw=0,
                 zorder=1)
axc.plot(xd[o2], ec[o2], color=C_MOIST, lw=2.0, marker=MK[1], ms=6.0,
         zorder=4, markerfacecolor='white', markeredgewidth=1.3,
         label='$|T_{num}-T_{ana}|$')
axc.plot([], [], ls='none', label='守恒残差 ≤%.1e' % float(mc['eps_max']))
axc.set_xscale('log')
axc.set_yscale('log')
axc.set_xlabel('时间步 $\\Delta t$ (s)')
axc.set_ylabel('中心温度绝对偏差 (℃)')
axc.set_ylim(1.2e-4, 1.4e-2)
hC, lC = axc.get_legend_handles_labels()
legax[2].legend(hC, lC, loc='center', ncol=1, frameon=False, fontsize=FS_S,
                handlelength=1.3, labelspacing=0.25, handletextpad=0.45,
                borderpad=0.1)
panel(axc, 'c')

out(fig, 'fig_verification')
