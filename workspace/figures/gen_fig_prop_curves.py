"""附录 2/3/4 物性对照：ρ、c_p、λ 随含水率变化 + 三附录在 C→0 的隔离核对。

四个面板逐一比较三套附录的密度、比热容、导热系数与扩散系数；附录 2 为常物性，
附录 3/4 随 C（c_p、λ 走湿基 C/(C+1)）变化。全部由 code/utils.Props 现算。
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _figcommon import *  # noqa: F401,F403
import utils as U

C = np.linspace(1e-6, C_INIT, 400)
T_ref = 50.0
P2, P3, P4 = U.Props(2), U.Props(3), U.Props(4)
cols = {2: C_AMB, 3: C_MOIST, 4: C_ACC2}
lss = {2: ':', 3: '-', 4: '--'}
names = {2: '附录 2（问题 1）', 3: '附录 3（问题 2、3）', 4: '附录 4（问题 4）'}

fig = plt.figure(figsize=(6.4, 4.6), layout='constrained')
gs = fig.add_gridspec(2, 2, hspace=0.06, wspace=0.06)
axs = [fig.add_subplot(gs[i // 2, i % 2]) for i in range(4)]

specs = [('rho', '密度 $\\rho$ (kg/m$^3$)', 1.0, 'a'),
         ('cp', '比热容 $c_p$ (J/(kg$\\cdot$K))', 1.0, 'b'),
         ('lam', '导热系数 $\\lambda$ (W/(m$\\cdot$K))', 1.0, 'c')]

for ax, (attr, ylab, sc, tag) in zip(axs[:3], specs):
    for k, pr in ((2, P2), (3, P3), (4, P4)):
        y = np.asarray(getattr(pr, attr)(C), dtype=float) * sc
        if y.size == 1:
            y = np.full_like(C, float(y))
        ax.plot(C, y, color=cols[k], lw=2.0, ls=lss[k], label=names[k])
    ax.set_ylabel(ylab)
    ax.set_xlim(0, C_INIT)
    ax.axvline(C_TH, color=C_SURF, lw=1.2, ls='-.', zorder=1)
    panel(ax, tag, dx=-0.24, dy=1.09)

axs[0].tick_params(labelbottom=False)
axs[1].tick_params(labelbottom=False)
axs[2].set_xlabel('干基含水率 $C$ (kg/kg)')
axs[0].annotate('$C_{th}$=0.15', xy=(C_TH, 700), xytext=(0.42, 690),
                fontsize=FS_S, color=C_TEXT, va='center',
                arrowprops=dict(arrowstyle='-|>', color=C_LINE, lw=1.0,
                                shrinkA=2, shrinkB=2))

axd = axs[3]
T_arr = np.full_like(C, T_ref)
for k, pr in ((2, P2), (3, P3), (4, P4)):
    y = np.asarray(pr.D(C, T_arr), dtype=float)
    axd.semilogy(C, y, color=cols[k], lw=2.0, ls=lss[k], label=names[k])
axd.axvline(C_TH, color=C_SURF, lw=1.2, ls='-.', zorder=1)
axd.set_xlabel('干基含水率 $C$ (kg/kg)')
axd.set_ylabel('扩散系数 $D$ (m$^2$/s)')
axd.set_xlim(0, C_INIT)
axd.set_ylim(1e-13, 3e-8)
axd.yaxis.set_major_locator(LogLocator(base=10.0, numticks=6))
panel(axd, 'd', dx=-0.24, dy=1.09)

d3 = float(np.asarray(P3.D(np.array([C_TH]), np.array([T_ref])))[0])
d4 = float(np.asarray(P4.D(np.array([C_TH]), np.array([T_ref])))[0])
axd.annotate('$C_{th}$ 处 $D$\n附录 3 %.2e\n附录 4 %.2e' % (d3, d4),
             xy=(C_TH, d4), xytext=(0.66, 1.2e-12), fontsize=FS_S,
             color=C_TEXT, va='center',
             arrowprops=dict(arrowstyle='-|>', color=C_LINE, lw=1.0,
                             shrinkA=2, shrinkB=2))

handles, labels = axs[0].get_legend_handles_labels()
fig.legend(handles, labels, loc='outside lower center', ncol=3, frameon=False,
           fontsize=FS_S)

out(fig, 'fig_prop_curves')
