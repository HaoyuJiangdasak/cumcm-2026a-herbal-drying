"""四模型对比：平均含水率轨迹 + 能力矩阵热图。

a：M1—M4 的 C̄(t) 实算轨迹，灰带为题面 S9 给出的 2—3 d 外部锚点窗口，
只有 M3、M4 的终点落在窗内；圆点为各模型自身判据下的 t_dry。
b：能力矩阵，逐格标记模型是否具备该能力（不适用者留"—"）。
数据来源 model_compare_results.json 的 series_by_model / table / biot_evidence。
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _figcommon import *  # noqa: F401,F403

dm = load('model_compare_results')
rows = {r['model']: r for r in dm['table']}
sm = dm['series_by_model']
win = [float(v) * 24.0 for v in dm['S9_window_days']]

ORDER = ['M1', 'M2', 'M3', 'M4']
COLS = {'M1': C_AMB, 'M2': C_ACC1, 'M3': C_MOIST, 'M4': C_ACC2}
STY = {'M1': ':', 'M2': '-.', 'M3': '-', 'M4': '--'}

fig = plt.figure(figsize=(6.4, 4.05), layout='constrained')
gs = fig.add_gridspec(2, 1, height_ratios=[1.0, 0.86], hspace=0.06)
ax = fig.add_subplot(gs[0, 0])
axb = fig.add_subplot(gs[1, 0])

ax.axvspan(win[0], win[1], color=C_GRN1, alpha=0.55, lw=0, zorder=0)
for m in ORDER:
    t = np.asarray(sm[m]['times_h'], dtype=float)
    cb = np.asarray(sm[m]['C_bar'], dtype=float)
    ax.plot(t, cb, color=COLS[m], lw=2.0, ls=STY[m], zorder=3,
            label='%s  %.2f h' % (m, float(rows[m]['t_dry_h'])))
    td = float(rows[m]['t_dry_h'])
    ax.plot([td], [float(np.interp(td, t, cb))], marker=MK[0], ms=6.0,
            color=COLS[m], markerfacecolor='white', markeredgewidth=1.4,
            ls='none', zorder=5)
ax.axhline(C_TH, color=C_LINE, lw=1.2, ls=':', zorder=2)
ax.set_xlabel('时间 $t$ (h)')
ax.set_ylabel('平均含水率 $\\bar{C}$ (kg/kg)')
ax.set_xlim(0, 82.0)
ax.set_ylim(0, 2.72)
ax.legend(loc='upper right', ncol=2, frameon=False, fontsize=FS_S,
          labelspacing=0.35, handlelength=1.7, columnspacing=0.9)
panel(ax, 'a')

# b：能力矩阵。1=具备，0=不具备，nan=不适用
# 首行标签必须与其余各行同为 6 个全角字：混排空格与 "S9" 会多占约 5 px，
# 触发 save 时的 y 标签溢出截断（实测被截成「落在 S9 锚…」）。锚点窗口出处见题注。
CAPS = [('落在锚点窗内', [float(rows[m]['in_S9_window']) for m in ORDER]),
        ('可给径向分布', [float(rows[m]['radial_capability']) for m in ORDER]),
        ('可用逐点判据', [float(rows[m]['can_deliver_pointwise']) for m in ORDER]),
        ('物性随场重算', [0.0, 0.0, 1.0, 1.0]),
        ('可给收缩半径', [0.0, 0.0, 0.0, 1.0])]
M = np.array([c[1] for c in CAPS], dtype=float)
cmap = LinearSegmentedColormap.from_list('cap', [NP['red_1'], NP['green_2']])
axb.imshow(M, cmap=cmap, vmin=0, vmax=1, aspect='auto', origin='upper')
for i in range(M.shape[0]):
    for j in range(M.shape[1]):
        axb.text(j, i, '有' if M[i, j] > 0.5 else '无', fontsize=FS_S,
                 color=C_TEXT, ha='center', va='center')
axb.set_xticks(range(4))
axb.set_xticklabels(['%s\n%s' % (m, rows[m]['dim']) for m in ORDER])
axb.set_yticks(range(len(CAPS)))
axb.set_yticklabels([c[0] for c in CAPS])
axb.set_xticks(np.arange(-0.5, 4, 1), minor=True)
axb.set_yticks(np.arange(-0.5, len(CAPS), 1), minor=True)
axb.grid(which='minor', color='white', linewidth=1.6)
axb.grid(which='major', visible=False)
axb.tick_params(which='minor', length=0)
panel(axb, 'b')

out(fig, 'fig_model_compare')
