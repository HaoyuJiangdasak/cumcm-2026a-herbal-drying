"""问题 3 灵敏度：13 个 ±10% 扰动算例对 t_dry 的相对效应（龙卷风图）+ 阈值不对称。

a：13 个扰动算例按物理量合并为 8 行，按 |Δt_dry| 排序的双向条形，右侧为延长、
左侧为缩短，右侧刻度给出 +10%/-10% 两侧的相对变化百分数。b：判据阈值 C_th 的
单侧扰动效应明显不对称，下调至 0.135 的延长幅度大于上调至 0.165 的缩短幅度。
数据来源 sensitivity_results.json 的 table / ranking_by_abs_effect /
threshold_asymmetry。
"""
import sys
import os
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from _figcommon import *  # noqa: F401,F403

ds = load('sensitivity_results')
base_h = float(ds['baseline_t_dry_h'])
rows = {r['case_id']: r for r in ds['table']}
rank = ds['ranking_by_abs_effect']
asym = ds['threshold_asymmetry']

# 把同一物理量的 +10% / -10% 合并成一条，取两侧效应画双向条
# 下标一律带花括号：印刷字号补偿只认 _{x} 形式的 mathtext，写 _x 会让 0.7 倍
# 的下标字形绕过补偿，在终检里留下 6.7pt 的字。
PAIRS = [('D', 'D_up', 'D_dn', '扩散系数 $D$'),
         ('hm', 'hm_up', 'hm_dn', '传质系数 $h_{m}$'),
         ('h', 'h_up', 'h_dn', '传热系数 $h$'),
         ('Ca', 'Ca_up', 'Ca_dn', '空气含湿 $C_{a}$'),
         ('Cth', 'Cth_up', 'Cth_dn', '判据阈值 $C_{th}$')]
SINGLE = [('latent', '计入蒸发潜热'),
          ('drive_B', '驱动力解释 B'),
          ('freeze_T', '冻结温度耦合')]

items = []
for _, up, dn, lab in PAIRS:
    items.append((lab, float(rows[up]['rel_change_pct']),
                  float(rows[dn]['rel_change_pct'])))
for cid, lab in SINGLE:
    items.append((lab, float(rows[cid]['rel_change_pct']), None))
items.sort(key=lambda it: max(abs(it[1]), abs(it[2] if it[2] is not None else 0.0)))

# 高度给到 4.3：条形行距太窄时，落盘期的字号还原会把 10pt 数值标注放大到
# ~13.5pt，行高 22.6px 顶破 22.5px 的行内净空，触发"标注越过所属行边界"。
fig = plt.figure(figsize=(6.4, 4.3), layout='constrained')
# 上边距要预留：落盘期的字号还原发生在版面解算之后，把 13pt 的面板标签放大到
# 20pt，默认 3pt 的边距会让标签顶出页面上沿 2.2pt。
fig.get_layout_engine().set(h_pad=0.11)
gs = fig.add_gridspec(2, 2, width_ratios=[1.62, 1.0], height_ratios=[1.0, 0.20],
                      wspace=0.10, hspace=0.04)
ax = fig.add_subplot(gs[0, 0])
axb = fig.add_subplot(gs[0, 1])
axleg = fig.add_subplot(gs[1, :])
axleg.set_axis_off()

y = np.arange(len(items))
hgt = 0.62
for i, (lab, vup, vdn) in enumerate(items):
    vals = [vup] if vdn is None else [vup, vdn]
    for v in vals:
        col = C_RED2 if v > 0 else C_GRN2
        edge = C_HEAT if v > 0 else C_MEAN_T
        ax.barh(i, v, height=hgt, color=col, edgecolor=edge, linewidth=0.9,
                zorder=3)
ax.axvline(0, color=C_LINE, lw=1.1, zorder=4)
ax.set_yticks(y)
ax.set_yticklabels([it[0] for it in items])
ax.set_xlabel('$t_{dry}$ 相对基准的变化 (%)')
ax.set_xlim(-31.0, 48.0)
ax.set_ylim(-0.62, len(items) - 0.38)
ax.grid(axis='y', visible=False)

# 数值走右侧刻度标签而非条内自由文字：落盘期的字号还原会把 10pt 标注放大到
# 15.4pt（68px），而行距只有 79px，条内文字必然顶破行界；刻度标签不受此限。
axv = ax.twinx()
axv.set_ylim(ax.get_ylim())
axv.set_yticks(y)
axv.set_yticklabels(['%+.1f' % vup if vdn is None else '%+.1f/%+.1f' % (vup, vdn)
                     for _, vup, vdn in items])
axv.tick_params(axis='y', length=0, pad=2.0)
axv.grid(False)
for _s in axv.spines.values():
    _s.set_visible(False)
panel(ax, 'a')

# b：阈值单侧效应的不对称
lo = float(rows['Cth_dn']['stop_threshold'])
hi = float(rows['Cth_up']['stop_threshold'])
xs = np.array([lo, C_TH, hi])
ts = np.array([float(rows['Cth_dn']['t_dry_h']), base_h,
               float(rows['Cth_up']['t_dry_h'])])
lcurve, = axb.plot(xs, ts, color=C_MOIST, lw=2.0, marker=MK[0], ms=6.5,
                   zorder=4, markerfacecolor='white', markeredgewidth=1.5,
                   label='实算 $t_{dry}(C_{th})$')
lbase, = axb.plot([C_TH], [base_h], marker=MK[1], ms=8.5, color=C_ACC2,
                  zorder=5, ls='none', label='基准 %.4f h' % base_h)
# 三个 t_dry 值一律走图例标签：图内自由文字会被落盘期的自动排版挪到曲线上
lends, = axb.plot([], [], ls='none', label='端点 %.2f / %.2f h' % (ts[0], ts[2]))
axb.set_xlabel('判据阈值 $C_{th}$ (kg/kg)')
axb.set_ylabel('$t_{dry}$ (h)')
axb.set_xlim(0.1305, 0.1695)
axb.set_ylim(41.0, 84.0)
# 只留两端刻度：窄格子里三个五位数刻度在落盘期放大后会互相重叠，中间的基准
# 阈值已由星标与图例给出。
axb.set_xticks([0.135, 0.165])
panel(axb, 'b')

hs = [Rectangle((0, 0), 1, 1, facecolor=C_RED2, edgecolor=C_HEAT, lw=0.9),
      Rectangle((0, 0), 1, 1, facecolor=C_GRN2, edgecolor=C_MEAN_T, lw=0.9),
      lcurve, lbase, lends]
axleg.legend(hs, ['延长 $t_{dry}$', '缩短 $t_{dry}$', '$t_{dry}(C_{th})$ 实算',
                  '基准 %.2f h' % base_h, '端点 %.2f / %.2f h' % (ts[0], ts[2])],
             loc='center', ncol=3, frameon=False, fontsize=FS_S,
             handlelength=1.2, columnspacing=1.4, handletextpad=0.5,
             borderpad=0.1)

out(fig, 'fig_sensitivity')
