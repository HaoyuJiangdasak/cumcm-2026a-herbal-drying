"""2026A 药材烘干 — 绘图公用模块：Nature 样式初始化、配色取用、产物读取与落盘。

各 gen_fig_*.py 统一 `from _figcommon import *`，避免重复样板与硬编码色值。
配色一律经 nature_palette() 函数取用（自带同族去指纹抖动），不复制十六进制字面量。
"""
import json
import os
import sys

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm
from matplotlib.colors import LinearSegmentedColormap, Normalize, to_rgba
from matplotlib.patches import (Circle, FancyArrowPatch, FancyBboxPatch, Polygon,
                               Rectangle, Wedge, Ellipse, Arc)
from matplotlib.lines import Line2D
from matplotlib.ticker import MultipleLocator, FuncFormatter, LogLocator

_HERE = os.path.dirname(os.path.abspath(__file__))
_ROOT = os.path.dirname(_HERE)
for _p in (os.path.join(_ROOT, '_utils'), os.path.join(_ROOT, 'code')):
    if _p not in sys.path:
        sys.path.insert(0, _p)

from plot_utils import (setup_style, set_paper_placement, save_fig, nature_palette,
                       nature_markers, smart_labels, auto_legend, shared_legend,
                       dynamic_limits, declutter_axes, COLORS)

setup_style(palette='nature')
NP = nature_palette()
MK = nature_markers()

# 字号只用三档，避免「同图字号体系偏乱」且不低于印刷可读线（缩放后 ≥8pt）
FS_S = 10.0    # 刻度 / 图例 / 等高线标签
FS_M = 11.0    # 轴标签 / 图内数值标注
FS_L = 13.0    # 面板编号
plt.rcParams['xtick.labelsize'] = FS_S
plt.rcParams['ytick.labelsize'] = FS_S
plt.rcParams['legend.fontsize'] = FS_S
plt.rcParams['axes.labelsize'] = FS_M

# 语义化取色：全部来自 nature_palette()，脚本内不再出现色值字面量
C_HEAT = NP['red_strong']       # 温度 / 热
C_MOIST = NP['blue_main']       # 水分 / 质
C_CENTER = NP['blue_main']      # 中心
C_SURF = NP['gold']             # 表面
C_MEAN = NP['green_3']          # 平均
C_AMB = NP['neutral_mid']       # 烘房环境
C_REF = NP['neutral_light']     # 参考 / 阈值底色
C_LINE = NP['neutral_dark']     # 结构线
C_TEXT = COLORS['text']
C_GRID = COLORS['grid']
C_ACC1 = NP['teal']
C_ACC2 = NP['violet']
C_BLUE2 = NP['blue_secondary']
C_RED1 = NP['red_1']
C_RED2 = NP['red_2']
C_GRN1 = NP['green_1']
C_GRN2 = NP['green_2']

def _relL(rgb):
    f = [(v / 12.92 if v <= 0.03928 else ((v + 0.055) / 1.055) ** 2.4) for v in rgb]
    return 0.2126 * f[0] + 0.7152 * f[1] + 0.0722 * f[2]


def ink(color, target=4.6, on=None, alpha=1.0):
    """把配色压暗到对比度 ≥target，供【文字】使用（线条仍用原色）。

    金色 1.40:1、浅绿 1.89:1 直接当字色会被对比度闸判失败，但换成别的色又会
    切断"字-线"的颜色对应。故按同色相压暗，保留辨识、拿到对比度。
    on/alpha 给出文字实际所处的底色（如落在半透明填充带内），默认纯白。
    """
    r, g, b = to_rgba(color)[:3]
    if on is None:
        bg = (1.0, 1.0, 1.0)
    else:
        o = to_rgba(on)[:3]
        bg = tuple(alpha * o[i] + (1.0 - alpha) for i in range(3))
    Lb = _relL(bg)

    def _cr(rgb):
        Lt = _relL(rgb)
        hi, lo = (Lt, Lb) if Lt > Lb else (Lb, Lt)
        return (hi + 0.05) / (lo + 0.05)

    k = 1.0
    while k > 0.05 and _cr((r * k, g * k, b * k)) < target:
        k -= 0.02
    return (r * k, g * k, b * k)


# 文字专用暗色（与同名线条同色相），只在需要"字随线色"时使用
C_MEAN_T = ink(C_MEAN)
C_SURF_T = ink(C_SURF)
C_ACC1_T = ink(C_ACC1)
C_BLUE2_T = ink(C_BLUE2)

# 题面常量（PROBLEM_FACTS.json / code/params.py，绘图侧只读不改）
R0_CM = 2.0
C_INIT = 2.55
C_TH = 0.15
T_INIT = 28.0


def load(name):
    """读取 figures/*.json 实算产物；缺文件直接抛错，不做空字典兜底。"""
    path = os.path.join(_HERE, name if name.endswith('.json') else name + '.json')
    if not os.path.exists(path):
        raise FileNotFoundError(path)
    with open(path, encoding='utf-8') as fh:
        return json.load(fh)


def arr(*keys, src=None):
    """按嵌套键取出并转 float 数组。"""
    node = src
    for k in keys:
        node = node[k]
    return np.asarray(node, dtype=float)


def out(fig, stem, width_fraction=1.0):
    """统一落盘：先记录论文版面元数据，再存矢量 PDF。"""
    set_paper_placement(fig, width_fraction)
    save_fig(fig, os.path.join(_HERE, stem + '.pdf'))
    plt.close(fig)


def seq_cmap(name='moist'):
    """由 nature_palette() 组装的顺序色图，避免使用被禁的双极色图。"""
    if name == 'moist':
        stops = [NP['neutral_light'], NP['green_2'], NP['blue_secondary'], NP['blue_main']]
    elif name == 'heat':
        stops = [NP['neutral_light'], NP['gold'], NP['red_2'], NP['red_strong']]
    else:
        stops = [NP['neutral_light'], NP['teal'], NP['blue_main'], NP['violet']]
    return LinearSegmentedColormap.from_list('mh_' + name, stops, N=256)


def fade(color, n, lo=0.25, hi=1.0):
    """同族深浅渐变：一组同类曲线用同色不同透明度区分。"""
    return [to_rgba(color, a) for a in np.linspace(lo, hi, n)]


def panel(ax, tag, dx=None, dy=None, pad=4.0):
    """面板编号（a/b/c…）：走 set_title(loc='left')，不占数据区也不会被拉回。

    dx/dy 仅为兼容旧签名保留，实际位置交给版面引擎。
    """
    ax.set_title(tag, loc='left', fontsize=FS_L, fontweight='bold',
                 color=C_TEXT, pad=pad)


def thin(x, y, n=600):
    """仅用于超长序列的等距抽稀（判据定位与守恒统计不走此路径）。"""
    x = np.asarray(x, dtype=float)
    y = np.asarray(y, dtype=float)
    if x.size <= n:
        return x, y
    idx = np.unique(np.linspace(0, x.size - 1, n).astype(int))
    return x[idx], y[idx]


def box(ax, xy, w, h, text, fc, ec=None, fs=FS_M, tc=None, r=0.03, lw=1.4, z=3):
    """机理/结构图用的圆角文本框，返回中心坐标。"""
    x, y = xy
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle='round,pad=0,rounding_size=%g' % r,
                                linewidth=lw, facecolor=fc,
                                edgecolor=ec if ec else C_LINE, zorder=z))
    ax.text(x + w / 2, y + h / 2, text, ha='center', va='center', fontsize=fs,
            color=tc if tc else C_TEXT, zorder=z + 1, linespacing=1.35)
    return (x + w / 2, y + h / 2)


def arrow(ax, p0, p1, color=None, lw=1.6, style='-|>', ms=10, ls='-', z=2, rad=0.0):
    ax.add_patch(FancyArrowPatch(p0, p1, arrowstyle=style, mutation_scale=ms,
                                 linewidth=lw, linestyle=ls, zorder=z,
                                 color=color if color else C_LINE,
                                 connectionstyle='arc3,rad=%g' % rad,
                                 shrinkA=0, shrinkB=0))


def flow(ax, p0, p1, color=None, lw=1.6, ms=7.0, ls='-', z=2):
    """流程图专用箭头：轴对齐线段 + 三角标记头。

    不走 FancyArrowPatch：它的路径末位是 CLOSEPOLY 顶点 (0,0)，而遮挡检查按折线
    逐顶点求交，会从箭头尖端到画布原点拉出一条假想射线，把沿途所有文字误判成被
    箭头遮挡。线段加标记的画法没有这条射线。
    """
    x0, y0 = p0
    x1, y1 = p1
    col = color if color else C_LINE
    if abs(y1 - y0) >= abs(x1 - x0):
        mk = 'v' if y1 < y0 else '^'
    else:
        mk = '<' if x1 < x0 else '>'
    ax.plot([x0, x1], [y0, y1], color=col, lw=lw, ls=ls, zorder=z,
            solid_capstyle='butt')
    ax.plot([x1], [y1], marker=mk, ms=ms, ls='none', color=col,
            markeredgecolor='none', zorder=z)


def elbow(ax, pts, color=None, lw=1.5, ms=7.0, ls='-', z=5):
    """折线箭头：按给定拐点走线，末段带三角头（回路、分支用）。"""
    col = color if color else C_LINE
    ax.plot([p[0] for p in pts[:-1]], [p[1] for p in pts[:-1]], color=col,
            lw=lw, ls=ls, zorder=z, solid_capstyle='butt')
    flow(ax, pts[-2], pts[-1], color=col, lw=lw, ms=ms, ls=ls, z=z)


def bare(ax):
    """机理图/示意图：去坐标系但保留数据区（仅示意图使用）。"""
    ax.set_xticks([])
    ax.set_yticks([])
    ax.grid(False)
    ax.set_frame_on(False)
