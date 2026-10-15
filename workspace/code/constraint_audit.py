# -*- coding: utf-8 -*-
"""约束复核 K01—K16：只从已落盘的产物重新算，不读求解器的任何通过标记。

设计铁律（与 _utils 的各道闸同一套哲学）：
  1. 每条约束都自己把数值算一遍，再和限值比。凡是 JSON 里写着 "status": "PASS"、
     "verdict"、"ok" 之类的字段，本文件一律不读——只读原始场、原始序列、原始标量。
  2. 数据来源优先级：xlsx 交付件 > result_grid 全分辨率场 > criterion_series 抽样序列。
     每条约束都登记自己用的是哪一档，抽样档不冒充全量档。
  3. 基准/对照算例按关键词自动发现，发现了但没被任何一条约束覆盖 → 直接 FAIL。
     这样以后新增分支（新的基准、新的对照）不补审计就过不了闸。
  4. K05 按坐标系分档：物质坐标必须守恒，实验室坐标必须违反（报告第 1.2 节的推翻论证
     靠的就是它违反）。「按设计违反」是显式判定，不是放行。
  5. 铁律 1 由运行期护栏保证，不靠正则自查：load() 把产物包成 NoFlags 视图，取用
     capability_verdict / status / verdict / ok 一类字段直接抛 AssertionError。
     写错了会在本文件里当场炸，而不是靠事后扫描猜。

运行：python _utils/run_compute.py code/constraint_audit.py
退出码：0=全部通过 1=存在真实违规
"""
from __future__ import annotations

import json
import re
import sys
from pathlib import Path

import numpy as np
import openpyxl

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_utils"))

from params import *  # noqa: F401,F403,E402
from utils import C_FLOOR, Props, ensure_dirs, jsonable, save_json  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
BASELINE_RE = re.compile(r"基准|baseline|参考态|对照|control|主模型|reference_model")
SELF_PRODUCT = "constraint_audit_results.json"     # 本文件自己的产物，不参与基准棘轮
FLAG_KEYS = frozenset({"capability_verdict", "status", "verdict", "ok", "passed", "all_pass"})
_CACHE: dict[str, "NoFlags"] = {}
TOUCHED: set[str] = set()          # 本次审计实际检查过的实体（file 或 file::branch）
GUARD_HITS: list[str] = []         # 护栏实际拦下的取用（自证探针会往里写）


class NoFlags(dict):
    """产物的只读视图：一取通过标记就抛。铁律 1 由此在运行期成立，而非靠正则事后猜。"""

    def __getitem__(self, key):
        if key in FLAG_KEYS:
            GUARD_HITS.append(str(key))
            raise AssertionError(
                f"约束审计禁止读取求解器的通过标记 {key!r}——每条约束必须自己现算")
        return super().__getitem__(key)

    def get(self, key, default=None):
        if key in FLAG_KEYS:
            return self[key]                       # 走 __getitem__，同样拦下
        return super().get(key, default)


def load(name: str) -> NoFlags:
    """读产物并包成 NoFlags（object_pairs_hook 让每一层嵌套字典都带护栏）。"""
    if name not in _CACHE:
        path = FIG_DIR / name
        assert path.exists(), f"缺 {name}，先跑对应的 code/*.py"
        _CACHE[name] = json.loads(path.read_text(encoding="utf-8"),
                                  object_pairs_hook=lambda kv: NoFlags(kv))
    return _CACHE[name]


def touch(entity: str) -> str:
    TOUCHED.add(entity)
    return entity


def grid(problem: int) -> dict:
    touch(f"problem_{problem}_results.json")
    return load(f"problem_{problem}_results.json")["result_grid"]


def field(problem: int, key: str) -> np.ndarray:
    """result_grid 的全分辨率场，问题 4 的域外空缺填 NaN（不填 0，避免污染极值）。"""
    rows = grid(problem)[key]
    return np.asarray([[np.nan if v is None else v for v in row] for row in rows], dtype=float)


def series(problem: int, key: str) -> tuple[np.ndarray, str]:
    """取时间序列：优先 mass_series（全分辨率），退到 criterion_series（600 点抽样）。"""
    data = load(f"problem_{problem}_results.json")
    touch(f"problem_{problem}_results.json")
    for block, tier in (("mass_series", "全分辨率"), ("criterion_series", "600 点抽样")):
        src = data.get(block, {})
        if key in src:
            return np.asarray(src[key], dtype=float), f"{block}（{tier}）"
    if problem == 4:
        src = data["main"]["criterion_series"]
        if key in src:
            touch("problem_4_results.json::main")
            return np.asarray(src[key], dtype=float), "main.criterion_series（600 点抽样）"
    raise AssertionError(f"P{problem} 找不到序列 {key}")


def t_dry_block(problem: int) -> dict:
    data = load(f"problem_{problem}_results.json")
    touch(f"problem_{problem}_results.json")
    return data["t_dry"] if problem == 3 else data["main"]["t_dry"]


def sens_cases() -> dict:
    """灵敏度的 14 个算例明细（基准 + 13 扰动）。基准算例也在里头，棘轮要求它被覆盖。"""
    data = load("sensitivity_results.json")
    touch("sensitivity_results.json")
    detail = data["detail_by_case"]
    base = data["baseline_case_id"]
    assert base in detail, f"灵敏度声明基准算例 {base!r} 但明细里没有"
    assert len(detail) == int(data["n_cases"]), "灵敏度算例数与明细条数不符"
    return detail


def compare_rows() -> tuple[list, dict]:
    """模型对比的 M1—M4 行 + M1/M2 的实算明细（含各自的跨阈值区间）。"""
    data = load("model_compare_results.json")
    touch("model_compare_results.json")
    rows = list(data["table"])
    assert [r["model"] for r in rows] == ["M1", "M2", "M3", "M4"], "模型对比表不是 M1—M4"
    return rows, {"M1": data["lumped_detail"], "M2": data["frozen_detail"],
                  "M3": data["m3_detail"], "M4": data["m4_detail"]}


def record(kid: str, statement: str, scope: str, method: str, measured: dict,
           limit, ok: bool, source: str, verdict: str | None = None) -> dict:
    return {"id": kid, "statement": statement, "scope": scope, "method": method,
            "measured": measured, "limit": limit, "source": source,
            "verdict": verdict or ("PASS" if ok else "FAIL"), "ok": bool(ok)}


def k01() -> dict:
    """0 ≤ C ≤ 2.55：四问题全分辨率场 + 灵敏度 14 算例 + 模型对比 M1/M2/M3 的剖面。"""
    per = {}
    for p in (1, 2, 3, 4):
        c = field(p, "C")
        per[f"P{p}"] = {"C_min": float(np.nanmin(c)), "C_max": float(np.nanmax(c)),
                        "n_cells": int(np.count_nonzero(~np.isnan(c)))}
    lo = min(v["C_min"] for v in per.values())
    hi = max(v["C_max"] for v in per.values())

    # 灵敏度的每个扰动算例也是本模型的输出，同样受 K01 约束
    sens = {}
    for cid, d in sens_cases().items():
        vals = [float(d[k]) for k in ("C_center_end", "C_surface_end", "C_bar_end", "C_max_end")]
        sens[cid] = {"C_min": min(vals), "C_max": max(vals)}
        lo, hi = min(lo, min(vals)), max(hi, max(vals))
    # 模型对比：M2 的末刻剖面（21 列）、M3 的末刻剖面、M1 的 C̄ 全序列
    rows, detail = compare_rows()
    cmp_probe = {}
    for tag in ("M2", "M3"):
        prof = np.asarray(detail[tag]["profile_end_C"], dtype=float)
        cmp_probe[f"{tag}_profile_end"] = {"C_min": float(prof.min()), "C_max": float(prof.max()),
                                           "n_points": int(prof.size)}
        lo, hi = min(lo, float(prof.min())), max(hi, float(prof.max()))
    m1 = np.asarray(detail["M1"]["series"]["C_bar"], dtype=float)
    cmp_probe["M1_C_bar_series"] = {"C_min": float(m1.min()), "C_max": float(m1.max()),
                                    "n_points": int(m1.size)}
    lo, hi = min(lo, float(m1.min())), max(hi, float(m1.max()))
    return record("K01", "0 ≤ C(r,t) ≤ 2.55", "四问题全场逐点 + 灵敏度 14 算例 + 模型对比 M1—M3",
                  "result_grid 的 C 场逐格取 min/max（问题 4 域外空缺记 NaN，不当 0）；"
                  "灵敏度取每算例的中心/表面/均值/最大值；模型对比取 M2、M3 末刻剖面与 M1 的 C̄ 序列",
                  {"C_min_all": lo, "C_max_all": hi, "per_problem": per,
                   "n_sensitivity_cases": len(sens), "sensitivity": sens,
                   "model_compare": cmp_probe, "n_compare_rows": len(rows)},
                  [0.0, C_INIT], -1e-9 <= lo and hi <= C_INIT + 1e-9,
                  "figures/problem_{1..4}_results.json 的 result_grid.C + sensitivity_results.json "
                  "的 detail_by_case + model_compare_results.json 的剖面/序列")


def k02() -> dict:
    """28 ≤ T ≤ max_t T_a + 0.5：问题 1、2 用全场，问题 3、4 用节点温度序列。"""
    per = {}
    for p in (1, 2):
        t = field(p, "T")
        per[f"P{p}"] = {"T_min": float(np.nanmin(t)), "T_max": float(np.nanmax(t)),
                        "tier": "result_grid.T 全分辨率场"}
    for p in (3, 4):
        data = load(f"problem_{p}_results.json")
        cs = data["criterion_series"] if p == 3 else data["main"]["criterion_series"]
        touch(f"problem_{p}_results.json")
        both = np.concatenate([np.asarray(cs["T_center"], float),
                               np.asarray(cs["T_surface"], float)])
        per[f"P{p}"] = {"T_min": float(both.min()), "T_max": float(both.max()),
                        "tier": "criterion_series 的中心/表面温度（600 点抽样，非全场）"}
    ta_max = max(load(f"problem_{p}_results.json")["diagnostics"]["Ta_max_used"]
                 for p in (1, 2, 3))
    limit = ta_max + 0.5
    lo = min(v["T_min"] for v in per.values())
    hi = max(v["T_max"] for v in per.values())
    return record("K02", "28 ≤ T(r,t) ≤ max_t T_a(t) + 0.5", "问题 1、2 全场；问题 3、4 节点序列",
                  "逐格/逐点取极值，上限由各问题实测 Ta_max_used 的最大值 +0.5 现算",
                  {"T_min_all": lo, "T_max_all": hi, "Ta_max_used": ta_max,
                   "per_problem": per},
                  [T_INIT, limit], T_INIT - 1e-9 <= lo and hi <= limit,
                  "result_grid.T（P1/P2）+ criterion_series.T_*（P3/P4）")


def k03() -> dict:
    """∂C/∂r ≤ 0：逐时刻沿半径差分，等价于 max_r C = C(0,t)。"""
    per = {}
    worst = -np.inf
    for p in (1, 2, 3, 4):
        c = field(p, "C")
        d = np.diff(c, axis=1)
        rise = float(np.nanmax(d))
        n_bad = int(np.count_nonzero(d > 1e-12))
        centre_is_max = bool(np.all(np.abs(c[:, 0] - np.nanmax(c, axis=1)) < 1e-12))
        per[f"P{p}"] = {"max_radial_rise": rise, "n_rising_pairs": n_bad,
                        "center_equals_max_all_rows": centre_is_max}
        worst = max(worst, rise)
    # 灵敏度与模型对比只留了中心/表面两点，退化为「表面不高于中心」的两点判据
    sens = {}
    for cid, d in sens_cases().items():
        gap = float(d["C_surface_end"] - d["C_center_end"])
        sens[cid] = {"surface_minus_center": gap}
        worst = max(worst, gap)
    rows, _ = compare_rows()
    cmp_probe = {}
    for r in rows:
        if not r["radial_capability"]:
            continue                                  # M1 是 0 维，没有径向可言
        gap = -float(r["C_center_minus_surface"])
        cmp_probe[r["model"]] = {"surface_minus_center": gap}
        worst = max(worst, gap)
    return record("K03", "∂C/∂r ≤ 0（等价 max_r C = C(0,t)）",
                  "四问题逐时刻径向差分 + 灵敏度/模型对比的中心-表面两点判据",
                  "对 result_grid.C 沿半径作一阶差分取最大值；逐行核对中心值即行最大值；"
                  "灵敏度与模型对比按 C(R)-C(0) ≤ 0 核验",
                  {"max_radial_rise_all": worst, "per_problem": per,
                   "n_sensitivity_cases": len(sens), "sensitivity": sens,
                   "model_compare": cmp_probe},
                  1e-12, worst <= 1e-12,
                  "figures/problem_{1..4}_results.json 的 result_grid.C + sensitivity/model_compare "
                  "的中心与表面末值")


def k04() -> dict:
    """C̄(t) 单调不增：P1/P2 用逐秒全序列，P3/P4 用 600 点判据序列 + P4 表面列全序列。"""
    per = {}
    worst = -np.inf
    for p in (1, 2, 3, 4):
        arr, tier = series(p, "C_bar")
        rise = float(np.max(np.diff(arr))) if arr.size > 1 else 0.0
        per[f"P{p}"] = {"C_bar_max_rise": rise, "n_points": int(arr.size), "tier": tier,
                        "C_bar_first": float(arr[0]), "C_bar_last": float(arr[-1])}
        worst = max(worst, rise)
    surf = np.asarray(grid(4)["C_surface"], dtype=float)
    surf_rise = float(np.max(np.diff(surf)))
    per["P4_surface_column"] = {"max_rise": surf_rise, "n_points": int(surf.size),
                                "tier": "result_grid.C_surface（全分辨率）"}
    worst = max(worst, surf_rise)
    return record("K04", "C̄(t) 单调不增", "四问题逐时刻 + 问题 4 表面列",
                  "对 C̄ 序列作一阶差分取最大值；上升量 > 1e-12 即判违规",
                  {"max_rise_all": worst, "per_problem": per}, 1e-12, worst <= 1e-12,
                  "mass_series.C_bar（P1/P2 全分辨率）+ criterion_series.C_bar（P3/P4）")


def k05() -> dict:
    """ε_mass < 1%：按坐标系分档。物质坐标必须守恒，实验室坐标必须违反（推翻论证的证据）。"""
    must_hold, must_fail = {}, {}
    for p in (1, 2, 3):
        d = load(f"problem_{p}_results.json")["diagnostics"]
        touch(f"problem_{p}_results.json")
        eps = float(d["eps_mass"])
        recomputed = abs((d["m0"] - d["m_end"]) - d["flux_total"]) / d["m0"]
        assert abs(recomputed - eps) < 1e-12 or recomputed < 1e-9, \
            f"P{p} ε_mass 与 |(m0-m_end)-flux|/m0 现算值不符：{eps:.3e} vs {recomputed:.3e}"
        must_hold[f"P{p}"] = {"eps_mass": eps, "recomputed_from_m0_mend_flux": recomputed}
    p4 = load("problem_4_results.json")
    for branch in ("main", "lab", "fixed"):
        d = p4[branch]["diagnostics"]
        touch(f"problem_4_results.json::{branch}")
        eps = float(d["eps_mass"])
        recomputed = abs((d["m0"] - d["m_end"]) - d["flux_total"]) / d["m0"]
        entry = {"eps_mass": eps, "recomputed_from_m0_mend_flux": recomputed,
                 "frame": p4[branch]["frame"]}
        (must_fail if branch == "lab" else must_hold)[f"P4-{branch}"] = entry
    # 灵敏度 14 个算例各自的守恒误差；扰动改的是物性/边界，不该动守恒
    sens = {cid: float(d["eps_mass"]) for cid, d in sens_cases().items()}
    worst_sens = max(sens.values())
    # 模型对比：M2 是本文件之外独立跑的冻结物性算例；M3/M4 的值须与来源问题一致
    rows, _ = compare_rows()
    cmp_eps = {r["model"]: (None if r["eps_mass"] is None else float(r["eps_mass"]))
               for r in rows}
    tie = {"M3_vs_P3": abs(cmp_eps["M3"] - must_hold["P3"]["eps_mass"]),
           "M4_vs_P4_main": abs(cmp_eps["M4"] - must_hold["P4-main"]["eps_mass"])}
    worst_cmp = max(v for v in cmp_eps.values() if v is not None)

    worst_hold = max([v["eps_mass"] for v in must_hold.values()] + [worst_sens, worst_cmp])
    lab = must_fail["P4-lab"]["eps_mass"]
    ok = (worst_hold < 1e-2 and lab > 1e-2 and cmp_eps["M1"] is None
          and max(tie.values()) < 1e-18)
    return record("K05", "ε_mass < 1%（物质坐标）；实验室坐标按设计违反",
                  "四问题 + 问题 4 三分支 + 灵敏度 14 算例 + 模型对比 M2—M4",
                  "由 m0、m_end、flux_total 现算 |(m0-m_end)-flux|/m0 并与记录值对齐；"
                  "灵敏度与模型对比取各算例 ε_mass；并核对 M3/M4 与来源问题同值",
                  {"worst_eps_material_frames": worst_hold, "eps_lab_frame": lab,
                   "must_hold": must_hold, "must_violate": must_fail,
                   "worst_eps_sensitivity": worst_sens, "n_sensitivity_cases": len(sens),
                   "sensitivity": sens, "worst_eps_model_compare": worst_cmp,
                   "model_compare": cmp_eps, "model_compare_vs_source_abs_diff": tie,
                   "note": "实验室坐标 ε_mass 必须 > 1%，否则第 1.2 节的推翻论证失去数值证据；"
                           "M1 是 0 维模型，无径向守恒可言，其 ε_mass 必须为空"},
                  1e-2, ok,
                  "各问题 diagnostics 的 m0/m_end/flux_total + sensitivity_results.json 的 "
                  "detail_by_case.eps_mass + model_compare_results.json 的 table.eps_mass")


def k06() -> dict:
    """D 求值时 Arrhenius 指数用的是绝对温度：本文件独立复算公式，再和 Props.D 对齐。

    Props.D 的温度实参是摄氏度，内部调 to_kelvin 换算。所以「用摄氏度调用会被拦下」是
    错的探针。真正要证的是三件事：
      (a) 独立按附录 3 公式 2.4e-3·exp(-0.45/C)·exp(-3850/T_K) 复算，与 Props.D 相符；
      (b) 若漏掉 +273.15（把摄氏温度当绝对温度塞进指数），结果会差出量级——即该换算
          确有实质影响，不是可有可无的形式；
      (c) 全场最低温换算后仍 > 273 K，且 to_kelvin 对非绝对温度输入确实会拦。
    """
    props = Props(APPENDIX_OF_PROBLEM[3])
    c_probe = np.array([C_INIT])
    d_code = float(props.D(c_probe, np.array([T_INIT]))[0])
    t_k = T_INIT + KELVIN
    d_ref = D0_A3 * np.exp(-DC_A3 / C_INIT) * np.exp(-DT_ARR_A3 / t_k)     # 独立复算
    d_wrong = D0_A3 * np.exp(-DC_A3 / C_INIT) * np.exp(-DT_ARR_A3 / T_INIT)  # 漏换算
    rel = abs(d_code - d_ref) / d_ref
    ratio = d_ref / d_wrong
    blocked = False
    try:
        props.D(c_probe, np.array([-300.0]))       # -300 ℃ → 负绝对温度，必须被拦
    except AssertionError:
        blocked = True
    t_min = min(min(load(f"problem_{p}_results.json")["diagnostics"]["T_range"])
                for p in (1, 2, 3))
    ok = rel < 1e-13 and ratio > 1e3 and blocked and (t_min + KELVIN) > 273.0
    return record("K06", "D 求值时 Arrhenius 指数用绝对温度（T_K > 273）", "每次调用",
                  "独立按附录 3 解析式复算 D 并与 Props.D 对齐；再算漏掉 273.15 会差多少倍；"
                  "并核对全场最低温换算后 > 273 K",
                  {"D_from_Props": d_code, "D_recomputed_here": float(d_ref),
                   "rel_diff": rel, "D_if_kelvin_omitted": float(d_wrong),
                   "ratio_correct_over_omitted": float(ratio),
                   "sub_zero_kelvin_call_blocked": blocked,
                   "T_min_celsius": float(t_min), "T_min_kelvin": float(t_min + KELVIN)},
                  273.0, ok,
                  "code/utils.py 的 Props.D 运行期探针 + params 的 D0_A3/DC_A3/DT_ARR_A3 解析式复算")


def _rho_at_dry(appendix: int) -> float:
    """C→0 时的湿密度探针（附录 2 为常数 820）。"""
    return float(Props(appendix).rho(np.array([0.0]))[0])


# ρ(C→0) 不等于裸常数：物性求值前有 C_FLOOR=1e-6 的下截断，故极限值是 a + b·C_FLOOR。
# 三套附录的干密度相差 ≥ 90，这点 O(1e-4) 的偏置远不影响可分辨性，但容差必须按真值给。
RHO_FLOORED = {2: RHO_A2, 3: RHO3_A + RHO3_B * C_FLOOR, 4: RHO4_A + RHO4_B * C_FLOOR}


def _appendix_isolation() -> dict:
    """三个附录的干密度探针：串用会立刻在这个数上暴露。"""
    got = {a: _rho_at_dry(a) for a in RHO_FLOORED}
    return {"expected_floored": {str(a): v for a, v in RHO_FLOORED.items()},
            "bare_constants": {"2": RHO_A2, "3": RHO3_A, "4": RHO4_A},
            "C_FLOOR": C_FLOOR,
            "probed_rho_at_C0": {str(a): got[a] for a in RHO_FLOORED},
            "max_abs_err_vs_floored": max(abs(got[a] - RHO_FLOORED[a]) for a in RHO_FLOORED),
            "min_gap_between_appendices": min(abs(got[a] - got[b]) for a in got for b in got
                                              if a < b),
            "assigned_appendix": {str(p): APPENDIX_OF_PROBLEM[p] for p in (1, 2, 3, 4)}}


def k07() -> dict:
    """问题 2、3 不得混入附录 2、4 常数：探针 ρ(C→0) 必须落在附录 3 的 650(+截断偏置) 上。"""
    iso = _appendix_isolation()
    used = iso["probed_rho_at_C0"][str(APPENDIX_OF_PROBLEM[3])]
    ok = (abs(used - RHO_FLOORED[3]) < 1e-9 and abs(used - RHO_A2) > 1.0
          and abs(used - RHO4_A) > 1.0
          and APPENDIX_OF_PROBLEM[2] == APPENDIX_OF_PROBLEM[3] == 3)
    return record("K07", "问题 2、3 的物性模块不含附录 2、4 常数", "运行期物性探针",
                  "对问题 2/3 指定的附录实例化 Props 求 ρ(C→0)，与 a+b·C_FLOOR 现算值比对，"
                  "并核验其与另两套附录常数相距 > 1",
                  {"rho_dry_used_by_P2P3": used, "appendix3_floored_reference": RHO_FLOORED[3],
                   "appendix3_bare_constant": RHO3_A, "abs_err": abs(used - RHO_FLOORED[3]),
                   "appendix2_reference": RHO_A2, "appendix4_reference": RHO4_A,
                   "isolation_probe": iso}, RHO_FLOORED[3], ok,
                  "code/utils.py 的 Props（按 APPENDIX_OF_PROBLEM 实例化）+ params 的 RHO3_A/RHO3_B")


def k08() -> dict:
    """问题 4 不得混入附录 3 常数：探针 ρ(C→0) 必须落在附录 4 的 760(+截断偏置) 上。"""
    iso = _appendix_isolation()
    used = iso["probed_rho_at_C0"][str(APPENDIX_OF_PROBLEM[4])]
    ok = (abs(used - RHO_FLOORED[4]) < 1e-9 and abs(used - RHO3_A) > 1.0
          and APPENDIX_OF_PROBLEM[4] == 4)
    return record("K08", "问题 4 的物性模块不含附录 3 常数", "运行期物性探针",
                  "对问题 4 指定的附录实例化 Props 求 ρ(C→0)，与 a+b·C_FLOOR 现算值比对",
                  {"rho_dry_used_by_P4": used, "appendix4_floored_reference": RHO_FLOORED[4],
                   "appendix4_bare_constant": RHO4_A, "abs_err": abs(used - RHO_FLOORED[4]),
                   "appendix3_reference": RHO3_A, "isolation_probe": iso},
                  RHO_FLOORED[4], ok,
                  "code/utils.py 的 Props（按 APPENDIX_OF_PROBLEM 实例化）+ params 的 RHO4_A/RHO4_B")


def k09() -> dict:
    """dR/dt ≤ 0：直接对 result4 的 R(t) 全分辨率序列作差分，不读 rdot_max。"""
    g = grid(4)
    r_cm = np.asarray(g["R_cm"], dtype=float)
    d = np.diff(r_cm)
    rise = float(d.max())
    prof = load("problem_4_results.json")["radius_profile"]
    r_prof = np.asarray(prof["R_cm"], dtype=float)
    rise_prof = float(np.max(np.diff(r_prof)))
    return record("K09", "dR/dt ≤ 0 全程", "问题 4 的 R(t) 序列逐点差分",
                  "对 result_grid.R_cm（3069 点）与 radius_profile.R_cm（600 点）分别作差分取最大值",
                  {"max_rise_result_grid_cm": rise, "max_rise_radius_profile_cm": rise_prof,
                   "R_first_cm": float(r_cm[0]), "R_last_cm": float(r_cm[-1]),
                   "n_points": int(r_cm.size)}, 1e-12,
                  rise <= 1e-12 and rise_prof <= 1e-12,
                  "figures/problem_4_results.json 的 result_grid.R_cm 与 radius_profile.R_cm")


def k10() -> dict:
    """t_dry < 259200 s：问题 4 必须落在附件 2 的半径覆盖窗内，不外推。"""
    t4 = t_dry_block(4)["seconds"]
    t3 = t_dry_block(3)["seconds"]
    return record("K10", "t_dry < 259200 s（附件 2 覆盖内）", "问题 4 结束时",
                  "取问题 4 主算例的 t_dry 与附件 2 的时间跨度 ATT2_T_END_S 直接比较",
                  {"t_dry_P4_s": float(t4), "t_dry_P4_h": float(t4 / SEC_PER_HOUR),
                   "t_dry_P3_s": float(t3), "t_dry_P3_h": float(t3 / SEC_PER_HOUR),
                   "attachment2_span_s": float(ATT2_T_END_S),
                   "margin_s": float(ATT2_T_END_S - t4)},
                  float(ATT2_T_END_S), t4 < ATT2_T_END_S,
                  "figures/problem_{3,4}_results.json 的 t_dry + params 的 ATT2_T_END_S")


def k11() -> dict:
    """t_dry 落在跨阈值输出步内，且区间两端真正跨过阈值。"""
    per = {}
    ok = True
    for p in (3, 4):
        td = t_dry_block(p)
        lo_t, hi_t = float(td["bracket_t_lo_s"]), float(td["bracket_t_hi_s"])
        lo_c, hi_c = float(td["bracket_Cmax_lo"]), float(td["bracket_Cmax_hi"])
        out_dt = {3: P3_OUT_DT_S, 4: P4_OUT_DT_S}[p]
        # 用区间两端重算一次线性插值，必须与记录的 t_dry 一致
        redo = lo_t + (lo_c - C_THRESHOLD) / (lo_c - hi_c) * (hi_t - lo_t)
        good = (lo_c >= C_THRESHOLD > hi_c and lo_t <= td["seconds"] <= hi_t
                and abs((hi_t - lo_t) - out_dt) < 1e-9
                and abs(redo - td["seconds"]) < 1e-6)
        per[f"P{p}"] = {"bracket_t_s": [lo_t, hi_t], "bracket_Cmax": [lo_c, hi_c],
                        "t_dry_s": float(td["seconds"]), "interp_recomputed_s": float(redo),
                        "interp_abs_diff_s": abs(redo - td["seconds"]),
                        "bracket_width_s": hi_t - lo_t, "out_dt_s": out_dt}
        ok = ok and good
    # 灵敏度 14 个算例：阈值扰动算例用它自己的 stop_threshold，不套用 0.15
    sens, bad_sens = {}, []
    worst_sens = 0.0
    for cid, d in sens_cases().items():
        lo_t, hi_t, lo_c, hi_c = (float(v) for v in d["bracket"])
        th = float(d["stop_threshold"])
        redo = lo_t + (lo_c - th) / (lo_c - hi_c) * (hi_t - lo_t)
        dev = abs(redo - float(d["t_dry_s"]))
        sens[cid] = {"threshold": th, "bracket_t_s": [lo_t, hi_t],
                     "interp_recomputed_s": redo, "interp_abs_diff_s": dev}
        worst_sens = max(worst_sens, dev)
        if not (lo_c >= th > hi_c and abs((hi_t - lo_t) - P3_OUT_DT_S) < 1e-9 and dev < 1e-6):
            bad_sens.append(cid)
    # 模型对比的 M1、M2 也各自实算了跨阈值区间（M1 的判据退化为 C̄ < 0.15）
    _, detail = compare_rows()
    cmp_probe, bad_cmp = {}, []
    for tag in ("M1", "M2"):
        lo_t, hi_t, lo_c, hi_c = (float(v) for v in detail[tag]["bracket"])
        redo = lo_t + (lo_c - C_THRESHOLD) / (lo_c - hi_c) * (hi_t - lo_t)
        dev = abs(redo - float(detail[tag]["t_dry_s"]))
        cmp_probe[tag] = {"bracket_t_s": [lo_t, hi_t], "bracket_C": [lo_c, hi_c],
                          "interp_recomputed_s": redo, "interp_abs_diff_s": dev}
        worst_sens = max(worst_sens, dev)
        if not (lo_c >= C_THRESHOLD > hi_c and abs((hi_t - lo_t) - P3_OUT_DT_S) < 1e-9
                and dev < 1e-6):
            bad_cmp.append(tag)
    ok = ok and not bad_sens and not bad_cmp
    return record("K11", "t_dry ∈ [t_k, t_{k+1}] 且两端跨阈值",
                  "问题 3、4 + 灵敏度 14 算例 + 模型对比 M1、M2",
                  "由区间两端的 max_r C（M1 为 C̄）重算线性插值，须与记录的 t_dry 在 1e-6 s 内一致；"
                  "阈值扰动算例用其自身 stop_threshold 复算",
                  {"threshold": C_THRESHOLD, "per_problem": per,
                   "n_sensitivity_cases": len(sens), "worst_interp_abs_diff_s": worst_sens,
                   "sensitivity": sens, "sensitivity_violations": bad_sens,
                   "model_compare": cmp_probe, "model_compare_violations": bad_cmp},
                  C_THRESHOLD, ok,
                  "figures/problem_{3,4}_results.json 的 t_dry.bracket_* + sensitivity_results.json "
                  "的 detail_by_case.bracket + model_compare_results.json 的 lumped/frozen bracket"
                  "（全部自行重算插值）")


def k12() -> dict:
    """网格加倍关键输出变化 < 0.1%：由 verification 的原始逐 N 算例现算相对差。"""
    gc = load("verification_results.json")["grid_convergence"]
    touch("verification_results.json")
    tol = float(gc["tol_pct"])
    per = {}
    worst = 0.0
    for name, blk in gc.items():
        if not isinstance(blk, dict) or "cases" not in blk:
            continue
        key = blk["key_output"]
        cases = {int(c["N"]): c for c in blk["cases"]}
        base = cases[N_GRID]
        finer = cases[N_GRID * 2]
        rel = abs(finer[key] - base[key]) / abs(base[key]) * 100.0
        per[name] = {"key_output": key, "N_base": N_GRID, "N_fine": N_GRID * 2,
                     "value_base": float(base[key]), "value_fine": float(finer[key]),
                     "rel_change_pct": rel}
        worst = max(worst, rel)
    assert per, "verification_results.json 未提供逐 N 原始算例，K12 无法现算"
    return record("K12", "网格加倍关键输出变化 < 0.1%", "收敛性算例",
                  f"取 N={N_GRID} 与 N={N_GRID * 2} 的 key_output 现算相对变化，不读收敛结论字段",
                  {"worst_rel_change_pct": worst, "per_case": per}, tol, worst < tol,
                  "figures/verification_results.json 的 grid_convergence.cases（原始逐 N 值）")


def _scan_xlsx() -> dict:
    """打开四个交付 xlsx 逐格扫描：小数位、空缺位置、时间列、表面列。只读一次，供 K13—K15 共用。"""
    out = {}
    for problem, name in sorted(RESULT_FILES.items()):
        path = ROOT / name
        assert path.exists(), f"缺交付件 {name}，先跑 code/export_excel.py"
        touch(name)
        r_cm = np.asarray(grid(4)["R_cm"], dtype=float) if problem == 4 else None
        book = openpyxl.load_workbook(path, read_only=True, data_only=True)
        try:
            sheets = {}
            for ws in book.worksheets:
                rows = ws.iter_rows(values_only=True)
                head = list(next(rows))
                radii = [v for v in head[1:] if isinstance(v, (int, float))]
                surf_col = isinstance(head[-1], str)
                n_val = n_blank = n_bad_round = n_bad_blank = 0
                worst_round = 0.0
                for k, raw in enumerate(rows):
                    cells = list(raw)[1:len(head)]
                    cells += [None] * (len(head) - 1 - len(cells))
                    body = cells[:len(radii)]
                    for j, v in enumerate(body):
                        if v is None:
                            n_blank += 1
                            if r_cm is None or radii[j] <= r_cm[k] + 1e-12:
                                n_bad_blank += 1     # 域内却空缺
                            continue
                        n_val += 1
                        err = abs(float(v) - round(float(v), DECIMALS))
                        worst_round = max(worst_round, err)
                        n_bad_round += err > 1e-12
                        if r_cm is not None and radii[j] > r_cm[k] + 1e-12:
                            n_bad_blank += 1         # 域外却有值
                    if surf_col:
                        v = cells[len(radii)]
                        assert v is not None, f"{name} 第 {k + 1} 行表面列为空"
                        n_val += 1
                        worst_round = max(worst_round,
                                          abs(float(v) - round(float(v), DECIMALS)))
                sheets[ws.title] = {"n_rows": k + 1, "n_radii_cols": len(radii),
                                    "has_surface_col": surf_col,
                                    "surface_header": head[-1] if surf_col else None,
                                    "n_values": n_val, "n_blanks": n_blank,
                                    "n_bad_rounding": int(n_bad_round),
                                    "n_blank_rule_violations": int(n_bad_blank),
                                    "worst_rounding_err": worst_round}
            out[name] = sheets
        finally:
            book.close()
    return out


def k13(scan: dict) -> dict:
    """四位小数、无空缺（问题 4 域外列除外）：直接扫交付 xlsx，不看 JSON。"""
    worst = max(s["worst_rounding_err"] for f in scan.values() for s in f.values())
    n_bad = sum(s["n_bad_rounding"] for f in scan.values() for s in f.values())
    n_val = sum(s["n_values"] for f in scan.values() for s in f.values())
    blanks_non_p4 = sum(s["n_blanks"] for n, f in scan.items() for s in f.values()
                        if n != RESULT_FILES[4])
    return record("K13", "result1—4 全部数值四位小数、无空缺（问题 4 域外列除外）", "文件校验",
                  "openpyxl 逐格回读交付件，现算 |v - round(v,4)| 的最大值与非问题 4 文件的空缺数",
                  {"n_values_scanned": n_val, "worst_rounding_err": worst,
                   "n_cells_not_4dp": n_bad, "n_blanks_outside_P4": blanks_non_p4,
                   "per_file": scan}, 1e-12,
                  worst <= 1e-12 and n_bad == 0 and blanks_non_p4 == 0,
                  "工作区根目录的 result1—4.xlsx（逐格扫描）")


def k14(scan: dict) -> dict:
    """result4 域外单元为空、不填 0：逐格核对「空缺 <=> r_j > R(t)」。"""
    name = RESULT_FILES[4]
    sheet = next(iter(scan[name].values()))
    declared = grid(4)["n_blank_cells"]
    return record("K14", "result4 中 r_j > R(t) 的单元为空，不填 0", "文件校验",
                  "以 result_grid.R_cm 为准逐格判定域内外，与 xlsx 的空缺位置做双向核对",
                  {"n_blanks_in_xlsx": sheet["n_blanks"], "n_blanks_declared_json": declared,
                   "n_rule_violations": sheet["n_blank_rule_violations"],
                   "n_values_in_xlsx": sheet["n_values"],
                   "surface_header": sheet["surface_header"]}, 0,
                  sheet["n_blank_rule_violations"] == 0 and sheet["n_blanks"] == declared
                  and sheet["n_blanks"] > 0,
                  f"{name}（逐格扫描）+ figures/problem_4_results.json 的 result_grid.R_cm")


def k15(scan: dict) -> dict:
    """论文表 6 末列取 s=1 节点值：与表体同分辨率的表面列逐时刻对齐。"""
    data = load("problem_4_results.json")
    t6 = data["table6_C"]
    surf = np.asarray(grid(4)["C_surface"], dtype=float)
    diffs = []
    for i, tt in enumerate(t6["times_s"][:-1]):          # 末行是 t_dry，不在 60 s 网格上
        k = int(round(tt / P4_OUT_DT_S)) - 1
        diffs.append(abs(t6["values"][i][-1] - float(surf[k])))
    tail = abs(t6["values"][-1][-1] - data["main"]["t_dry"]["C_surface_at_t_dry"])
    name = RESULT_FILES[4]
    header = next(iter(scan[name].values()))["surface_header"]
    worst = max(diffs + [tail])
    return record("K15", "论文表 6 末列取 s=1 节点值", "表格校验",
                  "把表 6 末列与 result_grid.C_surface 按时刻对齐逐行作差；末行与 t_dry 处剖面比对",
                  {"n_rows_checked": len(diffs) + 1, "worst_abs_diff": worst,
                   "t_dry_row_abs_diff": tail, "xlsx_surface_header": header,
                   "json_surface_column": grid(4)["surface_column"],
                   "rule": grid(4)["surface_rule"]}, 1e-12,
                  worst < 1e-12 and isinstance(header, str) and "药材表面" in header,
                  "figures/problem_4_results.json 的 table6_C 与 result_grid.C_surface + result4.xlsx 表头")


def k16() -> dict:
    """C_a 与 C 在代码中是独立命名对象：扫源码，禁止把空气侧含水率写回材料场变量。"""
    hits, bad = {}, []
    for path in sorted((ROOT / "code").rglob("*.py")):
        src = path.read_text(encoding="utf-8", errors="replace")
        body = "\n".join("" if ln.strip().startswith("#") else ln for ln in src.splitlines())
        n_ca = len(re.findall(r"\bc_amb\b|\bCa\b|\bscale_Ca\b|\bC_a\b", body))
        n_c = len(re.findall(r"\bC_field\b|\bC_final\b|\bc_bar\b|\bC_bar\b", body))
        for m in re.finditer(r"^\s*(C|C_field|C_final)\s*=\s*[^=].*?\b(c_amb|C_a|Ca)\b",
                             body, re.M):
            bad.append(f"{path.name}:{body.count(chr(10), 0, m.start()) + 1}")
        if n_ca or n_c:
            hits[path.name] = {"n_ambient_refs": n_ca, "n_material_refs": n_c}
    return record("K16", "C_a 与 C 在代码中为独立命名对象", "源码扫描",
                  "扫 code/*.py，禁止出现把空气侧含水率直接赋给材料场变量的语句",
                  {"n_files_referencing_both": len(hits), "violations": bad,
                   "per_file": hits}, 0, not bad,
                  "code/*.py 源码扫描（已剔除注释行）")


def discover_baselines() -> dict:
    """按关键词自动发现基准/对照算例，逐个核对是否已被本次审计覆盖。

    覆盖的判定：该基准所属的产物文件（或 文件::分支）出现在 TOUCHED 里。以后新增分支
    却不补审计，这里会直接报出来——这是防「审计随产物腐化」的棘轮。
    """
    found = {}
    for path in sorted(FIG_DIR.glob("*_results.json")):
        if path.name == SELF_PRODUCT:
            continue                       # 本文件自己的审计报告，扫它只会扫到本次审计的措辞
        data = json.loads(path.read_text(encoding="utf-8"))
        labels = []

        def walk(node, trail: str) -> None:
            if isinstance(node, dict):
                for key, val in node.items():
                    if BASELINE_RE.search(str(key)):
                        labels.append({"where": f"{trail}.{key}", "kind": "key",
                                       "value": val if isinstance(val, (str, int, float))
                                       else type(val).__name__})
                    walk(val, f"{trail}.{key}")
            elif isinstance(node, list):
                for item in node[:200]:
                    walk(item, f"{trail}[]")
            elif isinstance(node, str) and BASELINE_RE.search(node) and len(node) < 60:
                labels.append({"where": trail, "kind": "value", "value": node})

        walk(data, "")
        branches = [b for b in ("main", "lab", "fixed") if isinstance(data.get(b), dict)]
        entities = [f"{path.name}::{b}" for b in branches] or [path.name]
        covered = [e for e in entities if e in TOUCHED or path.name in TOUCHED]
        found[path.name] = {"n_baseline_mentions": len(labels),
                            "mentions": labels[:12], "entities": entities,
                            "covered_entities": covered,
                            "audited": bool(covered) if labels else None}
    # 只要一个文件出现了基准/对照字样，它就必须被至少一条约束真正读过
    gaps = [name for name, v in found.items() if v["n_baseline_mentions"] and not v["covered_entities"]]
    return {"keyword": BASELINE_RE.pattern, "n_files_scanned": len(found),
            "self_product_excluded": SELF_PRODUCT,
            "per_file": found, "uncovered_files": gaps,
            "touched_entities": sorted(TOUCHED),
            "ok": not gaps}


def guard_probe() -> dict:
    """自证：本文件读不到产物里的通过标记——试着读一次，必须被 NoFlags 护栏拦下。

    这不是正则自查（那会把本文件自己 record 字典的 ["ok"] 也算成违规），而是运行期证明：
    K01—K16 全程走的是同一个 load()，取用 capability_verdict / status 之类会直接抛。
    """
    probes = []
    for name in sorted(_CACHE):
        for key in ("capability_verdict", "all_pass"):
            if key not in dict.keys(_CACHE[name]):
                continue
            try:
                _CACHE[name][key]
                probes.append({"file": name, "key": key, "blocked": False})
            except AssertionError:
                probes.append({"file": name, "key": key, "blocked": True})
    nested = None
    if "problem_1_results.json" in _CACHE:
        cv = dict.__getitem__(_CACHE["problem_1_results.json"], "capability_verdict")
        inner = dict.__getitem__(cv, "P1-C1")           # 绕过护栏取到内层，再证内层也带护栏
        try:
            inner["status"]
            nested = False
        except AssertionError:
            nested = True
    ok = bool(probes) and all(p["blocked"] for p in probes) and nested is True
    return {"mechanism": "load() 以 object_pairs_hook 把每层字典包成 NoFlags，"
                         "__getitem__/get 命中标记字段即抛 AssertionError",
            "guarded_keys": sorted(FLAG_KEYS), "n_probes": len(probes), "probes": probes,
            "nested_level_guarded": nested, "n_guard_hits_this_run": len(GUARD_HITS),
            "ok": ok,
            "note": "K01—K16 的判定全部由本文件现算，产物里的 PASS 字段读不到，也没有被读"}


def _walk_eps(node, trail: str, out: list) -> None:
    """只收 eps_mass 数值，不碰任何通过标记字段。"""
    if isinstance(node, dict):
        for key, val in node.items():
            if key == "eps_mass" and isinstance(val, (int, float)):
                out.append((f"{trail}.{key}", float(val)))
            _walk_eps(val, f"{trail}.{key}", out)
    elif isinstance(node, list):
        for i, item in enumerate(node):
            _walk_eps(item, f"{trail}[{i}]", out)


def _cmp(want: float, got: float, tol: float = 1e-12) -> dict:
    dev = abs(float(got) - float(want))
    rel = dev / abs(float(want)) if float(want) != 0.0 else dev
    return {"want": float(want), "got": float(got), "abs_diff": dev,
            "rel_diff": rel, "ok": bool(rel <= tol)}


def audit_summary(n_constraints: int, n_failed: int) -> dict:
    """复核 figures/all_results.json：顶层标量与方向探针逐条从来源产物独立复算。

    汇总文件是下游 logic_audit 与论文写作直接引用的口径；它一旦与来源产物脱节，后面每一
    步都会引用错数。故这里不采信它的任何自述，只问「能否从来源复算出同一个数」，且刻意
    走与 main.py 不同的取值路径（如由 seconds 除 3600 复算 hours、由 detail_by_case 复算
    灵敏度基准、由 table 复算 M1/M2），使复算不是把同一行代码抄第二遍。

    新增顶层标量而不在此登记复算规则 → unregistered 非空 → 直接判 FAIL，防汇总口径腐化。
    """
    summ = load("all_results.json")
    touch("all_results.json")
    p1, p2, p3, p4 = (load(f"problem_{i}_results.json") for i in (1, 2, 3, 4))
    ver, sens, cmp_ = (load(f"{n}_results.json") for n in
                       ("verification", "sensitivity", "model_compare"))
    for name in ("verification_results.json", "sensitivity_results.json",
                 "model_compare_results.json"):
        touch(name)
    for branch in ("main", "lab", "fixed"):
        touch(f"problem_4_results.json::{branch}")
    td3, td4 = p3["t_dry"], p4["main"]["t_dry"]

    c_all = [field(i, "C") for i in (1, 2, 3, 4)]
    t_vals = [np.asarray(p["result_grid"]["T"], dtype=float).ravel() for p in (p1, p2)]
    for src in (p3["criterion_series"], p4["main"]["criterion_series"]):
        t_vals += [np.asarray(src["T_center"], dtype=float),
                   np.asarray(src["T_surface"], dtype=float)]
    eps_occ: list = []
    _walk_eps(ver, "", eps_occ)
    assert eps_occ, "verification_results.json 里没有 eps_mass，无法复算守恒极值"
    worst_grid = 0.0
    for blk in ver["grid_convergence"].values():
        if not isinstance(blk, dict) or "cases" not in blk:
            continue
        key = blk["key_output"]
        cases = {int(c["N"]): c for c in blk["cases"]}
        worst_grid = max(worst_grid, abs(cases[N_GRID * 2][key] - cases[N_GRID][key])
                         / abs(cases[N_GRID][key]) * 100.0)
    detail = sens["detail_by_case"]
    by_model = {r["model"]: float(r["t_dry_h"]) for r in cmp_["table"]}
    lab_h, main_h = float(p4["lab"]["t_dry"]["hours"]), float(td4["hours"])

    rules = {
        "C_min_all": min(float(np.nanmin(c)) for c in c_all),
        "C_max_all": max(float(np.nanmax(c)) for c in c_all),
        "T_min_all": min(float(v.min()) for v in t_vals),
        "T_max_all": max(float(v.max()) for v in t_vals),
        "t_dry_Q3_s": float(td3["seconds"]),
        "t_dry_Q3_h": float(td3["seconds"]) / SEC_PER_HOUR,
        "t_dry_Q3_days": float(td3["seconds"]) / SEC_PER_HOUR / HOURS_PER_DAY,
        "t_dry_Q4_s": float(td4["seconds"]),
        "t_dry_Q4_h": float(td4["seconds"]) / SEC_PER_HOUR,
        "C_bar_at_t_dry_Q3": float(td3["C_bar_at_t_dry"]),
        "C_bar_at_t_dry_Q4": float(td4["C_bar_at_t_dry"]),
        "R_ratio_Q4": float(td4["R_ratio"]),
        "mass_balance_error": max(float(x["diagnostics"]["eps_mass"])
                                  for x in (p1, p2, p3, p4["main"])),
        "mass_balance_error_lab_frame_Q4": float(p4["lab"]["diagnostics"]["eps_mass"]),
        "grid_doubling_max_change_pct": worst_grid,
        "mass_conservation_eps_max_all_cases": max(v for _, v in eps_occ),
        "t_dry_Q4_fixed_domain_h": float(p4["fixed"]["t_dry"]["hours"]),
        "frame_rel_diff_Q4": (lab_h - main_h) / main_h,
        "sensitivity_baseline_t_dry_h": float(detail[sens["baseline_case_id"]]["t_dry_h"]),
        "t_dry_M1_lumped_h": by_model["M1"],
        "t_dry_M2_frozen_h": by_model["M2"],
    }
    ints = {"n_subproblems": len(summ["per_problem"]),
            "sensitivity_n_cases": len(detail)}
    # 这两个计数来自本次审计自身，在此复算是循环论证；且 main.py 第一遍汇总刻意不写它们。
    exempt = {"n_constraints_checked": n_constraints, "n_constraints_failed": n_failed}
    numeric = {k for k, v in summ.items()
               if isinstance(v, (int, float)) and not isinstance(v, bool)}
    ints.update({k: v for k, v in exempt.items() if k in numeric})

    checks = {k: _cmp(want, summ[k]) for k, want in rules.items() if k in numeric}
    for k, want in ints.items():
        got = summ[k] if k in numeric else None
        checks[k] = {"want": int(want), "got": None if got is None else int(got),
                     "abs_diff": None, "rel_diff": None,
                     "ok": got is not None and int(got) == int(want)}
    known = set(rules) | set(ints)
    unregistered = sorted(numeric - known - set(exempt))
    missing = sorted(known - numeric)

    def sgn(x: float) -> int:
        return 0 if x == 0.0 else (1 if x > 0.0 else -1)

    def endsign(arr) -> int:
        a = np.asarray(arr, dtype=float)
        return sgn(float(a[-1] - a[0]))

    cs3 = p3["criterion_series"]
    hot = np.asarray(cs3["T_center"], dtype=float) > 49.9      # 近等温平台段，见 main.py 注释
    cp = np.asarray(cs3["C_max"], dtype=float)[hot]
    dp = np.asarray(cs3["D_center"], dtype=float)[hot]
    rel_by_case = {r["case_id"]: r["rel_change_pct"] for r in sens["table"]}
    want_signs = {
        ("t", "max_r C"): endsign(cs3["C_max"]),
        ("t", "volume_mean_C"): endsign(cs3["C_bar"]),
        ("r", "C"): endsign(p3["result_grid"]["C"][-1]),
        ("t", "R"): endsign(p4["result_grid"]["R_cm"]),
        ("C", "D"): sgn(float((dp[-1] - dp[0]) / (cp[-1] - cp[0]))),
        ("D_scale", "t_dry"): sgn(rel_by_case["D_up"]),
        ("C_threshold", "t_dry"): sgn(rel_by_case["Cth_up"]),
        ("h_m_scale", "t_dry"): sgn(rel_by_case["hm_up"]),
        ("C_a_scale", "t_dry"): sgn(rel_by_case["Ca_up"]),
        ("shrinkage_enabled", "t_dry"): sgn(main_h - float(p4["fixed"]["t_dry"]["hours"])),
        ("model_fidelity", "t_dry"): sgn(by_model["M3"] - by_model["M1"]),
    }
    probes = list(summ["logic_probes"]["monotonic"])
    probe_rows = []
    for p in probes:
        redone = want_signs.get((p["more"], p["then"]))
        decl, exp = int(p["observed_sign"]), int(p["expect_sign"])
        agree = (decl > 0) == (exp > 0)
        probe_rows.append({"more": p["more"], "then": p["then"], "declared_sign": decl,
                           "redone_sign": None if redone is None else int(redone),
                           "expect_sign": exp, "agrees_with_expect": agree,
                           "ok": bool(redone is not None and int(redone) == decl and agree)})
    n_probe_bad = sum(1 for r in probe_rows if not r["ok"])

    ok = (not unregistered and not missing and all(c["ok"] for c in checks.values())
          and n_probe_bad == 0 and len(probes) == len(want_signs))
    return {"file": "all_results.json", "n_scalars_rechecked": len(checks),
            "n_unregistered": len(unregistered), "unregistered": unregistered,
            "missing_from_summary": missing, "scalar_checks": checks,
            "audit_sourced_keys": sorted(exempt),
            "audit_sourced_note": ("n_constraints_checked/n_constraints_failed 出自本次审计自身，"
                                   "在此复算属循环论证；main.py 审计前那一遍汇总刻意不写这两个键，"
                                   "若已写（对终版文件单独跑审计）则与本次计数逐一核对"),
            "n_probes": len(probes), "n_probe_rules": len(want_signs),
            "n_probes_mismatch": n_probe_bad, "probes": probe_rows,
            "n_bounds_probes": len(summ["logic_probes"]["bounds"]),
            "eps_occurrences_walked": len(eps_occ), "tol_rel": 1e-12, "ok": ok,
            "source": ("figures/problem_1—4/verification/sensitivity/model_compare"
                       "_results.json 的原始块（刻意走与 main.py 不同的取值路径）")}


def main() -> int:
    ensure_dirs()
    scan = _scan_xlsx()
    results = [k01(), k02(), k03(), k04(), k05(), k06(), k07(), k08(), k09(), k10(),
               k11(), k12(), k13(scan), k14(scan), k15(scan), k16()]
    assert [r["id"] for r in results] == [f"K{i:02d}" for i in range(1, 17)], \
        "K01—K16 必须逐条出现且顺序一致"

    n_fail = sum(1 for r in results if not r["ok"])
    summary = audit_summary(len(results), n_fail)      # 须在棘轮之前跑，好把汇总记进 TOUCHED
    baselines = discover_baselines()
    selfscan = guard_probe()
    ok_all = n_fail == 0 and baselines["ok"] and selfscan["ok"] and summary["ok"]

    print(f"{'约束':<6}{'判定':<8}{'量测':<46}{'限值':>14}")
    for r in results:
        m = r["measured"]
        head = next((f"{k}={v:.4g}" if isinstance(v, float) else f"{k}={v}"
                     for k, v in m.items() if isinstance(v, (int, float))), "见 JSON")
        lim = r["limit"]
        lim_s = f"{lim:.6g}" if isinstance(lim, float) else str(lim)
        print(f"{r['id']:<6}{r['verdict']:<8}{head:<46}{lim_s:>14}")
    k05m = results[4]["measured"]
    print(f"[CA] K05 分档：物质坐标最差 ε_mass={k05m['worst_eps_material_frames']:.3e}；"
          f"实验室坐标 ε_mass={k05m['eps_lab_frame']:.3e}（按设计违反，是推翻论证的证据）")
    k11m = results[10]["measured"]
    print(f"[CA] 扩展覆盖：灵敏度 {k11m['n_sensitivity_cases']} 个算例（含基准）与模型对比 M1—M4 "
          f"已并入 K01/K03/K05/K11；跨阈值插值复算最大偏差 "
          f"{k11m['worst_interp_abs_diff_s']:.2e} s")
    worst_scalar = max((c["rel_diff"] for c in summary["scalar_checks"].values()
                        if c["rel_diff"] is not None), default=0.0)
    print(f"[CA] 汇总复核 all_results.json：{summary['n_scalars_rechecked']} 个标量由来源复算，"
          f"最大相对偏差 {worst_scalar:.2e}；未登记 {summary['n_unregistered']} 个；"
          f"方向探针 {summary['n_probes']} 条复算不符 {summary['n_probes_mismatch']} 条")
    print(f"[CA] 基准/对照自动发现：扫 {baselines['n_files_scanned']} 个产物，"
          f"未覆盖 {len(baselines['uncovered_files'])} 个；实际读过的实体 "
          f"{len(baselines['touched_entities'])} 个")
    print(f"[CA] 通过标记运行期护栏：{selfscan['n_probes']} 次探针"
          f"{'全部被拦' if selfscan['ok'] else '存在漏网'}，嵌套层"
          f"{'亦受护' if selfscan['nested_level_guarded'] else '未受护'}")
    print(f"[CA] K01—K16 共 {len(results)} 条，违规 {n_fail} 条")

    payload = {"problem": "cross", "n_constraints": len(results), "n_failed": n_fail,
               "all_pass": ok_all, "constraints": results,
               "baseline_discovery": baselines, "flag_guard": selfscan,
               "summary_recheck": summary,
               "provenance": ("每条约束由本文件从 result1—4.xlsx、result_grid、原始诊断标量现算；"
                              "产物中的 PASS/verdict 字段被运行期护栏挡住，读不到也没被读。")}
    save_json(FIG_DIR / "constraint_audit_results.json", jsonable(payload))
    if not ok_all:
        print("[CA] 存在未通过项，见 figures/constraint_audit_results.json")
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
