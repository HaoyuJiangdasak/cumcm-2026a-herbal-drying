# -*- coding: utf-8 -*-
"""数值解可信性验证（MODELING_REPORT 第 10 节、能力 PA-C1、约束 K12）。

四类独立证据，全部由本文求解器实算：
  1. 网格收敛性  ——  问题 1（固定 Δt 隔离空间误差）、问题 3、问题 4；
  2. 时间步收敛性 —— 问题 3 平台段 Δt3、问题 1 一阶精度序列 + Richardson 外推；
  3. 全局质量守恒 —— 每个算例的 ε_mass；
  4. 解析解交叉校验 —— 常物性、T_a 冻结的退化算例对 Bessel 级数解。
外加 10.5 节 V1—V6 六个结构性退化算例。

每个算例独立落盘检查点，可分阶段运行：
    python _utils/run_compute.py code/verification.py --stage q1grid
    ... q3grid q4grid q1time q3time bessel degenerate v6
    python _utils/run_compute.py code/verification.py --stage assemble
不带 --stage 时按顺序跑全部阶段（墙钟约 25 min）。
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_utils"))
from compute_checkpoint import CaseCheckpoint  # noqa: E402

from params import *  # noqa: F401,F403,E402
from utils import (Props, bessel_roots, bessel_temperature, ensure_dirs, jsonable,  # noqa: E402
                   make_ambient, run_case, save_json)

ROOT = Path(__file__).resolve().parents[1]
DEPS = ["code/utils.py", "code/params.py", "user_data/附件1.xlsx", "user_data/附件2.xlsx"]
STAGES = ["q1grid", "q3grid", "q4grid", "q1time", "q3time", "bessel", "degenerate", "v6"]


def _scalars(res: dict, *, want_tdry: bool) -> dict:
    """从 run_case 结果抽取标量指纹（检查点只存这些，改判据不必重算）。"""
    t_fin = np.asarray(res["T_final"], dtype=float)
    c_fin = np.asarray(res["C_final"], dtype=float)
    out = {
        "N": int(res["N"]), "n_steps": int(res["n_steps"]),
        "wall_s": float(res["wall_s"]), "eps_mass": float(res["eps_mass"]),
        "m0": float(res["m0"]), "m_end": float(res["m_end"]),
        "flux_total": float(res["flux_total"]),
        "picard_mean": float(res["picard_mean"]), "picard_max": int(res["picard_max"]),
        "T_center": float(t_fin[0]), "T_surface": float(t_fin[-1]),
        "C_center": float(c_fin[0]), "C_surface": float(c_fin[-1]),
        "C_max": float(c_fin.max()), "C_min": float(c_fin.min()),
        "C_bar": float(res["C_bar"][-1]), "t_end_reached": float(res["t_end_reached"]),
    }
    if want_tdry:
        assert res["t_dry"] is not None, "该算例未在搜索窗内达标"
        out["t_dry_s"] = float(res["t_dry"])
        out["t_dry_h"] = float(res["t_dry"] / SEC_PER_HOUR)
        out["bracket"] = [float(v) for v in res["bracket"]]
    return out


def _validate_scalars(data: dict) -> bool:
    """检查点校验：只核对已写出的标量自洽，不重跑求解器。"""
    assert data["N"] >= 20 and data["n_steps"] > 0, "网格或步数非法"
    assert np.isfinite(data["eps_mass"]) and data["eps_mass"] >= 0.0, "ε_mass 非法"
    assert data["C_min"] >= -1e-9 and data["C_max"] <= C_INIT + 1e-9, \
        f"K01 含水率越界 [{data['C_min']},{data['C_max']}]"
    assert data["T_center"] >= T_INIT - 1e-9 and data["T_surface"] >= T_INIT - 1e-9, \
        "K02 温度低于初值"
    assert data["m0"] > 0.0 and data["m_end"] > 0.0, "干物质量必须为正"
    assert data["picard_max"] <= PICARD_MAX, "Picard 迭代触顶，未收敛"
    if "t_dry_s" in data:
        lo, hi, c_lo, c_hi = data["bracket"]
        assert c_lo >= C_THRESHOLD > c_hi, "K11 区间未跨阈值"
        assert lo < data["t_dry_s"] < hi, "K11 t_dry 不在区间内"
    return True


def _run(key: str, cfg: dict, *, want_tdry: bool, extra=None) -> dict:
    """跑一个算例并落盘；已有可用检查点时直接复用。extra(res) 可追加算例专属标量。"""
    ckpt = CaseCheckpoint(
        ROOT, key, dependencies=list(DEPS),
        parameters={"cfg": cfg, "picard_tol": PICARD_TOL, "gauss_n": GAUSS_N,
                    "numpy": np.__version__, "want_tdry": bool(want_tdry),
                    "extra": getattr(extra, "__name__", None)},
    )
    cached = ckpt.load(_validate_scalars)
    if cached is not None:
        print(f"  [ckpt] 复用 {key}")
        return cached
    print(f"  [run ] {key}")
    res = run_case(cfg)
    data = _scalars(res, want_tdry=want_tdry)
    if extra is not None:
        data.update(extra(res))
    ckpt.save(data, _validate_scalars)
    return data


def _rel_pct(new: float, old: float) -> float:
    return abs(new - old) / abs(old) * 100.0


# ---------------------------------------------------------------- 10.1 网格收敛

def _q1(n: int, dt: float) -> dict:
    """问题 1 的 1800 s 终态算例；键按 (N, Δt) 规范化，网格与步长两条序列共用检查点。"""
    cfg = {"appendix": APPENDIX_OF_PROBLEM[1], "N": n, "out_dt": P1_WINDOW_S,
           "t_end": P1_WINDOW_S, "dt1": dt, "record": "last", "tag": f"V-q1-N{n}dt{dt:g}"}
    row = _run(f"2026A/verify/q1/N{n}dt{dt:g}", cfg, want_tdry=False)
    return dict(row, n_grid=n, dt_s=float(dt))


def _q3(n: int, dt3: float) -> dict:
    """问题 3 的 t_dry 算例；out_dt 不得小于内部步长，故 Δt3=120 s 时输出步同步放宽。"""
    out_dt = max(P3_OUT_DT_S, dt3)
    cfg = {"appendix": APPENDIX_OF_PROBLEM[3], "N": n, "out_dt": out_dt,
           "t_end": T_MAX_SEARCH_S, "stop_threshold": C_THRESHOLD, "dt1": DT_STAGE1,
           "dt2": DT_STAGE2, "dt3": dt3, "record": "last", "tag": f"V-q3-N{n}dt{dt3:g}"}
    row = _run(f"2026A/verify/q3/N{n}dt{dt3:g}", cfg, want_tdry=True)
    return dict(row, n_grid=n, dt3_s=float(dt3), out_dt_s=out_dt)


def _q4(n: int, dt3: float) -> dict:
    """问题 4（物质坐标 + 附件 2 收缩）的 t_dry 算例。"""
    cfg = {"appendix": APPENDIX_OF_PROBLEM[4], "N": n, "out_dt": P4_OUT_DT_S,
           "t_end": T_MAX_SEARCH_S, "stop_threshold": C_THRESHOLD, "dt1": DT_STAGE1,
           "dt2": DT_STAGE2, "dt3": dt3, "shrink": True, "frame": "material",
           "record": "last", "tag": f"V-q4-N{n}dt{dt3:g}"}
    row = _run(f"2026A/verify/q4/N{n}dt{dt3:g}", cfg, want_tdry=True)
    return dict(row, n_grid=n, dt3_s=float(dt3))


def stage_q1grid() -> dict:
    """问题 1：1800 s 终态，Δt 固定 0.1 s 以隔离空间误差，N=20/40/80/160。"""
    dt = min(P1_DT_CONV)
    return {"fixed_dt_s": dt, "t_eval_s": P1_WINDOW_S, "key_output": "C_surface",
            "cases": [_q1(n, dt) for n in GRID_CONV]}


def stage_q3grid() -> dict:
    """问题 3：生产时间步下 N=20/40/80 的 t_dry。"""
    return {"fixed_dt3_s": DT_STAGE3, "key_output": "t_dry_h",
            "cases": [_q3(n, DT_STAGE3) for n in GRID_CONV_TDRY]}


def stage_q4grid() -> dict:
    """问题 4：生产时间步下 N=20/40/80 的 t_dry。"""
    return {"fixed_dt3_s": DT_STAGE3, "key_output": "t_dry_h",
            "cases": [_q4(n, DT_STAGE3) for n in GRID_CONV_TDRY]}


# ---------------------------------------------------------------- 10.2 时间步收敛

def stage_q1time() -> dict:
    """问题 1（N=20 固定以隔离时间误差）：Δt=2/1/0.5/0.25/0.1 s，验证一阶精度。"""
    n = min(GRID_CONV)
    return {"fixed_N": n, "t_eval_s": P1_WINDOW_S, "key_output": "T_center",
            "cases": [_q1(n, dt) for dt in P1_DT_CONV]}


def stage_q3time() -> dict:
    """问题 3：平台段 Δt3=120/60/30/10 s 的 t_dry（前两段保持生产步长）。"""
    return {"fixed_N": N_GRID, "key_output": "t_dry_h",
            "cases": [_q3(N_GRID, dt3) for dt3 in DT_CONV]}


# ------------------------------------------------- 10.4 解析解交叉校验（含 V3）

def _bessel_ref(t_a: float) -> dict:
    """附录 2 常物性、Robin 边界、T_a 阶跃的 Bessel 级数解（前 BESSEL_TERMS 项）。"""
    props = Props(APPENDIX_OF_PROBLEM[1])
    zero = np.zeros(1)
    lam = float(props.lam(zero)[0])
    rho = float(props.rho(zero)[0])
    cp = float(props.cp(zero)[0])
    alpha = lam / (rho * cp)
    bi = H_HEAT * R0 / lam
    zeta = bessel_roots(bi, BESSEL_TERMS)
    kw = dict(bi=bi, alpha=alpha, t_amb=t_a, t_init=T_INIT, n_terms=BESSEL_TERMS)
    t_ana = bessel_temperature([0.0, 1.0], P1_WINDOW_S, **kw)
    return {"Ta_const": float(t_a), "Bi_h": bi, "alpha_m2_s": alpha, "lambda": lam,
            "zeta_first3": [float(v) for v in zeta[:3]], "n_terms": BESSEL_TERMS,
            "T_center_analytic": float(t_ana[0]), "T_surface_analytic": float(t_ana[1])}


def _bessel_case(n: int, dt: float, t_a: float) -> dict:
    cfg = {"appendix": APPENDIX_OF_PROBLEM[1], "N": n, "out_dt": P1_WINDOW_S,
           "t_end": P1_WINDOW_S, "dt1": dt, "const_Ta": t_a, "record": "last",
           "tag": f"V-bessel-N{n}dt{dt:g}"}
    key = f"2026A/verify/bessel/N{n}dt{dt:g}Ta{t_a:.4f}"
    return dict(_run(key, cfg, want_tdry=False), n_grid=n, dt_s=float(dt))


def stage_bessel() -> dict:
    """V3：T_a 冻结为平台实测均值、物性取附录 2，对 Bessel 解逐步收紧 Δt。

    附加固定 Δt 而加密网格的两个算例，检验残差是否只由时间离散贡献。
    """
    t_a = PLATEAU_TA
    ref = _bessel_ref(t_a)
    rows = []
    for dt in BESSEL_DT:
        row = _bessel_case(BESSEL_N, dt, t_a)
        row["err_center"] = abs(row["T_center"] - ref["T_center_analytic"])
        row["err_surface"] = abs(row["T_surface"] - ref["T_surface_analytic"])
        rows.append(row)
    grid_rows = []
    dt_fix = min(BESSEL_DT)
    for n in (N_GRID, BESSEL_N, max(GRID_CONV)):
        row = _bessel_case(n, dt_fix, t_a)
        row["err_center"] = abs(row["T_center"] - ref["T_center_analytic"])
        grid_rows.append(row)
    ref_round = _bessel_ref(round(t_a, 0))
    return {"analytic": ref, "cases": rows, "grid_independence": grid_rows,
            "fixed_dt_for_grid_scan_s": dt_fix,
            "analytic_at_rounded_Ta": {"Ta_const": ref_round["Ta_const"],
                                       "T_center": ref_round["T_center_analytic"],
                                       "T_surface": ref_round["T_surface_analytic"]}}


# ------------------------------------------------- 10.5 结构性退化算例 V1/V2/V4/V5

def _field_stats(res: dict) -> dict:
    """退化算例需要的场级统计（检查点只存标量，故在此压缩）。"""
    c = np.asarray(res["C_final"], dtype=float)
    t = np.asarray(res["T_final"], dtype=float)
    bar = np.asarray(res["C_bar"], dtype=float)
    return {"C_dev_from_init": float(np.max(np.abs(c - C_INIT))),
            "T_dev_from_init": float(np.max(np.abs(t - T_INIT))),
            "C_spread": float(c.max() - c.min()),
            "C_bar_max_rise": float(np.max(np.diff(bar))) if bar.size > 1 else 0.0}


def stage_degenerate() -> dict:
    """V1（h_m=0）、V2（T_a≡T_0、C_a≡C_0）、V4（R≡R_0 退化）、V5（D 与 C_a 恒定）。"""
    out = {}

    # V1：关掉传质边界通量 → 水分场必须全程冻结在初值，且守恒残差为 0。
    v1 = _run("2026A/verify/v1/hm0",
              {"appendix": APPENDIX_OF_PROBLEM[3], "N": N_GRID, "out_dt": P3_OUT_DT_S,
               "t_end": ATT1_T_END_S, "dt1": DT_COARSE_STAGE1, "dt2": DT_STAGE2,
               "scale_hm": 0.0, "record": "last", "tag": "V1-hm0"},
              want_tdry=False, extra=_field_stats)
    out["V1"] = dict(v1, setting="scale_hm=0", expect="C≡C_init 且 ε_mass<1e-12")

    # V2：把空气拉回药材初态 → 两场都是不动点（驱动力为零，源项不得有杂散贡献）。
    v2 = _run("2026A/verify/v2/fixedpoint",
              {"appendix": APPENDIX_OF_PROBLEM[3], "N": N_GRID, "out_dt": P3_OUT_DT_S,
               "t_end": ATT1_T_END_S, "dt1": DT_COARSE_STAGE1, "dt2": DT_STAGE2,
               "const_Ta": T_INIT, "const_Ca": C_INIT, "record": "last",
               "tag": "V2-fixedpoint"},
              want_tdry=False, extra=_field_stats)
    out["V2"] = dict(v2, setting=f"const_Ta={T_INIT}, const_Ca={C_INIT}",
                     expect="T≡T_init 且 C≡C_init（不动点）")

    # V4：问题 4 求解器令 R≡R_0 后，物质坐标与实验室坐标必须逐位一致（Ṙ=0 时对流项消失）。
    common = {"appendix": APPENDIX_OF_PROBLEM[4], "N": N_GRID, "out_dt": P4_OUT_DT_S,
              "t_end": ATT1_T_END_S, "dt1": DT_COARSE_STAGE1, "dt2": DT_STAGE2,
              "shrink": False, "record": "last"}
    v4m = _run("2026A/verify/v4/material",
               dict(common, frame="material", tag="V4-material"), want_tdry=False)
    v4l = _run("2026A/verify/v4/lab", dict(common, frame="lab", tag="V4-lab"),
               want_tdry=False)
    keys = ("T_center", "T_surface", "C_center", "C_surface", "C_bar", "eps_mass")
    out["V4"] = {"setting": "shrink=False，material vs lab", "material": v4m, "lab": v4l,
                 "max_abs_diff": {k: abs(v4m[k] - v4l[k]) for k in keys},
                 "expect": "Ṙ=0 时两坐标系逐位一致（对流项恒为 0）"}

    # V5：冻结物性（D 常数）+ 恒定 C_a → C 单调降至 C_a，无过冲。
    v5 = _run("2026A/verify/v5/equilibrium",
              {"appendix": APPENDIX_OF_PROBLEM[3], "N": N_GRID, "out_dt": DT_STAGE3,
               "t_end": T_MAX_SEARCH_S, "dt1": DT_COARSE_STAGE1, "dt2": DT_STAGE2,
               "dt3": DT_STAGE3, "frozen_state": [C_INIT, T_INIT],
               "const_Ca": PLATEAU_CA, "record": "last", "tag": "V5-equilibrium"},
              want_tdry=False, extra=_field_stats)
    d_frozen = float(Props(APPENDIX_OF_PROBLEM[3],
                           frozen_state=(C_INIT, T_INIT)).D(np.zeros(1), np.zeros(1))[0])
    out["V5"] = dict(v5, setting=f"frozen_state=({C_INIT},{T_INIT}), const_Ca={PLATEAU_CA}",
                     D_frozen_m2_s=d_frozen, C_a=float(PLATEAU_CA),
                     residual_frac=abs(v5["C_bar"] - PLATEAU_CA) / (C_INIT - PLATEAU_CA),
                     expect="C̄ 单调降至 C_a，无过冲（C_min ≥ C_a）")
    return out


def stage_v6() -> dict:
    """V6：N 加倍且 Δt 同时减半，关键输出变化必须 < CONV_TOL_PCT。"""
    base_q1 = _q1(N_GRID, DT_STAGE2)
    fine_q1 = _q1(2 * N_GRID, DT_STAGE2 / 2.0)
    base_q3 = _q3(N_GRID, DT_STAGE3)
    fine_q3 = _q3(2 * N_GRID, DT_STAGE3 / 2.0)
    return {"q1": {"base": base_q1, "fine": fine_q1, "key_output": "C_surface",
                   "delta_pct": _rel_pct(fine_q1["C_surface"], base_q1["C_surface"])},
            "q3": {"base": base_q3, "fine": fine_q3, "key_output": "t_dry_h",
                   "delta_pct": _rel_pct(fine_q3["t_dry_h"], base_q3["t_dry_h"])}}


# ---------------------------------------------------------------- 汇总与判据

def _refine_table(cases: list, x_key: str, y_key: str) -> list:
    """相邻加密档之间的相对变化（%）。"""
    out = []
    for a, b in zip(cases[:-1], cases[1:]):
        out.append({"from": a[x_key], "to": b[x_key], "y_from": a[y_key], "y_to": b[y_key],
                    "delta_pct": _rel_pct(b[y_key], a[y_key])})
    return out


def _observed_order(y: list, h: list) -> float:
    """三点等比加密的观测收敛阶 p = log((y1-y2)/(y2-y3)) / log(h1/h2)。"""
    d1, d2 = y[0] - y[1], y[1] - y[2]
    if abs(d2) < 1e-15:
        return float("inf")
    return float(np.log(abs(d1 / d2)) / np.log(h[0] / h[1]))


def assemble(parts: dict) -> dict:
    """把各阶段结果拼成 figures/verification_results.json 的内容。"""
    q1g, q3g, q4g = parts["q1grid"], parts["q3grid"], parts["q4grid"]
    q1t, q3t = parts["q1time"], parts["q3time"]
    bes, deg, v6 = parts["bessel"], parts["degenerate"], parts["v6"]

    grid = {
        "tol_pct": CONV_TOL_PCT,
        "problem1": {"fixed_dt_s": q1g["fixed_dt_s"], "t_eval_s": q1g["t_eval_s"],
                     "cases": q1g["cases"], "key_output": "C_surface",
                     "refinement": _refine_table(q1g["cases"], "n_grid", "C_surface")},
        "problem3": {"cases": q3g["cases"], "key_output": "t_dry_h",
                     "refinement": _refine_table(q3g["cases"], "n_grid", "t_dry_h")},
        "problem4": {"cases": q4g["cases"], "key_output": "t_dry_h",
                     "refinement": _refine_table(q4g["cases"], "n_grid", "t_dry_h")},
    }
    # 生产网格 N=40 加倍到 N=80 是 PA-C1 判据的落点；N=20→40 超差正是不取 N=20 的理由。
    prod = {}
    for tag, blk, key in (("problem1", grid["problem1"], "C_surface"),
                          ("problem3", grid["problem3"], "t_dry_h"),
                          ("problem4", grid["problem4"], "t_dry_h")):
        by_n = {c["n_grid"]: c[key] for c in blk["cases"]}
        prod[tag] = {"key_output": key, "N_prod": N_GRID, "N_double": 2 * N_GRID,
                     "y_prod": by_n[N_GRID], "y_double": by_n[2 * N_GRID],
                     "delta_pct": _rel_pct(by_n[2 * N_GRID], by_n[N_GRID])}
        if min(GRID_CONV) in by_n:
            prod[tag]["delta_pct_half_to_prod"] = _rel_pct(by_n[N_GRID], by_n[min(GRID_CONV)])
    grid["production_doubling"] = prod

    t_c = [c["T_center"] for c in q1t["cases"]]
    dts = [c["dt_s"] for c in q1t["cases"]]
    fit = np.polyfit(dts, t_c, 1)
    resid = np.asarray(t_c) - np.polyval(fit, dts)
    ss_tot = float(np.sum((np.asarray(t_c) - np.mean(t_c)) ** 2))
    time = {
        "tol_pct": CONV_TOL_PCT,
        "problem1_first_order": {
            "fixed_N": q1t["fixed_N"], "dt_s": dts, "T_center": t_c,
            "C_surface": [c["C_surface"] for c in q1t["cases"]],
            "deltas": [t_c[i] - t_c[i + 1] for i in range(len(t_c) - 1)],
            "observed_order_halving": _observed_order(t_c[1:4], dts[1:3]),
            "slope_per_s": float(fit[0]), "richardson_limit": float(fit[1]),
            "fit_r2": float(1.0 - float(np.sum(resid ** 2)) / ss_tot),
            "richardson_two_finest": float(
                t_c[-1] + (t_c[-1] - t_c[-2]) * dts[-1] / (dts[-2] - dts[-1])),
        },
        "problem3_plateau_dt": {
            "fixed_N": q3t["fixed_N"], "cases": q3t["cases"],
            "refinement": _refine_table(q3t["cases"], "dt3_s", "t_dry_h"),
            "prod_vs_finest_pct": _rel_pct(
                next(c["t_dry_h"] for c in q3t["cases"] if c["dt3_s"] == min(DT_CONV)),
                next(c["t_dry_h"] for c in q3t["cases"] if c["dt3_s"] == DT_STAGE3)),
        },
    }

    all_cases = (q1g["cases"] + q3g["cases"] + q4g["cases"] + q1t["cases"] + q3t["cases"]
                 + bes["cases"] + bes["grid_independence"]
                 + [deg["V1"], deg["V2"], deg["V4"]["material"], deg["V4"]["lab"], deg["V5"]]
                 + [v6["q1"]["fine"], v6["q3"]["fine"]])
    eps = [c["eps_mass"] for c in all_cases]
    mass = {"n_cases": len(eps), "eps_max": float(max(eps)), "eps_min": float(min(eps)),
            "tol": 1e-2, "identity": "|(m0 - m_end) - ∫ h_m/R (C_N - C_a) dt| / m0",
            "note": "控制体体积精确积分 + 界面通量成对相消 \u21d2 守恒是离散恒等式，故残差在机器精度量级"}

    errs = [c["err_center"] for c in bes["cases"]]
    bdt = [c["dt_s"] for c in bes["cases"]]
    analytic = {
        **bes["analytic"], "cases": bes["cases"],
        "order_between_halvings": [
            float(np.log(errs[i] / errs[i + 1]) / np.log(bdt[i] / bdt[i + 1]))
            for i in range(len(errs) - 1)],
        "grid_independence": [{"n_grid": c["n_grid"], "err_center": c["err_center"]}
                              for c in bes["grid_independence"]],
        "grid_err_ratio": (max(c["err_center"] for c in bes["grid_independence"])
                           / min(c["err_center"] for c in bes["grid_independence"])),
        "analytic_at_rounded_Ta": bes["analytic_at_rounded_Ta"],
    }
    return {"problem": "all", "tag": "PA-C1", "grid_convergence": grid,
            "time_convergence": time, "mass_conservation": mass,
            "analytic_cross_check": analytic, "degenerate_cases": deg,
            "simultaneous_refinement": v6}


def validate_constraints(data: dict) -> bool:
    """K12 与 V1—V6：全部从写出的结果复核，任一条不过即硬停。"""
    grid, time = data["grid_convergence"], data["time_convergence"]
    mass, ana = data["mass_conservation"], data["analytic_cross_check"]
    deg, v6 = data["degenerate_cases"], data["simultaneous_refinement"]

    # ---- K12 / PA-C1 前半：生产网格加倍后关键输出变化 < 0.1%。
    for tag, blk in grid["production_doubling"].items():
        assert blk["N_prod"] == N_GRID and blk["N_double"] == 2 * N_GRID, "加倍档位不对"
        assert blk["delta_pct"] < CONV_TOL_PCT, \
            f"K12 {tag} 网格加倍变化 {blk['delta_pct']:.4f}% 超 {CONV_TOL_PCT}%"
    # 反向登记：N=20→40 超差，正是不取 N=20 的量化理由（若它也达标则报告 10.1 的取舍失据）。
    half = grid["production_doubling"]["problem1"]["delta_pct_half_to_prod"]
    assert half > CONV_TOL_PCT, \
        f"问题 1 的 N=20→40 变化 {half:.4f}% 未超 {CONV_TOL_PCT}%，报告 10.1 取 N=40 的理由失据"

    # 网格加密必须单调收敛：相邻档变化量逐次减小。
    for tag in ("problem1", "problem3", "problem4"):
        d = [r["delta_pct"] for r in grid[tag]["refinement"]]
        assert all(d[i] > d[i + 1] for i in range(len(d) - 1)), \
            f"{tag} 网格加密变化量未逐次减小：{d}"

    # ---- PA-C1 后半：每个算例的质量守恒残差 < 1%。
    assert mass["eps_max"] < mass["tol"], f"PA-C1 ε_mass 最大 {mass['eps_max']:.3e} 超 1%"
    assert mass["n_cases"] >= 20, f"验证算例只有 {mass['n_cases']} 个，证据不足"

    # ---- 10.2 时间步：后向 Euler 必须表现一阶精度。
    fo = time["problem1_first_order"]
    assert fo["fit_r2"] > 0.999, f"T(Δt) 对 Δt 的线性拟合 R²={fo['fit_r2']:.6f}，非一阶"
    assert 0.8 < fo["observed_order_halving"] < 1.3, \
        f"观测收敛阶 {fo['observed_order_halving']:.3f} 不在一阶邻域"
    dd = [abs(v) for v in fo["deltas"]]
    assert all(dd[i] > dd[i + 1] for i in range(len(dd) - 1)), f"逐次加密变化量未减小：{dd}"
    assert time["problem3_plateau_dt"]["prod_vs_finest_pct"] < CONV_TOL_PCT, \
        f"问题 3 生产步长相对最细步长偏差 " \
        f"{time['problem3_plateau_dt']['prod_vs_finest_pct']:.4f}% 超 {CONV_TOL_PCT}%"

    # ---- 10.4 / V3：对 Bessel 解析解的偏差随 Δt 一阶下降，且不随空间加密改变。
    e = [c["err_center"] for c in ana["cases"]]
    assert all(e[i] > e[i + 1] for i in range(len(e) - 1)), f"解析解偏差未随 Δt 单调下降：{e}"
    for p in ana["order_between_halvings"]:
        assert 0.7 < p < 1.3, f"解析解偏差的观测阶 {p:.3f} 不在一阶邻域"
    assert ana["grid_err_ratio"] < 3.0, \
        f"固定 Δt 下偏差随 N 变化 {ana['grid_err_ratio']:.2f} 倍，空间离散有系统性偏差"
    assert 1.0 < ana["Bi_h"] < 2.0 and ana["zeta_first3"][0] > 1.0, "Bi_h 或特征根异常"

    # ---- V1：h_m=0 => 水分场冻结、守恒残差为 0。
    v1 = deg["V1"]
    assert v1["C_dev_from_init"] < 1e-12, f"V1 水分场偏离初值 {v1['C_dev_from_init']:.3e}"
    assert v1["eps_mass"] < 1e-12, f"V1 守恒残差 {v1['eps_mass']:.3e} 未达机器精度"
    assert v1["T_dev_from_init"] > 1.0, "V1 温度场必须仍在升温（只关传质，不关传热）"

    # ---- V2：空气回到药材初态 => 两场都是不动点。
    v2 = deg["V2"]
    assert v2["C_dev_from_init"] < 1e-12, f"V2 水分场离开不动点 {v2['C_dev_from_init']:.3e}"
    assert v2["T_dev_from_init"] < 1e-12, f"V2 温度场离开不动点 {v2['T_dev_from_init']:.3e}"

    # ---- V4：R≡R_0 时两坐标系必须逐位一致（Ṙ=0 => 对流项恒为 0）。
    for k, v in deg["V4"]["max_abs_diff"].items():
        assert v < 1e-14, f"V4 {k} 两坐标系差 {v:.3e}，坐标变换引入了多余项"

    # ---- V5：D 与 C_a 恒定 => C 单调降至 C_a，无过冲。
    v5 = deg["V5"]
    assert v5["C_bar_max_rise"] <= 1e-12, f"V5 平均含水率回升 {v5['C_bar_max_rise']:.3e}"
    assert v5["C_min"] >= v5["C_a"] - 1e-9, f"V5 出现过冲：C_min={v5['C_min']:.6e} < C_a"
    assert v5["residual_frac"] < 1e-2, f"V5 未趋于平衡，残差占初始驱动力 {v5['residual_frac']:.3e}"
    assert v5["C_spread"] < 1e-3, f"V5 平衡态仍有径向梯度 {v5['C_spread']:.3e}"

    # ---- V6：N 加倍且 Δt 减半，关键输出变化 < 0.1%。
    for tag in ("q1", "q3"):
        assert v6[tag]["delta_pct"] < CONV_TOL_PCT, \
            f"V6 {tag} 同时加密变化 {v6[tag]['delta_pct']:.4f}% 超 {CONV_TOL_PCT}%"
    return True


def validate_capability(data: dict) -> dict:
    """PA-C1 的可反驳检验：收敛性证据 + 质量守恒，两条都必须由实算数字支撑。"""
    grid = data["grid_convergence"]["production_doubling"]
    mass = data["mass_conservation"]
    time = data["time_convergence"]
    ana = data["analytic_cross_check"]
    v6 = data["simultaneous_refinement"]

    worst_grid = max(grid.values(), key=lambda b: b["delta_pct"])
    assert worst_grid["delta_pct"] < CONV_TOL_PCT, "PA-C1 网格加倍判据不成立"
    assert mass["eps_max"] < mass["tol"], "PA-C1 质量守恒判据不成立"
    fo = time["problem1_first_order"]
    q1 = grid["problem1"]
    return {"PA-C1": {
        "status": "PASS",
        "evidence": (
            f"网格加倍：以生产网格 N={N_GRID} 为基准加倍到 N={2 * N_GRID}，"
            f"问题 1 表面含水率 {q1['y_prod']:.6f}→{q1['y_double']:.6f}（{q1['delta_pct']:.4f}%）、"
            f"问题 3 t_dry {grid['problem3']['y_prod']:.4f}→{grid['problem3']['y_double']:.4f} h"
            f"（{grid['problem3']['delta_pct']:.4f}%）、"
            f"问题 4 t_dry {grid['problem4']['y_prod']:.4f}→{grid['problem4']['y_double']:.4f} h"
            f"（{grid['problem4']['delta_pct']:.4f}%），最大 {worst_grid['delta_pct']:.4f}% "
            f"< {CONV_TOL_PCT}%；N={min(GRID_CONV)}→{N_GRID} 为 "
            f"{q1['delta_pct_half_to_prod']:.4f}%（超差，故生产网格不取 {min(GRID_CONV)}）。"
            f"时间步：问题 1 在 Δt={max(P1_DT_CONV):g}→{min(P1_DT_CONV):g} s 上 T(0) 对 Δt 线性"
            f"（R²={fo['fit_r2']:.6f}，观测阶 {fo['observed_order_halving']:.3f}），"
            f"Richardson 极限 {fo['richardson_limit']:.4f} °C；问题 3 平台段 Δt3="
            f"{DT_STAGE3:g} s 相对 {min(DT_CONV):g} s 偏差 "
            f"{time['problem3_plateau_dt']['prod_vs_finest_pct']:.4f}%。"
            f"同时加密（N 加倍且 Δt 减半）变化 q1 {v6['q1']['delta_pct']:.4f}%、"
            f"q3 {v6['q3']['delta_pct']:.4f}%。"
            f"质量守恒：{mass['n_cases']} 个算例的 ε_mass ∈ "
            f"[{mass['eps_min']:.3e}, {mass['eps_max']:.3e}]，全部 < 1%。"
            f"解析解交叉校验：Bi_h={ana['Bi_h']:.4f}、ζ₁₋₃="
            f"{', '.join(f'{v:.4f}' for v in ana['zeta_first3'])}，"
            f"T(0) 解析 {ana['T_center_analytic']:.4f} °C，数值偏差随 Δt 从 "
            f"{ana['cases'][0]['err_center']:.2e} 降至 {ana['cases'][-1]['err_center']:.2e} °C"
            f"（观测阶 {', '.join(f'{p:.2f}' for p in ana['order_between_halvings'])}），"
            f"固定 Δt 下随 N 仅变化 {ana['grid_err_ratio']:.2f} 倍。"),
        "grid_delta_pct_max": worst_grid["delta_pct"], "eps_mass_max": mass["eps_max"],
    }}


def _report(data: dict) -> None:
    grid, time = data["grid_convergence"], data["time_convergence"]
    ana, deg, v6 = data["analytic_cross_check"], data["degenerate_cases"], data["simultaneous_refinement"]
    print("[VER] 10.1 网格收敛（生产网格加倍）")
    for tag, blk in grid["production_doubling"].items():
        extra = (f"  (N={min(GRID_CONV)}→{N_GRID}: "
                 f"{blk['delta_pct_half_to_prod']:.4f}%)" if "delta_pct_half_to_prod" in blk else "")
        print(f"      {tag}: {blk['key_output']} {blk['y_prod']:.6f}→{blk['y_double']:.6f}"
              f"  Δ={blk['delta_pct']:.4f}%{extra}")
    fo = time["problem1_first_order"]
    print(f"[VER] 10.2 时间步：问题 1 T(0)="
          f"{'/'.join(f'{v:.4f}' for v in fo['T_center'])}"
          f"  观测阶={fo['observed_order_halving']:.3f}  R²={fo['fit_r2']:.6f}"
          f"  Richardson={fo['richardson_limit']:.4f}")
    q3dt = time["problem3_plateau_dt"]
    print(f"      问题 3 t_dry="
          f"{'/'.join(format(c['t_dry_h'], '.4f') for c in q3dt['cases'])} h"
          f"  {DT_STAGE3:g}s vs {min(DT_CONV):g}s = {q3dt['prod_vs_finest_pct']:.4f}%")
    m = data["mass_conservation"]
    print(f"[VER] 10.3 质量守恒：{m['n_cases']} 算例 ε_mass ∈ "
          f"[{m['eps_min']:.3e},{m['eps_max']:.3e}]")
    print(f"[VER] 10.4 解析解：Ta={ana['Ta_const']:.4f} Bi_h={ana['Bi_h']:.4f} "
          f"T(0)解析={ana['T_center_analytic']:.4f} T(R0)解析={ana['T_surface_analytic']:.4f}")
    for c in ana["cases"]:
        print(f"      Δt={c['dt_s']:>5g}s  T(0)={c['T_center']:.4f} 偏差={c['err_center']:.2e}"
              f"  T(R0)={c['T_surface']:.4f} 偏差={c['err_surface']:.2e}")
    print(f"      观测阶={', '.join(f'{p:.3f}' for p in ana['order_between_halvings'])}"
          f"  固定Δt随N变化 {ana['grid_err_ratio']:.2f} 倍")
    print(f"[VER] V1 ΔC={deg['V1']['C_dev_from_init']:.3e} ΔT={deg['V1']['T_dev_from_init']:.4f} "
          f"ε={deg['V1']['eps_mass']:.3e}")
    print(f"[VER] V2 ΔC={deg['V2']['C_dev_from_init']:.3e} ΔT={deg['V2']['T_dev_from_init']:.3e}")
    print(f"[VER] V4 material vs lab 最大差="
          f"{max(deg['V4']['max_abs_diff'].values()):.3e}")
    print(f"[VER] V5 C̄={deg['V5']['C_bar']:.6e} → C_a={deg['V5']['C_a']:.6e}"
          f"  残差占比={deg['V5']['residual_frac']:.3e}  D冻结={deg['V5']['D_frozen_m2_s']:.4e}"
          f"  回升={deg['V5']['C_bar_max_rise']:.3e}")
    print(f"[VER] V6 q1 Δ={v6['q1']['delta_pct']:.4f}%  q3 Δ={v6['q3']['delta_pct']:.4f}%")


def main() -> int:
    ensure_dirs()
    args = sys.argv[1:]
    want = None
    if "--stage" in args:
        want = args[args.index("--stage") + 1]
    runners = {name: globals()[f"stage_{name}"] for name in STAGES}

    if want is not None and want != "assemble":
        assert want in runners, f"未知阶段 {want}，可选 {STAGES + ['assemble']}"
        print(f"[VER] 阶段 {want}")
        runners[want]()
        print(f"[VER] 阶段 {want} 完成，检查点已落盘")
        return 0

    parts = {}
    for name in STAGES:
        print(f"[VER] 阶段 {name}")
        parts[name] = runners[name]()
    data = assemble(parts)
    assert validate_constraints(data) is True, "验证阶段约束未通过"
    cap = validate_capability(data)
    data["capability_verdict"] = cap
    _report(data)
    for cid, v in cap.items():
        print(f"[VER][{cid}] {v['status']}")
    save_json(FIG_DIR / "verification_results.json", jsonable(data))
    print(f"[VER] 写出 figures/verification_results.json  "
          f"{(FIG_DIR / 'verification_results.json').stat().st_size / 1024:.1f} KB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
