# -*- coding: utf-8 -*-
"""模型对比与择优（MODELING_REPORT 第十一节、能力 PA-C2）。

四个候选模型的每一项指标都由本文代码实算，无文献值、无定性占位：
  M1 集中参数（薄层干燥）：dC̄/dt = -(2/R_0) h_m (C̄ - C_a(t))，本文件内独立积分；
  M2 常物性一维：物性冻结在初始态 (C_0, T_0)，附录 3 公式求值一次后固定；
  M3 变物性一维：问题 3 的生产算例（读 figures/problem_3_results.json）；
  M4 变物性 + 动边界：问题 4 的物质坐标主算例（读 figures/problem_4_results.json）。
判据统一为 max_r C < 0.15；M1 只有均值，其判据退化为 C̄ < 0.15——这正是它的缺陷。

运行：python _utils/run_compute.py code/model_compare.py
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_utils"))
from compute_checkpoint import CaseCheckpoint  # noqa: E402

from params import *  # noqa: F401,F403,E402
from utils import (Props, downsample_index, ensure_dirs, jsonable, make_ambient,  # noqa: E402
                   node_columns, round_out, run_case, save_json)

ROOT = Path(__file__).resolve().parents[1]
LUMPED_TAU_S = R0 / (2.0 * H_MASS)      # 集中参数时间常数 R_0/(2 h_m)
BI_LUMPED_MAX = 0.1                     # 集中参数适用上限（Bi < 0.1）
FROZEN_STATE = [C_INIT, T_INIT]         # M2 唯一在求解前已知的状态


def solve_lumped() -> dict:
    """M1：一维标量 ODE，后向 Euler + 与 PDE 相同的 60 s 输出与跨阈值插值。

    dC̄/dt = -(2/R_0) h_m (C̄ - C_a(t))，解析地看是一阶松弛，时间常数 R_0/(2 h_m)。
    """
    _, c_amb = make_ambient()
    k = 2.0 * H_MASS / R0
    dt = DT_STAGE2
    n_out = int(round(T_MAX_SEARCH_S / P3_OUT_DT_S))
    c_bar = C_INIT
    t_cur = 0.0
    rec_t, rec_c = [0.0], [C_INIT]
    t_dry, bracket = None, None
    for j in range(1, n_out + 1):
        t_target = j * P3_OUT_DT_S
        n_sub = max(1, int(round((t_target - t_cur) / dt)))
        h = (t_target - t_cur) / n_sub
        for _ in range(n_sub):
            t_new = t_cur + h
            # 后向 Euler：(1 + k h) C^{n+1} = C^n + k h C_a(t^{n+1})
            c_bar = (c_bar + k * h * c_amb(t_new)) / (1.0 + k * h)
            t_cur = t_new
        t_cur = float(t_target)
        rec_t.append(t_cur)
        rec_c.append(float(c_bar))
        if c_bar < C_THRESHOLD:
            t_prev, c_prev = rec_t[-2], rec_c[-2]
            assert c_prev >= C_THRESHOLD > c_bar, "M1 跨阈值区间定位失败"
            t_dry = t_prev + (c_prev - C_THRESHOLD) / (c_prev - c_bar) * (t_cur - t_prev)
            bracket = (t_prev, t_cur, c_prev, float(c_bar))
            break
    assert t_dry is not None, "M1 未在搜索窗内达标"
    # 恒定 C_a 的闭式解，用于独立核对数值积分（C_a 实测随时间变化，故只作量级校验）
    ca_plateau = PLATEAU_CA
    t_closed = -LUMPED_TAU_S * np.log((C_THRESHOLD - ca_plateau) / (C_INIT - ca_plateau))
    return {"model": "M1", "name": "集中参数（薄层干燥）",
            "ode": "dC_bar/dt = -(2/R0) h_m (C_bar - C_a(t))",
            "tau_s": LUMPED_TAU_S, "tau_h": LUMPED_TAU_S / SEC_PER_HOUR,
            "t_dry_s": float(t_dry), "t_dry_h": float(t_dry / SEC_PER_HOUR),
            "bracket": [float(v) for v in bracket],
            "criterion_applied": "C_bar < 0.15（无径向分布，逐点判据无法定义）",
            "t_dry_h_closed_form_const_Ca": float(t_closed / SEC_PER_HOUR),
            "radial_capability": False, "C_center_minus_surface": 0.0,
            "uses_D": False, "can_deliver_pointwise": False,
            "series": {"times_h": [float(v / SEC_PER_HOUR) for v in rec_t],
                       "C_bar": [float(v) for v in rec_c]}}


def solve_frozen(n_grid: int) -> dict:
    """M2：附录 3 公式在 (C_0, T_0) 求值一次后冻结，其余与 M3 完全相同。"""
    cfg = {"appendix": APPENDIX_OF_PROBLEM[3], "N": n_grid, "out_dt": P3_OUT_DT_S,
           "t_end": T_MAX_SEARCH_S, "stop_threshold": C_THRESHOLD, "dt1": DT_STAGE1,
           "dt2": DT_STAGE2, "dt3": DT_STAGE3, "frozen_state": FROZEN_STATE,
           "record": "full", "tag": f"M2-frozen-N{n_grid}"}
    res = run_case(cfg)
    assert res["t_dry"] is not None, "M2 未在搜索窗内达标"
    c_field = np.asarray(res["C_field"], dtype=float)
    t_field = np.asarray(res["T_field"], dtype=float)
    c_cols = node_columns(c_field, n_grid)
    t_out = np.asarray(res["t"], dtype=float)
    dsi = downsample_index(len(t_out), 400)
    return {"N": n_grid, "t_dry_s": float(res["t_dry"]),
            "t_dry_h": float(res["t_dry"] / SEC_PER_HOUR),
            "bracket": [float(v) for v in res["bracket"]],
            "C_center_end": float(c_field[-1, 0]), "C_surface_end": float(c_field[-1, -1]),
            "C_center_minus_surface": float(c_field[-1, 0] - c_field[-1, -1]),
            "C_bar_end": float(res["C_bar"][-1]), "eps_mass": float(res["eps_mass"]),
            "n_steps": int(res["n_steps"]), "wall_s": float(res["wall_s"]),
            "D_frozen": float(res["D_surface"][0]),
            "profile_end_cm": list(OUT_RADII_CM),
            "profile_end_C": [round_out(v) for v in c_cols[-1]],
            "series": {"times_h": [float(t_out[i] / SEC_PER_HOUR) for i in dsi],
                       "C_max": [float(res["C_max"][i]) for i in dsi],
                       "C_bar": [float(res["C_bar"][i]) for i in dsi],
                       "T_center": [float(t_field[i, 0]) for i in dsi]}}


def _load(problem: int) -> dict:
    path = FIG_DIR / f"problem_{problem}_results.json"
    assert path.exists(), f"缺 {path.name}，先跑 code/problem{problem}.py"
    return json.loads(path.read_text(encoding="utf-8"))


def biot_evidence() -> dict:
    """M1/M2/M3 的适用性判据实算：Bi_m 全程量程、Bi_h 量程、D 的跨度。"""
    p3 = _load(3)
    cs = p3["criterion_series"]
    props = Props(APPENDIX_OF_PROBLEM[3])
    bi_m = np.asarray(cs["Bi_m_surface"], dtype=float)
    c_surf = np.asarray(cs["C_surface"], dtype=float)
    c_cent = np.asarray(cs["C_max"], dtype=float)
    bi_h = props.bi_heat(np.concatenate([c_surf, c_cent]))
    d_all = np.asarray(cs["D_center"] + cs["D_surface"], dtype=float)

    # D 在四个「合理」参考状态上的取值：冻结状态的选择本身不自洽（报告 11 节 M2 段）。
    ref = {"(C0,T0)=(2.55,28)": (C_INIT, T_INIT),
           "(C0,Ta)=(2.55,50)": (C_INIT, PLATEAU_TA),
           "(Cth,Ta)=(0.15,50)": (C_THRESHOLD, PLATEAU_TA),
           "(Csurf,Ta)=(0.0526,50)": (float(c_surf[-1]), PLATEAU_TA)}
    d_ref = {k: float(props.D(np.array([c]), np.array([t]))[0]) for k, (c, t) in ref.items()}
    # 全程 Bi_m 的下界出现在 D 最大处（表面被加热而尚未变干），单独实算以免误报量程
    d_max_hot = float(props.D(np.array([C_INIT]), np.array([PLATEAU_TA]))[0])
    return {"Bi_m_surface_min": float(bi_m.min()), "Bi_m_surface_max": float(bi_m.max()),
            "Bi_m_at_D_max": float(H_MASS * R0 / d_max_hot),
            "Bi_m_lumped_limit": BI_LUMPED_MAX,
            "Bi_h_min": float(bi_h.min()), "Bi_h_max": float(bi_h.max()),
            "D_min": float(d_all.min()), "D_max": float(d_all.max()),
            "D_span_orders": float(np.log10(d_all.max() / d_all.min())),
            "D_at_reference_states": d_ref,
            "D_ref_spread_ratio": max(d_ref.values()) / min(d_ref.values()),
            "source": "problem_3_results.json 的 criterion_series（M3 生产算例实测）"}


def m3_indicators() -> dict:
    """M3：读问题 3 的生产算例，不重算（同一求解器、同一设置）。"""
    p3 = _load(3)
    td, diag, cs = p3["t_dry"], p3["diagnostics"], p3["criterion_series"]
    grid = p3["result_grid"]
    c_last = np.asarray(grid["C"][-1], dtype=float)
    dsi = downsample_index(len(cs["times_h"]), 400)
    return {"t_dry_h": float(td["hours"]),
            "C_center_end": float(c_last[0]), "C_surface_end": float(c_last[-1]),
            "C_center_minus_surface": float(c_last[0] - c_last[-1]),
            "C_bar_end": float(td["C_bar_at_t_dry"]), "eps_mass": float(diag["eps_mass"]),
            "n_steps": int(p3["schedule"]["n_steps"]), "wall_s": float(diag["wall_s"]),
            "profile_end_cm": list(grid["radii_cm"]),
            "profile_end_C": [round_out(v) for v in c_last],
            "series": {"times_h": [float(cs["times_h"][i]) for i in dsi],
                       "C_max": [float(cs["C_max"][i]) for i in dsi],
                       "C_bar": [float(cs["C_bar"][i]) for i in dsi],
                       "T_center": [float(cs["T_center"][i]) for i in dsi]},
            "source": "figures/problem_3_results.json（问题 3 生产算例）"}


def m4_indicators() -> dict:
    """M4：读问题 4 的物质坐标主算例，不重算。"""
    p4 = _load(4)
    main, td = p4["main"], p4["main"]["t_dry"]
    diag, cs = main["diagnostics"], main["criterion_series"]
    dsi = downsample_index(len(cs["times_h"]), 400)
    return {"t_dry_h": float(td["hours"]),
            "C_center_end": float(td["C_center_at_t_dry"]),
            "C_surface_end": float(td["C_surface_at_t_dry"]),
            "C_center_minus_surface": float(td["C_center_at_t_dry"] - td["C_surface_at_t_dry"]),
            "C_bar_end": float(td["C_bar_at_t_dry"]), "eps_mass": float(diag["eps_mass"]),
            "n_steps": int(diag["n_steps"]), "wall_s": float(diag["wall_s"]),
            "R_at_t_dry_cm": float(td["R_at_t_dry_cm"]), "R_ratio": float(td["R_ratio"]),
            "series": {"times_h": [float(cs["times_h"][i]) for i in dsi],
                       "C_max": [float(cs["C_max"][i]) for i in dsi],
                       "C_bar": [float(cs["C_bar"][i]) for i in dsi],
                       "T_center": [float(cs["T_center"][i]) for i in dsi],
                       "R_cm": [float(cs["R_cm"][i]) for i in dsi]},
            "source": "figures/problem_4_results.json 的 main 分支（物质坐标）"}


def assemble(m1: dict, m2: dict, m3: dict, m4: dict, biot: dict) -> dict:
    """五指标对比表：t_dry、径向分布能力、中心-表面落差、质量守恒、可交付性。"""
    ref_h = m3["t_dry_h"]
    rows = [
        {"model": "M1", "name": "集中参数（薄层干燥）", "dim": "0 维（ODE）",
         "properties": "常物性（不进入方程）", "domain": "固定",
         "t_dry_h": m1["t_dry_h"], "radial_capability": False,
         "C_center_minus_surface": 0.0, "eps_mass": None,
         "criterion_applied": m1["criterion_applied"], "can_deliver_pointwise": False,
         "rel_vs_M3_pct": (m1["t_dry_h"] - ref_h) / ref_h * 100.0,
         "n_steps": None, "wall_s": None,
         "evidence": "本文件 solve_lumped() 实算（后向 Euler，输出步与 PDE 一致）"},
        {"model": "M2", "name": "常物性一维（物性冻结于初始态）", "dim": "1 维径向",
         "properties": f"冻结在 (C,T)=({FROZEN_STATE[0]}, {FROZEN_STATE[1]})", "domain": "固定",
         "t_dry_h": m2["t_dry_h"], "radial_capability": True,
         "C_center_minus_surface": m2["C_center_minus_surface"], "eps_mass": m2["eps_mass"],
         "criterion_applied": "max_r C < 0.15", "can_deliver_pointwise": True,
         "rel_vs_M3_pct": (m2["t_dry_h"] - ref_h) / ref_h * 100.0,
         "n_steps": m2["n_steps"], "wall_s": m2["wall_s"],
         "evidence": "本文件 solve_frozen() 实算（frozen_state 传入同一求解器）"},
        {"model": "M3", "name": "变物性一维（问题 3 主模型）", "dim": "1 维径向",
         "properties": "附录 3 随 C、T 逐层重算", "domain": "固定",
         "t_dry_h": m3["t_dry_h"], "radial_capability": True,
         "C_center_minus_surface": m3["C_center_minus_surface"], "eps_mass": m3["eps_mass"],
         "criterion_applied": "max_r C < 0.15", "can_deliver_pointwise": True,
         "rel_vs_M3_pct": 0.0, "n_steps": m3["n_steps"], "wall_s": m3["wall_s"],
         "evidence": m3["source"]},
        {"model": "M4", "name": "变物性 + 动边界（问题 4 主模型）", "dim": "1 维径向 + 动边界",
         "properties": "附录 4 随 C、T 逐层重算", "domain": f"收缩 R/R0={m4['R_ratio']:.4f}",
         "t_dry_h": m4["t_dry_h"], "radial_capability": True,
         "C_center_minus_surface": m4["C_center_minus_surface"], "eps_mass": m4["eps_mass"],
         "criterion_applied": "max_r C < 0.15", "can_deliver_pointwise": True,
         "rel_vs_M3_pct": (m4["t_dry_h"] - ref_h) / ref_h * 100.0,
         "n_steps": m4["n_steps"], "wall_s": m4["wall_s"],
         "evidence": m4["source"]},
    ]
    for r in rows:
        r["t_dry_h_r4"] = round_out(r["t_dry_h"])
        r["t_dry_days"] = r["t_dry_h"] / HOURS_PER_DAY
        r["in_S9_window"] = bool(S9_DAYS_MIN <= r["t_dry_days"] <= S9_DAYS_MAX)

    # 择优依据：S9 外部锚点 + Bi 判据 + 逐点交付能力，三条都由实算数值支撑。
    reject = {
        "M1": (f"Bi_m 全程 [{biot['Bi_m_at_D_max']:.4g}, {biot['Bi_m_surface_max']:.4g}]，"
               f"下界已是集中参数适用上限 {BI_LUMPED_MAX} 的 "
               f"{biot['Bi_m_at_D_max'] / BI_LUMPED_MAX:.1f} 倍，内部扩散阻力不可忽略；"
               f"且模型无径向坐标，题面「各处」的逐点判据无法定义，只能退化为 C̄<0.15，"
               f"实算 {rows[0]['t_dry_h']:.4f} h = {rows[0]['t_dry_days']:.2f} d，"
               f"落在题面 S9「{S9_SENTENCE}」之外，也无法交付 21 列的水分场。"),
        "M2": (f"物性冻结的参考状态无法自洽选择：D 在四个同样「合理」的参考状态间相差 "
               f"{biot['D_ref_spread_ratio']:.3g} 倍（全程实测跨 "
               f"{biot['D_span_orders']:.2f} 个数量级），冻结于初始态高估了末期扩散能力，"
               f"实算 {rows[1]['t_dry_h']:.4f} h = {rows[1]['t_dry_days']:.2f} d，"
               f"同样落在 S9 区间之外。"),
        "M4": (f"M4 与 M3 不是竞争关系：附录与域设置由题面按问题指定（问题 3 固定域附录 3、"
               f"问题 4 收缩域附录 4），实算 {rows[3]['t_dry_h']:.4f} h 比 M3 快 "
               f"{abs(rows[3]['rel_vs_M3_pct']):.2f}%，两者各自是对应问题的选定模型。"),
    }
    chosen = {
        "problem_1_2_3": "M3（变物性一维）",
        "problem_4": "M4（变物性 + 动边界）",
        "reasons": [
            f"外部锚点：M3 实算 {rows[2]['t_dry_h']:.4f} h = {rows[2]['t_dry_days']:.2f} d，"
            f"是四个候选中唯一落在题面 S9「{S9_SENTENCE}」区间内的固定域模型。",
            f"Bi 判据：Bi_m 由 {biot['Bi_m_surface_min']:.4g} 升至 "
            f"{biot['Bi_m_surface_max']:.4g}、Bi_h ∈ "
            f"[{biot['Bi_h_min']:.4f}, {biot['Bi_h_max']:.4f}]，两者都在分布参数区间，"
            f"排除集中参数与表面瞬时平衡两种简化。",
            f"交付能力：题面要求每 0.1 cm 的逐点水分场与「各处」判据，只有 1 维及以上模型"
            f"能给出 C(0)-C(R0)={rows[2]['C_center_minus_surface']:.4f} kg/kg 的径向落差。",
            f"守恒精度：M2/M3/M4 的 ε_mass 分别为 {rows[1]['eps_mass']:.3e}、"
            f"{rows[2]['eps_mass']:.3e}、{rows[3]['eps_mass']:.3e}，均在机器精度量级。",
        ],
    }
    return {
        "problem": "cross", "capability": "PA-C2", "n_models": len(rows),
        "reference_model": "M3", "indicators": ["t_dry_h", "radial_capability",
                                               "C_center_minus_surface", "eps_mass",
                                               "can_deliver_pointwise"],
        "criterion_note": ("统一判据 max_r C < 0.15；M1 无径向分布，其判据只能退化为 "
                          "C̄ < 0.15——这正是它被排除的原因之一。"),
        "table": rows, "biot_evidence": biot,
        "rejected_reasons": reject, "chosen": chosen,
        "S9_window_days": [S9_DAYS_MIN, S9_DAYS_MAX],
        "lumped_detail": m1, "frozen_detail": m2,
        "m3_detail": {k: v for k, v in m3.items() if k != "series"},
        "m4_detail": {k: v for k, v in m4.items() if k != "series"},
        "series_by_model": {"M1": m1["series"], "M2": m2["series"],
                           "M3": m3["series"], "M4": m4["series"]},
        "all_values_self_computed": True,
        "provenance": ("M1、M2 由本文件实算（solve_lumped / solve_frozen）；M3、M4 取自"
                      "问题 3、问题 4 的生产算例 JSON。无文献值、无定性占位。"),
    }


def validate_compare(data: dict) -> bool:
    """从写出的表复核：五指标无空缺（M1 的 ε_mass 除外，它没有径向守恒可言）、量级关系正确。"""
    rows = data["table"]
    assert len(rows) == 4 and [r["model"] for r in rows] == ["M1", "M2", "M3", "M4"], \
        "对比表必须是 M1—M4 四行"
    for r in rows:
        assert isinstance(r["t_dry_h"], float) and r["t_dry_h"] > 0.0, \
            f"{r['model']} 的 t_dry 不是实算正数"
        assert isinstance(r["evidence"], str) and r["evidence"].strip(), \
            f"{r['model']} 缺来源说明（PA-C2 要求每个指标可溯源）"
        assert isinstance(r["radial_capability"], bool), f"{r['model']} 径向能力标记缺失"
    m1, m2, m3, m4 = rows
    assert not m1["radial_capability"] and m1["C_center_minus_surface"] == 0.0, \
        "M1 不应有径向分布"
    assert m1["eps_mass"] is None, "M1 是 0 维模型，不应报径向质量守恒误差"
    for r in (m2, m3, m4):
        assert r["eps_mass"] is not None and r["eps_mass"] < 1e-2, \
            f"{r['model']} 质量守恒误差超 1%"
        assert r["C_center_minus_surface"] > 0.0, f"{r['model']} 中心应高于表面"
    assert m1["t_dry_h"] < m2["t_dry_h"] < m3["t_dry_h"], \
        "简化程度越高应越早（错误地）判达标：M1 < M2 < M3"
    assert m3["in_S9_window"] and not m1["in_S9_window"] and not m2["in_S9_window"], \
        "只有 M3 应落在 S9 的 2—3 天区间内"
    assert m4["t_dry_h"] < m3["t_dry_h"], "收缩应加快达标"
    b = data["biot_evidence"]
    assert b["Bi_m_at_D_max"] > b["Bi_m_lumped_limit"], "Bi_m 下界必须超过集中参数适用上限"
    assert b["D_span_orders"] > 1.0, "D 的全程跨度不足一个数量级则冻结物性不算简化"
    assert data["all_values_self_computed"] is True
    return True


def main() -> int:
    ensure_dirs()
    ckpt = CaseCheckpoint(
        ROOT, "2026A/model_compare/M1M2M3M4",
        dependencies=["code/utils.py", "code/params.py", "code/model_compare.py",
                      "figures/problem_3_results.json", "figures/problem_4_results.json",
                      "user_data/附件1.xlsx", "user_data/附件2.xlsx"],
        parameters={"N": N_GRID, "frozen_state": FROZEN_STATE, "out_dt": P3_OUT_DT_S,
                    "t_end": T_MAX_SEARCH_S, "picard_tol": PICARD_TOL,
                    "numpy": np.__version__},
    )
    cached = ckpt.load(validate_compare)
    if cached is None:
        print("[MC] 无可用检查点，开始实算 M1（集中参数）与 M2（冻结物性）")
        m1 = solve_lumped()
        print(f"[MC] M1 t_dry={m1['t_dry_h']:.4f} h  时间常数 τ={m1['tau_h']:.4f} h"
              f"（恒定 C_a 闭式解 {m1['t_dry_h_closed_form_const_Ca']:.4f} h，同量级）")
        m2 = solve_frozen(N_GRID)
        print(f"[MC] M2 t_dry={m2['t_dry_h']:.4f} h  D 冻结值={m2['D_frozen']:.4e}"
              f"  ε_mass={m2['eps_mass']:.3e}")
        data = assemble(m1, m2, m3_indicators(), m4_indicators(), biot_evidence())
        ckpt.save(data, validate_compare)
        print("[MC] 检查点已落盘（校验通过后写入）")
    else:
        data = cached
        print("[MC] 复用已验证检查点")

    print(f"\n{'模型':<6}{'维度':<16}{'t_dry (h)':>11}{'天':>7}{'S9内':>6}"
          f"{'C(0)-C(R)':>11}{'ε_mass':>11}{'逐点交付':>9}")
    for r in data["table"]:
        eps = "—" if r["eps_mass"] is None else f"{r['eps_mass']:.2e}"
        print(f"{r['model']:<6}{r['dim']:<16}{r['t_dry_h']:>11.4f}{r['t_dry_days']:>7.2f}"
              f"{'是' if r['in_S9_window'] else '否':>6}{r['C_center_minus_surface']:>11.4f}"
              f"{eps:>11}{'是' if r['can_deliver_pointwise'] else '否':>9}")
    b = data["biot_evidence"]
    print(f"[MC] Bi_m ∈ [{b['Bi_m_at_D_max']:.4g}, {b['Bi_m_surface_max']:.4g}]"
          f"（集中参数上限 {b['Bi_m_lumped_limit']}）  Bi_h ∈ "
          f"[{b['Bi_h_min']:.4f}, {b['Bi_h_max']:.4f}]")
    print(f"[MC] D 全程跨 {b['D_span_orders']:.2f} 个数量级；四个参考状态间相差 "
          f"{b['D_ref_spread_ratio']:.3g} 倍（M2 冻结态无法自洽选择）")
    print(f"[MC] 择优：问题 1—3 用 {data['chosen']['problem_1_2_3']}，"
          f"问题 4 用 {data['chosen']['problem_4']}")
    save_json(FIG_DIR / "model_compare_results.json", jsonable(data))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
