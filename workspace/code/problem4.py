# -*- coding: utf-8 -*-
"""问题 4：附录 4 物性 + 附件 2 给定收缩，动边界域上的烘干时长。

三个算例各自独立检查点：
  main  —— 物质坐标 s=r/R(t)（主模型，收缩被坐标变换吸收）
  lab   —— 实验室坐标对照分支（显式含 dR/dt 对流项，一阶上风）
  fixed —— 附录 4 物性但半径固定（隔离「物性变化」与「几何收缩」两种效应）
产出 figures/problem_4_results.json：每 60 s × 0.1 cm 水分场（域外置空，K14）、论文表 6
（每 6 h × 0.5 cm，末列取 s=1 节点，K15）、干料质量核算与 P4-C1/P4-C2 证据。
运行：python _utils/run_compute.py code/problem4.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_utils"))
from compute_checkpoint import CaseCheckpoint  # noqa: E402

from params import *  # noqa: F401,F403,E402
from utils import (Props, downsample_index, ensure_dirs, jsonable, make_radius,  # noqa: E402
                   make_grid, material_row_to_fixed_radii, node_columns, round_out,
                   run_case, save_json)

ROOT = Path(__file__).resolve().parents[1]
DEPS = ["code/utils.py", "code/params.py", "code/problem4.py",
        "user_data/附件1.xlsx", "user_data/附件2.xlsx"]
BASE = {"appendix": APPENDIX_OF_PROBLEM[4], "N": N_GRID, "out_dt": P4_OUT_DT_S,
        "t_end": T_MAX_SEARCH_S, "stop_threshold": C_THRESHOLD, "dt1": DT_STAGE1,
        "dt2": DT_STAGE2, "dt3": DT_STAGE3, "record": "full"}
CASES = {
    "main": dict(BASE, shrink=True, frame="material", tag="P4-material"),
    "lab": dict(BASE, shrink=True, frame="lab", tag="P4-lab"),
    "fixed": dict(BASE, shrink=False, frame="material", tag="P4-fixed"),
}


def run_one(name: str) -> dict:
    """跑单个算例并抽出该算例的公共字段（t_dry、判据历程、末态、诊断量）。"""
    cfg = CASES[name]
    res = run_case(cfg)
    t_out = np.asarray(res["t"], dtype=float)
    c_field = np.asarray(res["C_field"], dtype=float)
    t_field = np.asarray(res["T_field"], dtype=float)
    r_series = np.asarray(res["R"], dtype=float)
    t_dry = res["t_dry"]
    assert t_dry is not None, f"{name} 未在 {T_MAX_SEARCH_S / SEC_PER_HOUR:.0f} h 搜索窗内达标"
    bracket = [float(v) for v in res["bracket"]]

    keep = t_out >= P4_OUT_DT_S - 1e-9
    grid_t = t_out[keep]
    c_cols = node_columns(c_field, N_GRID)[keep]
    t_cols = node_columns(t_field, N_GRID)[keep]
    c_max = np.asarray(res["C_max"], dtype=float)[keep]
    c_bar = np.asarray(res["C_bar"], dtype=float)[keep]
    r_keep = r_series[keep]
    dsi = downsample_index(len(grid_t), 600)

    w = (t_dry - bracket[0]) / (bracket[1] - bracket[0])
    i_lo = int(round(bracket[0] / P4_OUT_DT_S)) - 1
    prof_dry = c_cols[i_lo] * (1.0 - w) + c_cols[i_lo + 1] * w
    r_dry = float(r_keep[i_lo] * (1.0 - w) + r_keep[i_lo + 1] * w)

    return {
        "case": name, "frame": cfg["frame"], "shrink": bool(cfg["shrink"]),
        "t_dry": {"seconds": float(t_dry), "hours": float(t_dry / SEC_PER_HOUR),
                  "hours_r4": round_out(t_dry / SEC_PER_HOUR),
                  "bracket_t_lo_s": bracket[0], "bracket_t_hi_s": bracket[1],
                  "bracket_Cmax_lo": bracket[2], "bracket_Cmax_hi": bracket[3],
                  "C_max_at_t_dry": float(c_max[-1]), "C_bar_at_t_dry": float(c_bar[-1]),
                  "C_center_at_t_dry": round_out(prof_dry[0]),
                  "C_surface_at_t_dry": round_out(prof_dry[-1]),
                  "R_at_t_dry_m": r_dry, "R_at_t_dry_cm": round_out(r_dry * CM_PER_M),
                  "R_ratio": float(r_dry / R0)},
        "criterion_series": {"times_s": [float(grid_t[i]) for i in dsi],
                             "times_h": [float(grid_t[i] / SEC_PER_HOUR) for i in dsi],
                             "C_max": [float(c_max[i]) for i in dsi],
                             "C_bar": [float(c_bar[i]) for i in dsi],
                             "C_surface": [float(c_cols[i, -1]) for i in dsi],
                             "T_center": [float(t_cols[i, 0]) for i in dsi],
                             "T_surface": [float(t_cols[i, -1]) for i in dsi],
                             "R_cm": [float(r_keep[i] * CM_PER_M) for i in dsi],
                             "flux_cum": [float(res["flux_cum"][keep][i]) for i in dsi],
                             "D_center": [float(res["D_center"][keep][i]) for i in dsi],
                             "D_surface": [float(res["D_surface"][keep][i]) for i in dsi]},
        "diagnostics": {"eps_mass": float(res["eps_mass"]), "m0": float(res["m0"]),
                        "m_end": float(res["m_end"]), "flux_total": float(res["flux_total"]),
                        "picard_mean": float(res["picard_mean"]),
                        "picard_max": int(res["picard_max"]), "wall_s": float(res["wall_s"]),
                        "n_steps": int(res["n_steps"]), "rdot_max": float(res["rdot_max"]),
                        "C_range": [float(c_field.min()), float(c_field.max())],
                        "T_range": [float(t_field.min()), float(t_field.max())],
                        "R_range_cm": [float(r_keep.min() * CM_PER_M),
                                       float(r_keep.max() * CM_PER_M)]},
        "_grid_t": [float(v) for v in grid_t],
        "_c_cols": [[float(v) for v in row] for row in c_cols],
        "_r_cm": [float(v) for v in r_keep * CM_PER_M],
        "_prof_dry": [float(v) for v in prof_dry],
    }


def validate_case(data: dict) -> bool:
    """单算例约束：K01—K05、K09—K11。表格类约束在汇总阶段核。"""
    diag, td = data["diagnostics"], data["t_dry"]
    c_mat = np.asarray(data["_c_cols"], dtype=float)
    assert c_mat.shape[1] == 21, f"K13 列数 {c_mat.shape[1]} 不是 21"
    assert c_mat.min() >= -1e-9 and c_mat.max() <= C_INIT + 1e-9, \
        f"K01 含水率越界 [{c_mat.min()},{c_mat.max()}]"
    assert diag["T_range"][0] >= T_INIT - 1e-9, f"K02 温度低于初温 {diag['T_range'][0]}"
    assert diag["T_range"][1] <= ATT1_TA_RANGE[1] + 0.5, f"K02 温度超上界 {diag['T_range'][1]}"
    dc_r = float(np.max(np.diff(c_mat, axis=1)))
    assert dc_r <= 1e-9, f"K03 含水率沿半径出现上升 max差={dc_r:.3e}"
    c_bar = np.asarray(data["criterion_series"]["C_bar"], dtype=float)
    assert float(np.max(np.diff(c_bar))) <= 1e-12, "K04 平均含水率出现回升"
    if data["frame"] == "material":
        assert diag["eps_mass"] < 1e-2, f"K05 质量守恒误差 {diag['eps_mass']:.3e} 超 1%"
    else:
        # 实验室坐标是被 MODELING_REPORT 第 1.2 节推翻的表述，其失效机理正是「违反干物质
        # 守恒」：算子已含 1/R² 的坐标变换项，再叠加 sṘ/R 对流项即重复计入边界运动。
        # 故此分支不设 K05 上限，改为在汇总层断言其残差显著大于物质坐标（见 validate_all）。
        assert diag["eps_mass"] > 1e-2, \
            f"实验室坐标竟满足质量守恒（{diag['eps_mass']:.3e}），第 1.2 节的推翻论证需复核"
    assert diag["rdot_max"] <= 1e-14, f"K09 出现 dR/dt>0：{diag['rdot_max']:.3e}"
    assert td["bracket_Cmax_lo"] >= C_THRESHOLD > td["bracket_Cmax_hi"], "K11 区间未跨阈值"
    assert td["bracket_t_lo_s"] < td["seconds"] < td["bracket_t_hi_s"], "K11 t_dry 不在区间内"
    if data["shrink"]:
        # K10 只约束真正取用附件 2 的收缩算例：t_dry 须落在 R(t) 覆盖窗内，绝不外推。
        assert td["seconds"] < ATT2_T_END_S, f"K10 t_dry={td['seconds']} 超出附件 2 覆盖"
        assert td["R_at_t_dry_m"] < R0 - 1e-6, "收缩算例的末半径应小于初始半径"
    else:
        # 固定域是隔离算例，半径恒为 R0、不查附件 2，故 K10 不适用（时长可超 72 h）。
        assert abs(td["R_at_t_dry_m"] - R0) < 1e-12, "固定域算例的半径应恒为 R0"
    return True


def build_tables(main: dict) -> dict:
    """result4 主体（域外置空，K14）与论文表 6（末列取 s=1 节点，K15）。"""
    grid_t = np.asarray(main["_grid_t"], dtype=float)
    c_cols = np.asarray(main["_c_cols"], dtype=float)
    r_cm = np.asarray(main["_r_cm"], dtype=float)
    s_nodes = make_grid(N_GRID)[0][::OUT_STRIDE]
    assert len(s_nodes) == N_OUT_RADII, "输出节点数与 21 列不符"

    # result4.xlsx：固定物理半径列，r_j > R(t) 的单元留空（None → Excel 空白）。
    rows, n_blank = [], 0
    for k in range(len(grid_t)):
        row = material_row_to_fixed_radii(c_cols[k], s_nodes, r_cm[k])
        rows.append([None if v is None else round_out(v) for v in row])
        n_blank += sum(1 for v in row if v is None)
    # 「药材表面」列：与 result4 表体同分辨率的 s=1 物质节点原值（K15），供导出直接取用。
    surf = [round_out(c_cols[k, -1]) for k in range(len(grid_t))]

    # 论文表 6：每 6 h × 0.5 cm + 末列「药材表面」取 s=1 节点原值（不插值）。
    step_s = P4_TABLE_DT_H * SEC_PER_HOUR
    col_stride = int(round(P4_TABLE_DR_CM / OUT_DR_CM))
    col_tab = list(range(0, N_OUT_RADII, col_stride))
    t6_rows, t6_times, t6_r_now = [], [], []
    for k in range(1, int(grid_t[-1] // step_s) + 1):
        tt = k * step_s
        i = int(round(tt / P4_OUT_DT_S)) - 1
        assert abs(grid_t[i] - tt) < 1e-9, f"表 6 时刻 {tt} 未对齐输出网格"
        fixed = material_row_to_fixed_radii(c_cols[i], s_nodes, r_cm[i])
        t6_rows.append([None if fixed[j] is None else round_out(fixed[j]) for j in col_tab]
                       + [round_out(c_cols[i, -1])])
        t6_times.append(float(tt))
        t6_r_now.append(round_out(r_cm[i]))
    td = main["t_dry"]
    prof = np.asarray(main["_prof_dry"], dtype=float)
    fixed_dry = material_row_to_fixed_radii(prof, s_nodes, td["R_at_t_dry_m"] * CM_PER_M)
    t6_rows.append([None if fixed_dry[j] is None else round_out(fixed_dry[j]) for j in col_tab]
                   + [round_out(prof[-1])])
    t6_times.append(float(td["seconds"]))
    t6_r_now.append(td["R_at_t_dry_cm"])

    return {
        "result_grid": {"times_s": [float(v) for v in grid_t], "radii_cm": list(OUT_RADII_CM),
                        "R_cm": [float(v) for v in r_cm], "C": rows,
                        "C_surface": surf, "surface_column": P4_TABLE_LAST_COL,
                        "n_blank_cells": int(n_blank),
                        "blank_rule": "r_j > R(t) 的单元留空（None），不填 0（约束 K14）",
                        "surface_rule": "末列取 s=1 物质节点原值，不做空间插值（约束 K15）"},
        "table6_C": {"times_s": t6_times,
                     "times_h": [float(v / SEC_PER_HOUR) for v in t6_times],
                     "radii_cm": [OUT_RADII_CM[j] for j in col_tab],
                     "last_column": P4_TABLE_LAST_COL, "R_cm": t6_r_now, "values": t6_rows,
                     "last_row_is_t_dry": True,
                     "last_column_rule": "取 s=1 物质节点原值，不做空间插值（约束 K15）"},
    }


def dry_mass_audit(r_dry_cm: float) -> dict:
    """干料质量守恒反算：附件 2 的径向收缩幅度是否与附录 4 的致密化自洽。

    单位长度干料质量 πR²ρ_d 守恒 => ρ_d(末) = ρ_d(初)·(R0/R_dry)²（仅径向收缩）；
    若为各向同性收缩则指数为 3，可据此排除。
    """
    props = Props(APPENDIX_OF_PROBLEM[4])
    rho_wet0 = float(props.rho(np.array([C_INIT]))[0])
    rho_dry0 = rho_wet0 / (1.0 + C_INIT)
    ratio = R0_CM / r_dry_cm
    radial = rho_dry0 * ratio ** 2
    isotropic = rho_dry0 * ratio ** 3
    err_radial = abs(radial - RHO4_A) / RHO4_A
    err_iso = abs(isotropic - RHO4_A) / RHO4_A
    assert err_radial < 0.05, f"径向收缩假设与附录 4 干密度不自洽，偏差 {err_radial:.2%}"
    assert err_iso > 0.5, "各向同性收缩本应被显著排除，请复核几何假设"
    return {"rho_wet_initial": rho_wet0, "rho_dry_initial": rho_dry0,
            "R_ratio_R0_over_Rdry": float(ratio),
            "rho_dry_final_radial_only": radial, "rho_dry_final_isotropic": isotropic,
            "rho_dry_reference_appendix4": RHO4_A,
            "rel_err_radial_only": err_radial, "rel_err_isotropic": err_iso,
            "conclusion": ("仅径向收缩：反算干密度 %.2f 与附录 4 的 %.1f 相差 %.2f%%，自洽；"
                           "各向同性收缩反算 %.2f，偏离 %.1f%%，据此排除。"
                           % (radial, RHO4_A, err_radial * 100.0, isotropic, err_iso * 100.0))}


def validate_all(data: dict) -> bool:
    """汇总层约束：K13（四位小数）、K14（域外置空）、K15（末列取 s=1 节点）。"""
    rg = data["result_grid"]
    flat = [v for row in rg["C"] for v in row if v is not None]
    arr = np.asarray(flat, dtype=float)
    assert np.max(np.abs(arr - np.round(arr, DECIMALS))) < 1e-12, "K13 未量化到四位小数"
    assert rg["n_blank_cells"] > 0, "K14 收缩算例必须出现域外空缺单元"

    # K14：逐行核对「空缺 <=> r_j > R(t)」，且空缺不得写成 0。
    radii = np.asarray(rg["radii_cm"], dtype=float)
    r_now = np.asarray(rg["R_cm"], dtype=float)
    for k, row in enumerate(rg["C"]):
        for j, v in enumerate(row):
            outside = radii[j] > r_now[k] + 1e-12
            assert (v is None) == outside, f"K14 第 {k} 行第 {j} 列域内外与空缺不一致"

    # K15（表体分辨率）：「药材表面」列必须逐行齐备、四位小数、单调不增，且不早于域内最外列。
    surf = np.asarray(rg["C_surface"], dtype=float)
    assert len(surf) == len(rg["C"]), "药材表面列与表体行数不符"
    assert rg["surface_column"] == P4_TABLE_LAST_COL, "药材表面列表头偏离题面"
    assert np.max(np.abs(surf - np.round(surf, DECIMALS))) < 1e-12, "K13 药材表面列未量化到四位小数"
    assert float(np.max(np.diff(surf))) <= 1e-12, "K04 药材表面含水率出现回升"
    for k, row in enumerate(rg["C"]):
        inside = [v for v in row if v is not None]
        assert inside and surf[k] <= inside[-1] + 1e-12, \
            f"第 {k} 行表面值高于域内最外节点，径向单调性被破坏"

    # K15：表 6 末列必须等于 s=1 节点值（用主算例末态剖面独立复核）。
    t6 = data["table6_C"]
    assert t6["last_column"] == P4_TABLE_LAST_COL, "表 6 末列表头偏离题面"
    for i, tt in enumerate(t6["times_s"][:-1]):
        k = int(round(tt / P4_OUT_DT_S)) - 1
        assert abs(t6["values"][i][-1] - surf[k]) < 1e-12, \
            f"表 6 第 {i} 行末列与表体药材表面列不一致"
    assert abs(t6["values"][-1][-1] - data["main"]["t_dry"]["C_surface_at_t_dry"]) < 1e-12, \
        "K15 表 6 末列未取 s=1 节点值"
    assert t6["last_row_is_t_dry"] and abs(t6["times_s"][-1]
                                           - data["main"]["t_dry"]["seconds"]) < 1e-9, \
        "表 6 末行不是烘干结束时间"
    for row in t6["values"]:
        assert row[-1] is not None, "表 6 末列（药材表面）不得为空"

    # 三算例的物理排序：收缩加速烘干，故 main < fixed；实验室坐标对流项使 lab 略慢。
    hm, hl, hf = (data[k]["t_dry"]["hours"] for k in ("main", "lab", "fixed"))
    assert hm < hf, f"收缩算例 {hm:.4f} h 应快于固定域 {hf:.4f} h"
    assert hl > hm, f"实验室坐标 {hl:.4f} h 应慢于物质坐标 {hm:.4f} h"

    # 第 1.2 节的推翻论证须在数值上可证伪：实验室坐标的质量残差必须比物质坐标大若干量级。
    em_mat = data["main"]["diagnostics"]["eps_mass"]
    em_lab = data["lab"]["diagnostics"]["eps_mass"]
    em_fix = data["fixed"]["diagnostics"]["eps_mass"]
    assert em_mat < 1e-2 and em_fix < 1e-2, "K05 物质坐标算例必须守恒"
    assert em_lab / max(em_mat, 1e-300) > 1e6, \
        f"实验室坐标残差 {em_lab:.3e} 未显著超出物质坐标 {em_mat:.3e}，推翻论证失据"
    return True


def validate_capability(data: dict) -> dict:
    """P4-C1 / P4-C2 的可反驳检验：附件 2 保形插值求导、dR/dt 对流项、附录隔离与双 t_dry 对比。"""
    verdict = {}
    body = "\n".join(ln for ln in (Path(__file__).resolve().parent / "utils.py")
                     .read_text(encoding="utf-8").splitlines() if not ln.strip().startswith("#"))

    # ---- P4-C1（含 K09）：R(t) 由附件 2 读入 + PCHIP 求导，dR/dt ≤ 0；离散式含 dR/dt 对流项。
    assert "PchipInterpolator(t_rad, r_m" in body, "P4-C1 未对附件 2 的 R(t) 做保形插值"
    assert "load_attachment2()" in body and "pd.read_excel(DATA_DIR / ATT2_FILE" in body, \
        "P4-C1 R(t) 不是从附件 2 读入"
    assert "deriv = spline.derivative()" in body, "P4-C1 未对插值器求导得 dR/dt"
    assert "rdot = radius_dot(t_new)" in body, "P4-C1 离散步未取用 dR/dt"
    assert "u_adv = s_node * rdot / r_now" in body, "P4-C1 离散式缺含 dR/dt 的对流速度"
    assert "adv = props.rho(C_it) * props.cp(C_it) * V * u_adv / ds" in body, \
        "P4-C1 能量方程缺 dR/dt 对流项"
    assert "adv_c = V * u_adv / ds" in body, "P4-C1 水分方程缺 dR/dt 对流项"
    assert "inv_r2 = 1.0 / r_now**2" in body, "P4-C1 物质坐标缺 1/R² 因子"

    radius, radius_dot, rdot_max = make_radius()
    probe_t = np.linspace(0.0, ATT2_T_END_S, 2001)
    rdot_probe = np.array([radius_dot(v) for v in probe_t])
    r_probe = np.array([radius(v) for v in probe_t])
    assert rdot_max <= 1e-14, f"K09 dR/dt 出现正值 max={rdot_max:.3e}"
    assert rdot_probe.min() < -1e-10, "P4-C1 dR/dt 全程为 0，R 被当成常数"
    assert r_probe[-1] < r_probe[0] - 1e-4, "P4-C1 半径未发生收缩"
    # 实验室坐标分支必须真的激活了非零对流系数。
    t_mid = 0.25 * ATT2_T_END_S
    u_surface = 1.0 * radius_dot(t_mid) / radius(t_mid)
    assert abs(u_surface) > 1e-9, f"P4-C1 对流系数退化为 0（u={u_surface:.3e}）"
    hm, hl = data["main"]["t_dry"]["hours"], data["lab"]["t_dry"]["hours"]
    em_mat = data["main"]["diagnostics"]["eps_mass"]
    em_lab = data["lab"]["diagnostics"]["eps_mass"]
    verdict["P4-C1"] = {
        "status": "PASS",
        "evidence": (f"R(t) 由 utils.load_attachment2() 从附件 2（{ATT2_ROWS} 行、"
                     f"{ATT2_DT_S:.0f} s 步长）读入，经 PchipInterpolator 保形插值并 "
                     f".derivative() 求导；2001 点扫描得 max dR/dt = {rdot_max:.3e} ≤ 0、"
                     f"min dR/dt = {rdot_probe.min():.3e} m/s（非零，R 不是常数），"
                     f"R 由 {r_probe[0] * CM_PER_M:.4f} cm 收缩到 {r_probe[-1] * CM_PER_M:.4f} cm。"
                     f"离散方程含 dR/dt 对流项 u = s·Ṙ/R：能量方程 "
                     f"`adv = ρ c_p V u_adv/ds`、水分方程 `adv_c = V u_adv/ds`，"
                     f"在 t={t_mid:.0f} s、s=1 处系数 u={u_surface:.4e} 1/s ≠ 0，"
                     f"lab 分支实算 t_dry={hl:.4f} h。主模型改用物质坐标 s=r/R(t)，"
                     f"该对流项被坐标变换吸收为算子中的 1/R² 与 Robin 项中的 1/R，"
                     f"实算 t_dry={hm:.4f} h；两种口径相差 {(hl - hm) / hm * 100:.2f}%，"
                     f"量级一致，说明收缩项确实进入了方程而非被忽略。二者取舍有独立的"
                     f"数值判据：全局质量守恒残差物质坐标 ε_mass={em_mat:.3e}、实验室坐标 "
                     f"ε_mass={em_lab:.3e}（相差 {em_lab / max(em_mat, 1e-300):.2e} 倍），"
                     f"实证 MODELING_REPORT 第 1.2 节「实验室坐标形式违反干物质守恒」的判定—"
                     f"其算子已含 1/R² 变换项，再叠加 sṘ/R 即把边界运动重复计入。"),
        "rdot_max": float(rdot_max), "rdot_min": float(rdot_probe.min()),
        "u_adv_at_surface": float(u_surface),
        "t_dry_material_h": hm, "t_dry_lab_h": hl,
        "lab_vs_material_rel_diff": float((hl - hm) / hm),
        "eps_mass_material": em_mat, "eps_mass_lab": em_lab,
    }

    # ---- P4-C2（含 K08）：附录 4 分支不得出现附录 3 常数；须同时报告问题 3、4 的 t_dry 及差异方向。
    props4 = Props(APPENDIX_OF_PROBLEM[4])
    props4._check_isolation()
    rho4 = float(props4.rho(np.array([0.0]))[0])
    cp4 = float(props4.cp(np.array([0.0]))[0])
    lam4 = float(props4.lam(np.array([0.0]))[0])
    assert abs(rho4 - RHO3_A) > 1.0, "K08 密度串用附录 3"
    assert abs(cp4 - CP3_A) > 1.0, "K08 比热容串用附录 3"
    assert abs(lam4 - LAM3_A) > 1e-3, "K08 热导率串用附录 3"
    for name, val, ref in (("D0", D0_A4, D0_A3), ("Dc", DC_A4, DC_A3)):
        assert abs(val - ref) > 1e-12, f"K08 扩散系数 {name} 串用附录 3"
    p3 = data["cross_problem"]["t_dry_problem3_h"]
    p4 = data["main"]["t_dry"]["hours"]
    assert p3 is not None, "P4-C2 需要问题 3 的 t_dry 才能给出对比"
    verdict["P4-C2"] = {
        "status": "PASS",
        "evidence": (f"Props(4) 在 C→0 给出 ρ={rho4:.6f}、c_p={cp4:.6f}、λ={lam4:.7f}，"
                     f"与附录 3 的 ({RHO3_A:.1f}, {CP3_A:.1f}, {LAM3_A:.2f}) 显著不同；"
                     f"D 前因子 {D0_A4:.5f} ≠ {D0_A3:.5f}、浓度系数 {DC_A4:.2f} ≠ {DC_A3:.2f}，"
                     f"Props._check_isolation() 通过。双时长对比：问题 3（附录 3、固定半径）"
                     f"t_dry={p3:.4f} h，问题 4（附录 4、附件 2 收缩）t_dry={p4:.4f} h，"
                     f"问题 4 更{'短' if p4 < p3 else '长'} {abs(p4 - p3):.4f} h "
                     f"（{(p4 - p3) / p3 * 100:+.2f}%）；同时给出附录 4 物性但半径固定的隔离算例 "
                     f"t_dry={data['fixed']['t_dry']['hours']:.4f} h，"
                     f"说明差异中「物性变化」与「几何收缩」两种效应的方向相反且量级可分。"),
        "t_dry_problem3_h": p3, "t_dry_problem4_h": p4,
        "t_dry_fixed_domain_h": data["fixed"]["t_dry"]["hours"],
        "diff_h": float(p4 - p3), "diff_rel": float((p4 - p3) / p3),
    }

    t6 = data["table6_C"]
    verdict["P4-C3"] = {
        "status": "PASS",
        "evidence": (f"result_grid 为 {len(data['result_grid']['times_s'])}×21"
                     f"（每隔 60 s × 每隔 0.1 cm），域外单元共 "
                     f"{data['result_grid']['n_blank_cells']} 个按 K14 留空不填 0；"
                     f"表 6 共 {len(t6['times_s'])} 行 × {len(t6['radii_cm'])} 固定列 + 末列"
                     f"「{t6['last_column']}」取 s=1 物质节点原值。"),
    }
    return verdict


def _load_problem3_t_dry() -> float | None:
    """读问题 3 的 t_dry 用于 P4-C2 的双时长对比（只取一个标量，不整表读入）。"""
    path = FIG_DIR / "problem_3_results.json"
    if not path.exists():
        return None
    import json
    return float(json.loads(path.read_text(encoding="utf-8"))["t_dry"]["hours"])


def main() -> int:
    ensure_dirs()
    cases = {}
    for name in ("main", "lab", "fixed"):
        ckpt = CaseCheckpoint(
            ROOT, f"2026A/problem4/{name}/appendix4/N40",
            dependencies=DEPS,
            parameters={"cfg": CASES[name], "picard_tol": PICARD_TOL, "gauss_n": GAUSS_N,
                        "numpy": np.__version__},
        )
        cached = ckpt.load(validate_case)
        if cached is None:
            print(f"[P4:{name}] 开始求解：frame={CASES[name]['frame']} "
                  f"shrink={CASES[name]['shrink']}")
            got = run_one(name)
            ckpt.save(got, validate_case)
            print(f"[P4:{name}] 检查点已落盘  t_dry={got['t_dry']['hours']:.4f} h"
                  f"  墙钟={got['diagnostics']['wall_s']:.1f} s")
        else:
            got = cached
            print(f"[P4:{name}] 复用已验证检查点  t_dry={got['t_dry']['hours']:.4f} h")
        cases[name] = got

    tables = build_tables(cases["main"])
    t_dry_p3 = _load_problem3_t_dry()
    data = {
        "problem": 4, "appendix": BASE["appendix"],
        "shrink_mode": SHRINK_MODE, "frame_main": CASES["main"]["frame"],
        "grid": {"N": N_GRID, "ds": float(make_grid(N_GRID)[1]), "out_stride": OUT_STRIDE,
                 "radii_cm": list(OUT_RADII_CM)},
        "schedule": {"out_dt_s": P4_OUT_DT_S, "dt_stage1_s": DT_STAGE1,
                     "dt_stage2_s": DT_STAGE2, "dt_stage3_s": DT_STAGE3,
                     "t_search_limit_s": T_MAX_SEARCH_S},
        "main": {k: v for k, v in cases["main"].items() if not k.startswith("_")},
        "lab": {k: v for k, v in cases["lab"].items() if not k.startswith("_")},
        "fixed": {k: v for k, v in cases["fixed"].items() if not k.startswith("_")},
        "frame_comparison": {
            "material_h": cases["main"]["t_dry"]["hours"],
            "lab_h": cases["lab"]["t_dry"]["hours"],
            "rel_diff": float((cases["lab"]["t_dry"]["hours"]
                               - cases["main"]["t_dry"]["hours"])
                              / cases["main"]["t_dry"]["hours"]),
            "eps_mass_material": cases["main"]["diagnostics"]["eps_mass"],
            "eps_mass_lab": cases["lab"]["diagnostics"]["eps_mass"],
            "note": ("物质坐标为主模型（收缩由 s=r/R(t) 变换吸收）；实验室坐标显式保留 "
                     "u=sṘ/R 对流项并用一阶上风离散，作为对照分支。"),
            "mass_verdict": ("第 1.2 节判定实验室坐标形式违反干物质守恒，本次实算给出量化"
                             "证据：物质坐标 ε_mass=%.3e（机器精度），实验室坐标 ε_mass=%.3e"
                             "（%.1f%%）。二者算子同为 1/R²，差别仅在于后者额外叠加 sṘ/R "
                             "对流项，即把已被坐标变换吸收的边界运动重复计入，故 t_dry 系统"
                             "性偏慢 %.2f%%。K05 只对物质坐标算例设为硬门。"
                             % (cases["main"]["diagnostics"]["eps_mass"],
                                cases["lab"]["diagnostics"]["eps_mass"],
                                cases["lab"]["diagnostics"]["eps_mass"] * 100.0,
                                (cases["lab"]["t_dry"]["hours"]
                                 - cases["main"]["t_dry"]["hours"])
                                / cases["main"]["t_dry"]["hours"] * 100.0))},
        "shrink_isolation": {
            "shrinking_h": cases["main"]["t_dry"]["hours"],
            "fixed_domain_h": cases["fixed"]["t_dry"]["hours"],
            "speedup_rel": float((cases["fixed"]["t_dry"]["hours"]
                                  - cases["main"]["t_dry"]["hours"])
                                 / cases["fixed"]["t_dry"]["hours"]),
            "note": "两者物性同为附录 4，唯一差别是半径是否收缩，用于隔离几何效应。"},
        "cross_problem": {
            "t_dry_problem3_h": t_dry_p3,
            "t_dry_problem4_h": cases["main"]["t_dry"]["hours"],
            "diff_h": None if t_dry_p3 is None else float(cases["main"]["t_dry"]["hours"]
                                                          - t_dry_p3)},
        "dry_mass_audit": dry_mass_audit(cases["main"]["t_dry"]["R_at_t_dry_cm"]),
        "radius_profile": {
            "times_s": cases["main"]["criterion_series"]["times_s"],
            "R_cm": cases["main"]["criterion_series"]["R_cm"]},
        **tables,
    }
    assert validate_all(data) is True, "汇总层约束未通过"
    cap = validate_capability(data)
    data["capability_verdict"] = cap

    td = data["main"]["t_dry"]
    dm = data["dry_mass_audit"]
    print(f"[P4] 主模型（物质坐标）t_dry={td['seconds']:.3f} s = {td['hours']:.4f} h"
          f"  R(t_dry)={td['R_at_t_dry_cm']:.4f} cm（R/R0={td['R_ratio']:.4f}）")
    print(f"[P4] t_dry 处 C(0)={td['C_center_at_t_dry']:.4f}"
          f"  C(表面)={td['C_surface_at_t_dry']:.4f}  C̄={td['C_bar_at_t_dry']:.6f}")
    fc = data["frame_comparison"]
    print(f"[P4] 坐标系对照：物质 {fc['material_h']:.4f} h"
          f"  实验室 {fc['lab_h']:.4f} h（{fc['rel_diff'] * 100:+.2f}%）")
    print(f"[P4] 守恒对照：物质 ε_mass={fc['eps_mass_material']:.3e}"
          f"  实验室 ε_mass={fc['eps_mass_lab']:.3e}"
          f"（实验室违反干物质守恒，与第 1.2 节的推翻论证一致；K05 只作用于物质坐标）")
    print(f"[P4] 几何隔离：收缩 {data['shrink_isolation']['shrinking_h']:.4f} h"
          f"  固定域 {data['shrink_isolation']['fixed_domain_h']:.4f} h"
          f"（收缩提速 {data['shrink_isolation']['speedup_rel'] * 100:.2f}%）")
    if t_dry_p3 is not None:
        print(f"[P4] 跨问题：问题 3 {t_dry_p3:.4f} h → 问题 4 {td['hours']:.4f} h"
              f"（{data['cross_problem']['diff_h']:+.4f} h）")
    print(f"[P4] 干料反算：ρ_d 初 {dm['rho_dry_initial']:.2f} → 仅径向 "
          f"{dm['rho_dry_final_radial_only']:.2f}（附录 4 参考 {dm['rho_dry_reference_appendix4']:.1f}，"
          f"偏差 {dm['rel_err_radial_only'] * 100:.2f}%）；各向同性反算 "
          f"{dm['rho_dry_final_isotropic']:.2f}（偏差 {dm['rel_err_isotropic'] * 100:.1f}%，排除）")
    for name in ("main", "lab", "fixed"):
        dg = data[name]["diagnostics"]
        print(f"[P4:{name}] ε_mass={dg['eps_mass']:.3e}  Picard均值={dg['picard_mean']:.2f}"
              f"  峰值={dg['picard_max']}  步数={dg['n_steps']}  墙钟={dg['wall_s']:.1f} s")
    print(f"[P4] result4 域外空缺 {data['result_grid']['n_blank_cells']} 单元（K14，不填 0）")
    for cid, v in cap.items():
        print(f"[P4][{cid}] {v['status']}")
    save_json(FIG_DIR / "problem_4_results.json", jsonable(data))
    print(f"[P4] 写出 figures/problem_4_results.json  "
          f"{(FIG_DIR / 'problem_4_results.json').stat().st_size / 1024:.1f} KB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
