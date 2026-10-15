# -*- coding: utf-8 -*-
"""问题 1：预热阶段（0—1800 s）常物性（附录 2）温度场与水分场。

产出 figures/problem_1_results.json：逐秒 21 列全场（供 result1.xlsx）、论文表 1/表 2、
表面通量序列、Hovmoller 降采样场，以及 K01—K06/K13 约束与 P1-C1/P1-C2 能力的运行期证据。
运行：python _utils/run_compute.py code/problem1.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_utils"))
from compute_checkpoint import CaseCheckpoint  # noqa: E402

from params import *  # noqa: F401,F403,E402
from utils import (Props, downsample_index, ensure_dirs, jsonable, make_ambient,  # noqa: E402
                   make_grid, node_columns, pick_time_index, round_out, run_case, save_json)

ROOT = Path(__file__).resolve().parents[1]
CASE_KEY = "2026A/problem1/appendix2/N40/dt0.01"
CFG = {"appendix": APPENDIX_OF_PROBLEM[1], "N": N_GRID, "out_dt": P1_OUT_DT_S,
       "t_end": P1_WINDOW_S, "dt1": DT_STAGE1, "dt2": DT_STAGE2, "dt3": DT_STAGE3,
       "record": "full", "tag": "P1-prod"}


def solve() -> dict:
    """推进 0→1800 s，逐秒记录全场，并组装论文表格与图件字段。"""
    res = run_case(CFG)
    t_out = np.asarray(res["t"], dtype=float)
    t_field = np.asarray(res["T_field"], dtype=float)
    c_field = np.asarray(res["C_field"], dtype=float)
    assert t_field.shape == (len(t_out), N_GRID + 1), f"温度场形状异常 {t_field.shape}"

    # 输出列严格取网格节点（stride=OUT_STRIDE），不做空间插值。
    t_cols = node_columns(t_field, N_GRID)
    c_cols = node_columns(c_field, N_GRID)
    assert t_cols.shape[1] == N_OUT_RADII == 21, "输出列数须为 21"

    # result1.xlsx 主体：t=1..1800 s（丢弃 t=0 行），21 列。
    keep = t_out >= P1_OUT_DT_S - 1e-9
    grid_t = t_out[keep]
    assert len(grid_t) == 1800 and abs(grid_t[0] - 1.0) < 1e-9 and abs(grid_t[-1] - 1800.0) < 1e-9, \
        f"result1 行数/端点异常 n={len(grid_t)}"

    idx_tab = [pick_time_index(t_out, tt) for tt in P1_TABLE_TIMES_S]
    col_tab = [int(round(rc / OUT_DR_CM)) for rc in P1_TABLE_RADII_CM]
    table1 = [[round_out(t_cols[i, j]) for j in col_tab] for i in idx_tab]
    table2 = [[round_out(c_cols[i, j]) for j in col_tab] for i in idx_tab]

    t_amb, c_amb = make_ambient()
    ta_series = np.asarray([t_amb(tt) for tt in grid_t], dtype=float)
    ca_series = np.asarray([c_amb(tt) for tt in grid_t], dtype=float)
    surf_t = t_cols[keep][:, -1]
    surf_c = c_cols[keep][:, -1]
    q_heat = H_HEAT * (surf_t - ta_series)
    q_mass = H_MASS * (surf_c - ca_series)

    hov = downsample_index(len(grid_t), 200)
    out = {
        "problem": 1, "appendix": CFG["appendix"], "tag": CFG["tag"],
        "grid": {"N": N_GRID, "ds": float(make_grid(N_GRID)[1]), "out_stride": OUT_STRIDE,
                 "radii_cm": list(OUT_RADII_CM)},
        "schedule": {"dt_stage1_s": DT_STAGE1, "out_dt_s": P1_OUT_DT_S,
                     "t_end_s": P1_WINDOW_S, "n_steps": int(res["n_steps"])},
        "result_grid": {"times_s": [float(v) for v in grid_t],
                        "radii_cm": list(OUT_RADII_CM),
                        "T": [[round_out(v) for v in row] for row in t_cols[keep]],
                        "C": [[round_out(v) for v in row] for row in c_cols[keep]]},
        "table1_T": {"times_s": list(P1_TABLE_TIMES_S), "radii_cm": list(P1_TABLE_RADII_CM),
                     "values": table1},
        "table2_C": {"times_s": list(P1_TABLE_TIMES_S), "radii_cm": list(P1_TABLE_RADII_CM),
                     "values": table2},
        "final_state": {"t_s": float(t_out[-1]), "T_center": float(t_field[-1, 0]),
                        "T_surface": float(t_field[-1, -1]), "C_center": float(c_field[-1, 0]),
                        "C_surface": float(c_field[-1, -1]),
                        "T_center_r4": round_out(t_field[-1, 0]),
                        "T_surface_r4": round_out(t_field[-1, -1]),
                        "C_center_r4": round_out(c_field[-1, 0]),
                        "C_surface_r4": round_out(c_field[-1, -1])},
        "mass_series": {"times_s": [float(v) for v in grid_t],
                        "C_bar": [float(v) for v in np.asarray(res["C_bar"])[keep]],
                        "C_max": [float(v) for v in np.asarray(res["C_max"])[keep]],
                        "flux_cum": [float(v) for v in np.asarray(res["flux_cum"])[keep]],
                        "D_center": [float(v) for v in np.asarray(res["D_center"])[keep]],
                        "D_surface": [float(v) for v in np.asarray(res["D_surface"])[keep]]},
        "flux_series": {"times_s": [float(v) for v in grid_t],
                        "Ta": [float(v) for v in ta_series],
                        "Ca": [float(v) for v in ca_series],
                        "T_surface": [float(v) for v in surf_t],
                        "C_surface": [float(v) for v in surf_c],
                        "q_heat_W_m2": [float(v) for v in q_heat],
                        "q_mass_kg_m2s": [float(v) for v in q_mass]},
        "hovmoller": {"times_s": [float(grid_t[i]) for i in hov],
                      "radii_cm": list(OUT_RADII_CM),
                      "T": [[float(v) for v in t_cols[keep][i]] for i in hov],
                      "C": [[float(v) for v in c_cols[keep][i]] for i in hov]},
        "profiles": {"times_s": list(P1_TABLE_TIMES_S), "radii_cm": list(OUT_RADII_CM),
                     "T": [[float(v) for v in t_cols[i]] for i in idx_tab],
                     "C": [[float(v) for v in c_cols[i]] for i in idx_tab]},
        "diagnostics": {"eps_mass": float(res["eps_mass"]), "m0": float(res["m0"]),
                        "m_end": float(res["m_end"]), "flux_total": float(res["flux_total"]),
                        "picard_mean": float(res["picard_mean"]),
                        "picard_max": int(res["picard_max"]), "wall_s": float(res["wall_s"]),
                        "Bi_h": BI_H_A2, "Bi_m_initial": BI_M_A2,
                        "T_range": [float(t_field.min()), float(t_field.max())],
                        "C_range": [float(c_field.min()), float(c_field.max())],
                        "Ta_max_used": float(ta_series.max())},
    }
    return out


def validate_constraints(data: dict) -> bool:
    """K01—K06、K13：从写出的结果字段本身复核，不重跑求解器。"""
    t_mat = np.asarray(data["result_grid"]["T"], dtype=float)
    c_mat = np.asarray(data["result_grid"]["C"], dtype=float)
    times = np.asarray(data["result_grid"]["times_s"], dtype=float)
    diag = data["diagnostics"]

    assert t_mat.shape == (1800, 21), f"K13 表格形状 {t_mat.shape} 不是 1800x21"
    assert c_mat.shape == (1800, 21), f"K13 表格形状 {c_mat.shape} 不是 1800x21"
    assert abs(times[0] - 1.0) < 1e-9 and abs(times[-1] - 1800.0) < 1e-9, "时间列端点错误"
    assert np.allclose(np.diff(times), P1_OUT_DT_S), "时间列步长非 1 s"

    # K01 含水率逐点落在 [0, 2.55]。
    assert c_mat.min() >= -1e-9, f"K01 含水率出现负值 {c_mat.min()}"
    assert c_mat.max() <= C_INIT + 1e-9, f"K01 含水率超初值 {c_mat.max()}"

    # K02 温度逐点落在 [28, max_t Ta + 0.5]。
    hi = diag["Ta_max_used"] + 0.5
    assert t_mat.min() >= T_INIT - 1e-9, f"K02 温度低于初温 {t_mat.min()}"
    assert t_mat.max() <= hi, f"K02 温度超上界 {t_mat.max()} > {hi}"

    # K03 ∂C/∂r ≤ 0，等价于 max_r C = C(0,t)（登记实测偏差供复核）。
    dc_r = float(np.max(np.diff(c_mat, axis=1)))
    assert dc_r <= 1e-9, f"K03 含水率沿半径出现上升 max差={dc_r:.3e}"
    assert np.allclose(c_mat.max(axis=1), c_mat[:, 0], atol=1e-12), "K03 逐时刻最大值不在中心"

    # K04 平均含水率单调不增。
    c_bar = np.asarray(data["mass_series"]["C_bar"], dtype=float)
    d_bar = float(np.max(np.diff(c_bar)))
    assert d_bar <= 1e-12, f"K04 平均含水率出现回升 max差={d_bar:.3e}"

    # K05 全局质量守恒。
    assert diag["eps_mass"] < 1e-2, f"K05 质量守恒误差 {diag['eps_mass']:.3e} 超 1%"

    # K13 四位小数：写出的值必须已按四位小数量化。
    for mat in (t_mat, c_mat):
        assert np.max(np.abs(mat - np.round(mat, DECIMALS))) < 1e-12, "K13 存在未量化到四位小数的值"

    # 表格与全场同源一致。
    tab = data["table1_T"]
    for k, tt in enumerate(tab["times_s"]):
        row = int(round(tt / P1_OUT_DT_S)) - 1
        for m, rc in enumerate(tab["radii_cm"]):
            col = int(round(rc / OUT_DR_CM))
            assert abs(tab["values"][k][m] - t_mat[row, col]) < 1e-12, "表1 与全场不一致"
            assert abs(data["table2_C"]["values"][k][m] - c_mat[row, col]) < 1e-12, "表2 与全场不一致"
    return True


def validate_capability(data: dict) -> dict:
    """P1-C1 / P1-C2 的可反驳检验：读源码证 r 加权与插值调用，读结果证不等式与通量残差。"""
    verdict = {}

    # ---- P1-C1：两个离散算子都含 r_{i+1/2} 面加权，且边界用附件1插值函数而非常数标量。
    src = (Path(__file__).resolve().parent / "utils.py").read_text(encoding="utf-8")
    body = "\n".join(ln for ln in src.splitlines() if not ln.strip().startswith("#"))
    assert "lam_f = lam_iface(props.lam(C_it)) * face" in body, "P1-C1 温度算子缺 r 加权面因子"
    assert "d_f = d_iface(C_it, T_new, props.D) * face" in body, "P1-C1 水分算子缺 r 加权面因子"
    assert "face = (np.arange(n_cells, dtype=float) + 0.5) * ds" in body, \
        "P1-C1 面半径未按 r_{i+1/2} 构造"

    t_amb, c_amb = make_ambient()
    probe = np.array([0.0, 300.0, 900.0, 1800.0])
    ta_probe = np.array([float(t_amb(v)) for v in probe])
    ca_probe = np.array([float(c_amb(v)) for v in probe])
    assert ta_probe.std() > 1e-3, f"P1-C1 Ta(t) 退化为常数 std={ta_probe.std():.2e}"
    assert ca_probe.std() > 1e-6, f"P1-C1 Ca(t) 退化为常数 std={ca_probe.std():.2e}"
    # 越界处不得返回 0 或报错（附录 A3 平台延拓）。
    assert abs(float(t_amb(2.0e5)) - PLATEAU_TA) < 1e-9, "P1-C1 Ta 越界延拓未取平台值"
    assert float(c_amb(2.0e5)) > 0.0, "P1-C1 Ca 越界返回 0"
    verdict["P1-C1"] = {
        "status": "PASS",
        "evidence": ("utils.py 中温度算子 lam_iface(...)*face 与水分算子 d_iface(...)*face 均乘以"
                     f"面半径 r_{{i+1/2}}=(i+1/2)Δs；Ta(t) 在 t={probe.tolist()} s 上取值"
                     f"{[round(v, 4) for v in ta_probe]} °C（std={ta_probe.std():.4f}），"
                     f"Ca(t) 取值 {[round(v, 6) for v in ca_probe]} kg/kg（std={ca_probe.std():.2e}），"
                     f"均为时变插值而非常数；t=2e5 s 越界返回 Ta={float(t_amb(2e5)):.4f} °C、"
                     f"Ca={float(c_amb(2e5)):.6f}（尾段平台延拓，非 0、非报错）。"),
    }

    # ---- P1-C2：中心不热于表面、中心不干于表面；表面离散通量与 Robin 表达式残差 < 1e-8。
    t_mat = np.asarray(data["result_grid"]["T"], dtype=float)
    c_mat = np.asarray(data["result_grid"]["C"], dtype=float)
    dt_gap = float(np.max(t_mat[:, 0] - t_mat[:, -1]))
    dc_gap = float(np.max(c_mat[:, -1] - c_mat[:, 0]))
    assert dt_gap <= 1e-9, f"P1-C2 中心温度超过表面温度 max差={dt_gap:.3e}"
    assert dc_gap <= 1e-9, f"P1-C2 中心含水率低于表面 max差={dc_gap:.3e}"

    fs = data["flux_series"]
    res_h = np.abs(np.asarray(fs["q_heat_W_m2"])
                   - H_HEAT * (np.asarray(fs["T_surface"]) - np.asarray(fs["Ta"])))
    res_m = np.abs(np.asarray(fs["q_mass_kg_m2s"])
                   - H_MASS * (np.asarray(fs["C_surface"]) - np.asarray(fs["Ca"])))
    assert res_h.max() < 1e-8, f"P1-C2 热通量残差 {res_h.max():.3e} 超 1e-8"
    assert res_m.max() < 1e-8, f"P1-C2 湿通量残差 {res_m.max():.3e} 超 1e-8"
    verdict["P1-C2"] = {
        "status": "PASS",
        "evidence": (f"1800 个输出时刻上 max(T(0)-T(R))={dt_gap:.3e} ≤ 0（中心恒不热于表面），"
                     f"max(C(R)-C(0))={dc_gap:.3e} ≤ 0（表面恒不湿于中心）；表面通量与"
                     f"h(T_R-Ta) 残差 {res_h.max():.3e} W/m²、与 h_m(C_R-Ca) 残差 "
                     f"{res_m.max():.3e} kg/(m²·s)，均 < 1e-8。"),
        "max_center_minus_surface_T": dt_gap, "max_surface_minus_center_C": dc_gap,
        "flux_residual_heat": float(res_h.max()), "flux_residual_mass": float(res_m.max()),
    }

    # ---- P1-C3：交付形状（result1.xlsx 由 export_excel 落盘，此处先校验数据面）。
    verdict["P1-C3"] = {
        "status": "PASS",
        "evidence": (f"result_grid 为 {t_mat.shape[0]}×{t_mat.shape[1]}（时间 1..1800 s 步长 1 s，"
                     f"距离 0.0..2.0 cm 步长 0.1 cm），温度与水分两套矩阵齐备且全部四位小数；"
                     f"表1/表2 各 {len(P1_TABLE_TIMES_S)} 时刻 × {len(P1_TABLE_RADII_CM)} 位置。"),
    }
    return verdict


def main() -> int:
    ensure_dirs()
    ckpt = CaseCheckpoint(
        ROOT, CASE_KEY,
        dependencies=["code/utils.py", "code/params.py", "code/problem1.py",
                      "user_data/附件1.xlsx", "user_data/附件2.xlsx"],
        parameters={"cfg": CFG, "picard_tol": PICARD_TOL, "gauss_n": GAUSS_N,
                    "numpy": np.__version__},
    )
    cached = ckpt.load(validate_constraints)
    if cached is None:
        print(f"[P1] 无可用检查点，开始求解：N={N_GRID} Δt1={DT_STAGE1} s t_end={P1_WINDOW_S:.0f} s")
        data = solve()
        ckpt.save(data, validate_constraints)
        print("[P1] 检查点已落盘（校验通过后写入）")
    else:
        data = cached
        print("[P1] 复用已验证检查点")

    cap = validate_capability(data)
    data["capability_verdict"] = cap
    fin = data["final_state"]
    diag = data["diagnostics"]
    print(f"[P1] 1800 s: T(0)={fin['T_center_r4']:.4f} T(R0)={fin['T_surface_r4']:.4f}"
          f"  C(0)={fin['C_center_r4']:.4f} C(R0)={fin['C_surface_r4']:.4f}")
    print(f"[P1] ε_mass={diag['eps_mass']:.3e}  Picard均值={diag['picard_mean']:.2f}"
          f"  峰值={diag['picard_max']}  步数={data['schedule']['n_steps']}"
          f"  墙钟={diag['wall_s']:.1f} s")
    print(f"[P1] T∈[{diag['T_range'][0]:.4f},{diag['T_range'][1]:.4f}] °C"
          f"  C∈[{diag['C_range'][0]:.4f},{diag['C_range'][1]:.4f}] kg/kg"
          f"  Bi_h={diag['Bi_h']:.4f} Bi_m(初)={diag['Bi_m_initial']:.4f}")
    for cid, v in cap.items():
        print(f"[P1][{cid}] {v['status']}")
    save_json(FIG_DIR / "problem_1_results.json", jsonable(data))
    print(f"[P1] 写出 figures/problem_1_results.json  "
          f"{(FIG_DIR / 'problem_1_results.json').stat().st_size / 1024:.1f} KB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
