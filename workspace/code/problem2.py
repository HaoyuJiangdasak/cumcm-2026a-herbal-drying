# -*- coding: utf-8 -*-
"""问题 2：变物性（附录 3）耦合模型，0—10800 s（3 h）温度场与水分场。

产出 figures/problem_2_results.json：逐秒 21 列全场（供 result2.xlsx）、论文表 3/表 4、
0.5 h 最大径向温差、双轴与山脊图字段，以及 K01—K06/K13 与 P2-C1/P2-C2 的运行期证据。
运行：python _utils/run_compute.py code/problem2.py
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
CASE_KEY = "2026A/problem2/appendix3/N40/dt0.01"
CFG = {"appendix": APPENDIX_OF_PROBLEM[2], "N": N_GRID, "out_dt": P2_OUT_DT_S,
       "t_end": P2_WINDOW_S, "dt1": DT_STAGE1, "dt2": DT_STAGE2, "dt3": DT_STAGE3,
       "record": "full", "tag": "P2-prod"}
HALF_HOUR_S = 0.5 * SEC_PER_HOUR


def solve() -> dict:
    """推进 0→10800 s，逐秒记录全场；表 3/表 4 取 P2_TABLE_TIMES_H 的 6 个时刻。"""
    res = run_case(CFG)
    t_out = np.asarray(res["t"], dtype=float)
    t_field = np.asarray(res["T_field"], dtype=float)
    c_field = np.asarray(res["C_field"], dtype=float)
    assert t_field.shape == (len(t_out), N_GRID + 1), f"温度场形状异常 {t_field.shape}"

    t_cols = node_columns(t_field, N_GRID)
    c_cols = node_columns(c_field, N_GRID)

    # result2.xlsx 主体：t=1..10800 s 步长 1 s（题面「每隔 1 s」，窗口取解释 (a)）。
    keep = t_out >= P2_OUT_DT_S - 1e-9
    grid_t = t_out[keep]
    assert len(grid_t) == 10800 and abs(grid_t[-1] - P2_WINDOW_S) < 1e-9, \
        f"result2 行数异常 n={len(grid_t)}"

    times_tab_s = [th * SEC_PER_HOUR for th in P2_TABLE_TIMES_H]
    idx_tab = [pick_time_index(t_out, tt) for tt in times_tab_s]
    col_tab = [int(round(rc / OUT_DR_CM)) for rc in P2_TABLE_RADII_CM]
    table3 = [[round_out(t_cols[i, j]) for j in col_tab] for i in idx_tab]
    table4 = [[round_out(c_cols[i, j]) for j in col_tab] for i in idx_tab]

    # 径向温差历程：每个输出时刻的 T(R)-T(0)，并定位全窗口峰值与 0.5 h 处的值。
    dt_radial = t_cols[keep][:, -1] - t_cols[keep][:, 0]
    i_peak = int(np.argmax(dt_radial))
    i_half = int(round(HALF_HOUR_S / P2_OUT_DT_S)) - 1
    assert abs(grid_t[i_half] - HALF_HOUR_S) < 1e-9, "0.5 h 时刻未对齐输出网格"

    t_amb, c_amb = make_ambient()
    ta_series = np.asarray([t_amb(tt) for tt in grid_t], dtype=float)
    ca_series = np.asarray([c_amb(tt) for tt in grid_t], dtype=float)

    props = Props(CFG["appendix"])
    ridge_idx = [pick_time_index(t_out, tt) for tt in times_tab_s]
    dual = downsample_index(len(grid_t), 400)
    out = {
        "problem": 2, "appendix": CFG["appendix"], "tag": CFG["tag"],
        "window_interpretation": "a: [1,10800] s（MODELING_REPORT 第 22 节裁定）",
        "grid": {"N": N_GRID, "ds": float(make_grid(N_GRID)[1]), "out_stride": OUT_STRIDE,
                 "radii_cm": list(OUT_RADII_CM)},
        "schedule": {"dt_stage1_s": DT_STAGE1, "dt_stage2_s": DT_STAGE2,
                     "out_dt_s": P2_OUT_DT_S, "t_end_s": P2_WINDOW_S,
                     "n_steps": int(res["n_steps"])},
        "result_grid": {"times_s": [float(v) for v in grid_t], "radii_cm": list(OUT_RADII_CM),
                        "T": [[round_out(v) for v in row] for row in t_cols[keep]],
                        "C": [[round_out(v) for v in row] for row in c_cols[keep]]},
        "table3_T": {"times_h": list(P2_TABLE_TIMES_H), "times_s": times_tab_s,
                     "radii_cm": list(P2_TABLE_RADII_CM), "values": table3},
        "table4_C": {"times_h": list(P2_TABLE_TIMES_H), "times_s": times_tab_s,
                     "radii_cm": list(P2_TABLE_RADII_CM), "values": table4},
        "mass_series": {"times_s": [float(v) for v in grid_t],
                        "C_bar": [float(v) for v in np.asarray(res["C_bar"])[keep]],
                        "C_max": [float(v) for v in np.asarray(res["C_max"])[keep]],
                        "flux_cum": [float(v) for v in np.asarray(res["flux_cum"])[keep]]},
        "radial_gap": {"times_s": [float(v) for v in grid_t],
                       "dT_surface_minus_center": [float(v) for v in dt_radial],
                       "peak_dT": float(dt_radial[i_peak]), "peak_t_s": float(grid_t[i_peak]),
                       "peak_t_h": float(grid_t[i_peak] / SEC_PER_HOUR),
                       "at_half_hour": {"t_s": float(grid_t[i_half]),
                                        "dT": float(dt_radial[i_half]),
                                        "T_center": round_out(t_cols[keep][i_half, 0]),
                                        "T_surface": round_out(t_cols[keep][i_half, -1])}},
        "final_state": {"t_s": float(t_out[-1]), "t_h": float(t_out[-1] / SEC_PER_HOUR),
                        "T_center": float(t_field[-1, 0]), "T_surface": float(t_field[-1, -1]),
                        "C_center": float(c_field[-1, 0]), "C_surface": float(c_field[-1, -1]),
                        "T_center_r4": round_out(t_field[-1, 0]),
                        "T_surface_r4": round_out(t_field[-1, -1]),
                        "C_center_r4": round_out(c_field[-1, 0]),
                        "C_surface_r4": round_out(c_field[-1, -1])},
        "dual_axis": {"times_s": [float(grid_t[i]) for i in dual],
                      "Ta": [float(ta_series[i]) for i in dual],
                      "Ca": [float(ca_series[i]) for i in dual],
                      "T_center": [float(t_cols[keep][i, 0]) for i in dual],
                      "T_surface": [float(t_cols[keep][i, -1]) for i in dual],
                      "C_center": [float(c_cols[keep][i, 0]) for i in dual],
                      "C_surface": [float(c_cols[keep][i, -1]) for i in dual],
                      "C_bar": [float(res["C_bar"][keep][i]) for i in dual]},
        "ridgeline": {"times_s": times_tab_s, "times_h": list(P2_TABLE_TIMES_H),
                      "radii_cm": list(OUT_RADII_CM),
                      "T": [[float(v) for v in t_cols[i]] for i in ridge_idx],
                      "C": [[float(v) for v in c_cols[i]] for i in ridge_idx]},
        "property_evolution": {
            "times_s": [float(grid_t[i]) for i in dual],
            "D_center": [float(res["D_center"][keep][i]) for i in dual],
            "D_surface": [float(res["D_surface"][keep][i]) for i in dual],
            "rho_center": [float(props.rho(np.array([c_cols[keep][i, 0]]))[0]) for i in dual],
            "rho_surface": [float(props.rho(np.array([c_cols[keep][i, -1]]))[0]) for i in dual],
            "cp_surface": [float(props.cp(np.array([c_cols[keep][i, -1]]))[0]) for i in dual],
            "lam_surface": [float(props.lam(np.array([c_cols[keep][i, -1]]))[0]) for i in dual],
            "Bi_m_surface": [float(props.bi_mass(np.array([c_cols[keep][i, -1]]),
                                                 np.array([t_cols[keep][i, -1]]))[0])
                             for i in dual]},
        "diagnostics": {"eps_mass": float(res["eps_mass"]), "m0": float(res["m0"]),
                        "m_end": float(res["m_end"]), "flux_total": float(res["flux_total"]),
                        "picard_mean": float(res["picard_mean"]),
                        "picard_max": int(res["picard_max"]), "wall_s": float(res["wall_s"]),
                        "T_range": [float(t_field.min()), float(t_field.max())],
                        "C_range": [float(c_field.min()), float(c_field.max())],
                        "Ta_max_used": float(ta_series.max()),
                        "D_center_range": [float(res["D_center"].min()),
                                           float(res["D_center"].max())],
                        "D_surface_range": [float(res["D_surface"].min()),
                                            float(res["D_surface"].max())],
                        "C_bar_end": float(res["C_bar"][-1])},
    }
    return out


def validate_constraints(data: dict) -> bool:
    """K01—K06、K13：从写出的结果矩阵复核，不重跑求解器。"""
    t_mat = np.asarray(data["result_grid"]["T"], dtype=float)
    c_mat = np.asarray(data["result_grid"]["C"], dtype=float)
    times = np.asarray(data["result_grid"]["times_s"], dtype=float)
    diag = data["diagnostics"]

    assert t_mat.shape == (10800, 21) and c_mat.shape == (10800, 21), \
        f"K13 表格形状 {t_mat.shape}/{c_mat.shape} 不是 10800x21"
    assert np.allclose(np.diff(times), P2_OUT_DT_S), "时间列步长非 1 s"
    assert abs(times[0] - 1.0) < 1e-9 and abs(times[-1] - P2_WINDOW_S) < 1e-9, "时间列端点错误"

    # K01 含水率逐点落在 [0, 2.55]。
    assert c_mat.min() >= -1e-9 and c_mat.max() <= C_INIT + 1e-9, \
        f"K01 含水率越界 [{c_mat.min()},{c_mat.max()}]"

    # K02 温度逐点落在 [28, max_t Ta + 0.5]。
    hi = diag["Ta_max_used"] + 0.5
    assert t_mat.min() >= T_INIT - 1e-9, f"K02 温度低于初温 {t_mat.min()}"
    assert t_mat.max() <= hi, f"K02 温度超上界 {t_mat.max()} > {hi}"

    # K03 ∂C/∂r ≤ 0，等价于 max_r C = C(0,t)。
    dc_r = float(np.max(np.diff(c_mat, axis=1)))
    assert dc_r <= 1e-9, f"K03 含水率沿半径出现上升 max差={dc_r:.3e}"
    assert np.allclose(c_mat.max(axis=1), c_mat[:, 0], atol=1e-12), "K03 逐时刻最大值不在中心"

    # K04 平均含水率单调不增；K05 全局质量守恒。
    c_bar = np.asarray(data["mass_series"]["C_bar"], dtype=float)
    d_bar = float(np.max(np.diff(c_bar)))
    assert d_bar <= 1e-12, f"K04 平均含水率出现回升 max差={d_bar:.3e}"
    assert diag["eps_mass"] < 1e-2, f"K05 质量守恒误差 {diag['eps_mass']:.3e} 超 1%"

    for mat in (t_mat, c_mat):
        assert np.max(np.abs(mat - np.round(mat, DECIMALS))) < 1e-12, "K13 存在未量化到四位小数的值"

    for key, mat in (("table3_T", t_mat), ("table4_C", c_mat)):
        tab = data[key]
        for k, tt in enumerate(tab["times_s"]):
            row = int(round(tt / P2_OUT_DT_S)) - 1
            for m, rc in enumerate(tab["radii_cm"]):
                col = int(round(rc / OUT_DR_CM))
                assert abs(tab["values"][k][m] - mat[row, col]) < 1e-12, f"{key} 与全场不一致"

    # 变物性必须真的在变（否则等价于常物性模型）：3 h 窗口内表面 D 的量程比登记备查。
    lo, hi_d = diag["D_surface_range"]
    assert hi_d / max(lo, 1e-300) > 1.2, f"表面 D 未随变干下降 ratio={hi_d / lo:.3f}"
    return True


def validate_capability(data: dict) -> dict:
    """P2-C1 / P2-C2 的可反驳检验：物性数组量重算、界面平均、平台延拓、开尔文与附录隔离。"""
    verdict = {}
    src_dir = Path(__file__).resolve().parent
    body = "\n".join(ln for ln in (src_dir / "utils.py").read_text(encoding="utf-8").splitlines()
                     if not ln.strip().startswith("#"))

    # ---- P2-C1：物性每层由当前 C、T 场重算（数组），界面 D 取相邻控制体平均。
    assert "props.lam(C_it)" in body and "props.rho(C_it)" in body, "P2-C1 物性未用当前 C 场重算"
    assert "d_iface(C_it, T_new, props.D)" in body, "P2-C1 界面 D 未由相邻控制体值平均"
    props = Props(APPENDIX_OF_PROBLEM[2])
    c_probe = np.array([C_INIT, 1.0, C_THRESHOLD])
    t_probe = np.array([T_INIT, 40.0, PLATEAU_TA])
    rho_v, cp_v, lam_v = props.rho(c_probe), props.cp(c_probe), props.lam(c_probe)
    d_v = props.D(c_probe, t_probe)
    for name, arr in (("rho", rho_v), ("cp", cp_v), ("lam", lam_v), ("D", d_v)):
        assert isinstance(arr, np.ndarray) and arr.shape == c_probe.shape, \
            f"P2-C1 {name} 未返回与场同形的数组"
        assert float(np.ptp(arr)) > 0.0, f"P2-C1 {name} 退化为常数标量"
    # 界面平均必须落在两端点之间（积分平均的必要条件）。
    from utils import d_iface
    c_pair = np.array([C_INIT, 1.0])
    t_pair = np.array([PLATEAU_TA, PLATEAU_TA])
    iface = float(d_iface(c_pair, t_pair, props.D)[0])
    d_lo = float(props.D(np.array([1.0]), np.array([PLATEAU_TA]))[0])
    d_hi = float(props.D(np.array([C_INIT]), np.array([PLATEAU_TA]))[0])
    assert d_lo < iface < d_hi, f"P2-C1 界面平均 {iface:.4e} 未落在两端点 ({d_lo:.4e},{d_hi:.4e}) 之间"

    t_amb, c_amb = make_ambient()
    beyond = ATT1_T_END_S + 1000.0
    ta_beyond, ca_beyond = float(t_amb(beyond)), float(c_amb(beyond))
    assert abs(ta_beyond - PLATEAU_TA) < 1e-9, f"P2-C1 t>14400 s 未返回平台实测值 {ta_beyond}"
    assert abs(ca_beyond - PLATEAU_CA) < 1e-12 and ca_beyond > 0.0, \
        f"P2-C1 t>14400 s 的 C_a 延拓错误 {ca_beyond}"
    verdict["P2-C1"] = {
        "status": "PASS",
        "evidence": (f"utils.py 每个 Picard 层调用 props.rho/cp/lam(C_it)、props.D(C_it,T_new) 重算物性；"
                     f"在 C={c_probe.tolist()} 上 ρ={np.round(rho_v, 2).tolist()}、"
                     f"c_p={np.round(cp_v, 1).tolist()}、λ={np.round(lam_v, 4).tolist()}、"
                     f"D={['%.3e' % v for v in d_v]}，均为随场变化的数组量而非常数；"
                     f"界面 D 由 d_iface 对相邻控制体 (C_i,C_{{i+1}})=(2.55,1.00) 积分平均得 "
                     f"{iface:.4e}，严格落在端点值 ({d_lo:.4e}, {d_hi:.4e}) 之间；"
                     f"t={beyond:.0f} s（>14400 s）激励返回 T_a={ta_beyond:.4f} °C、"
                     f"C_a={ca_beyond:.6f}（附件 1 尾段实测平台均值，非 0、非报错）。"),
    }

    # ---- P2-C2（含 K06、K07）：调用 D 前断言 T_K>273；附录 3 不得混入附录 2/4 常数。
    assert "assert np.all(t_k > 273.0)" in body, "K06 缺开尔文正值断言"
    assert "to_kelvin" in body, "K06 未走开尔文换算"
    guard_fired = False
    try:
        props.D(np.array([1.0]), np.array([-300.0]))
    except AssertionError:
        guard_fired = True   # 反证探针：故意传 <0 K，必须被 K06 断言拦住
    assert guard_fired, "K06 传入 -300 °C（低于 0 K）竟未触发开尔文断言"

    # K07 源码扫描：问题 2/3 的物性分支不得出现附录 2、附录 4 的常数字面量。
    props._check_isolation()
    rho3 = float(props.rho(np.array([0.0]))[0])
    cp3 = float(props.cp(np.array([0.0]))[0])
    lam3 = float(props.lam(np.array([0.0]))[0])
    assert abs(rho3 - RHO_A2) > 1.0 and abs(rho3 - RHO4_A) > 1.0, "K07 密度串用了其他附录"
    assert abs(cp3 - CP_A2) > 1.0 and abs(cp3 - CP4_A) > 1.0, "K07 比热容串用了其他附录"
    assert abs(lam3 - LAM_A2) > 1e-3 and abs(lam3 - LAM4_A) > 1e-3, "K07 热导率串用了其他附录"
    verdict["P2-C2"] = {
        "status": "PASS",
        "evidence": (f"utils.to_kelvin 内含 assert np.all(t_k > 273.0)，扩散系数只经该函数取绝对温度；"
                     f"反证探针传入 −300 °C（−26.85 K）时断言按预期抛出。附录隔离实测：Props(3) 在 C→0 "
                     f"给出 ρ={rho3:.3f}、c_p={cp3:.3f}、λ={lam3:.5f}，与附录 2 的 "
                     f"({RHO_A2:.1f}, {CP_A2:.1f}, {LAM_A2:.4f})、附录 4 的 "
                     f"({RHO4_A:.1f}, {CP4_A:.1f}, {LAM4_A:.4f}) 均显著不同，"
                     f"Props._check_isolation() 通过。"),
        "rho_appendix3_at_C0": rho3, "cp_appendix3_at_C0": cp3, "lam_appendix3_at_C0": lam3,
    }

    # ---- P2-C3：交付形状。
    t_mat = np.asarray(data["result_grid"]["T"], dtype=float)
    verdict["P2-C3"] = {
        "status": "PASS",
        "evidence": (f"result_grid 为 {t_mat.shape[0]}×{t_mat.shape[1]}（每隔 1 s × 每隔 0.1 cm，"
                     f"窗口按第 22 节裁定取 [1,10800] s），温度与水分双矩阵齐备且全部四位小数；"
                     f"表 3/表 4 各 {len(P2_TABLE_TIMES_H)} 时刻 × {len(P2_TABLE_RADII_CM)} 位置。"),
    }
    return verdict


def main() -> int:
    ensure_dirs()
    ckpt = CaseCheckpoint(
        ROOT, CASE_KEY,
        dependencies=["code/utils.py", "code/params.py", "code/problem2.py",
                      "user_data/附件1.xlsx", "user_data/附件2.xlsx"],
        parameters={"cfg": CFG, "picard_tol": PICARD_TOL, "gauss_n": GAUSS_N,
                    "numpy": np.__version__},
    )
    cached = ckpt.load(validate_constraints)
    if cached is None:
        print(f"[P2] 无可用检查点，开始求解：附录 {CFG['appendix']} N={N_GRID} "
              f"t_end={P2_WINDOW_S:.0f} s")
        data = solve()
        ckpt.save(data, validate_constraints)
        print("[P2] 检查点已落盘（校验通过后写入）")
    else:
        data = cached
        print("[P2] 复用已验证检查点")

    cap = validate_capability(data)
    data["capability_verdict"] = cap
    fin, diag, gap = data["final_state"], data["diagnostics"], data["radial_gap"]
    print(f"[P2] 3 h: T(0)={fin['T_center_r4']:.4f} T(R0)={fin['T_surface_r4']:.4f}"
          f"  C(0)={fin['C_center_r4']:.4f} C(R0)={fin['C_surface_r4']:.4f}")
    print(f"[P2] 径向温差峰值={gap['peak_dT']:.4f} °C @ t={gap['peak_t_h']:.4f} h"
          f"；0.5 h 处 ΔT={gap['at_half_hour']['dT']:.4f} °C"
          f"（{gap['at_half_hour']['T_center']:.4f}→{gap['at_half_hour']['T_surface']:.4f}）")
    print(f"[P2] ε_mass={diag['eps_mass']:.3e}  Picard均值={diag['picard_mean']:.2f}"
          f"  峰值={diag['picard_max']}  步数={data['schedule']['n_steps']}"
          f"  墙钟={diag['wall_s']:.1f} s")
    print(f"[P2] D_surface∈[{diag['D_surface_range'][0]:.4e},{diag['D_surface_range'][1]:.4e}]"
          f"  D_center∈[{diag['D_center_range'][0]:.4e},{diag['D_center_range'][1]:.4e}]"
          f"  C̄(末)={diag['C_bar_end']:.4f}")
    for cid, v in cap.items():
        print(f"[P2][{cid}] {v['status']}")
    save_json(FIG_DIR / "problem_2_results.json", jsonable(data))
    print(f"[P2] 写出 figures/problem_2_results.json  "
          f"{(FIG_DIR / 'problem_2_results.json').stat().st_size / 1024:.1f} KB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
