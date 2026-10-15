# -*- coding: utf-8 -*-
"""问题 3：附录 3 物性、固定半径，求最早全域达标时刻 t_dry = inf{t : max_r C < 0.15}。

产出 figures/problem_3_results.json：每 60 s × 0.1 cm 的水分场（供 result3.xlsx）、论文表 5
（每 6 h × 0.5 cm，末行为 t_dry）、判据历程与跨阈值区间，以及 K01—K05/K10—K13 与
P3-C1/P3-C2 的运行期证据。
运行：python _utils/run_compute.py code/problem3.py
"""
from __future__ import annotations

import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_utils"))
from compute_checkpoint import CaseCheckpoint  # noqa: E402

from params import *  # noqa: F401,F403,E402
from utils import (Props, downsample_index, ensure_dirs, jsonable, make_ambient,  # noqa: E402
                   make_grid, node_columns, round_out, run_case, save_json)

ROOT = Path(__file__).resolve().parents[1]
CASE_KEY = "2026A/problem3/appendix3/N40/stop0.15"
CFG = {"appendix": APPENDIX_OF_PROBLEM[3], "N": N_GRID, "out_dt": P3_OUT_DT_S,
       "t_end": T_MAX_SEARCH_S, "stop_threshold": C_THRESHOLD, "dt1": DT_STAGE1,
       "dt2": DT_STAGE2, "dt3": DT_STAGE3, "record": "full", "tag": "P3-prod"}


def solve() -> dict:
    """推进到 max_r C 首次跌破 0.15；输出步长 60 s，t_dry 用跨阈值线性插值。"""
    res = run_case(CFG)
    t_out = np.asarray(res["t"], dtype=float)
    c_field = np.asarray(res["C_field"], dtype=float)
    t_field = np.asarray(res["T_field"], dtype=float)
    t_dry = res["t_dry"]
    assert t_dry is not None, "未在搜索窗内达标，需放宽 T_MAX_SEARCH_S"
    bracket = [float(v) for v in res["bracket"]]

    c_cols = node_columns(c_field, N_GRID)
    t_cols = node_columns(t_field, N_GRID)
    c_max = np.asarray(res["C_max"], dtype=float)
    c_bar = np.asarray(res["C_bar"], dtype=float)

    # result3.xlsx 主体：t=60,120,... 直到烘干结束时间（含首次跌破的那一行）。
    keep = t_out >= P3_OUT_DT_S - 1e-9
    grid_t = t_out[keep]
    assert abs(grid_t[0] - P3_OUT_DT_S) < 1e-9, "输出起点应为 60 s"
    assert abs(grid_t[-1] - bracket[1]) < 1e-9, "末行应为首次跌破 0.15 的输出时刻"

    # 论文表 5：每 6 h × 0.5 cm，末行补 t_dry（该时刻用跨阈值两端线性插值得到剖面）。
    step_s = P3_TABLE_DT_H * SEC_PER_HOUR
    tab_times = [k * step_s for k in range(1, int(grid_t[-1] // step_s) + 1)]
    col_tab = list(range(0, N_OUT_RADII, int(round(P3_TABLE_DR_CM / OUT_DR_CM))))
    rows, row_times = [], []
    for tt in tab_times:
        i = int(round(tt / P3_OUT_DT_S)) - 1
        assert abs(grid_t[i] - tt) < 1e-9, f"表 5 时刻 {tt} 未对齐输出网格"
        rows.append([round_out(c_cols[keep][i, j]) for j in col_tab])
        row_times.append(float(tt))
    w = (t_dry - bracket[0]) / (bracket[1] - bracket[0])
    i_lo = int(round(bracket[0] / P3_OUT_DT_S)) - 1
    prof_dry = c_cols[keep][i_lo] * (1.0 - w) + c_cols[keep][i_lo + 1] * w
    rows.append([round_out(prof_dry[j]) for j in col_tab])
    row_times.append(float(t_dry))

    t_amb, c_amb = make_ambient()
    dsi = downsample_index(len(grid_t), 600)
    props = Props(CFG["appendix"])
    out = {
        "problem": 3, "appendix": CFG["appendix"], "tag": CFG["tag"],
        "criterion": {"expr": "t_dry = inf{t : max_r C(r,t) < 0.15}",
                      "threshold": C_THRESHOLD, "scope": CRITERION_SCOPE,
                      "applied_to": "径向数组的最大值（等价于中心值 C(0,t)）"},
        "grid": {"N": N_GRID, "ds": float(make_grid(N_GRID)[1]), "out_stride": OUT_STRIDE,
                 "radii_cm": list(OUT_RADII_CM)},
        "schedule": {"out_dt_s": P3_OUT_DT_S, "dt_stage1_s": DT_STAGE1,
                     "dt_stage2_s": DT_STAGE2, "dt_stage3_s": DT_STAGE3,
                     "t_search_limit_s": T_MAX_SEARCH_S, "n_steps": int(res["n_steps"])},
        "t_dry": {"seconds": float(t_dry), "hours": float(t_dry / SEC_PER_HOUR),
                  "hours_r4": round_out(t_dry / SEC_PER_HOUR),
                  "bracket_t_lo_s": bracket[0], "bracket_t_hi_s": bracket[1],
                  "bracket_Cmax_lo": bracket[2], "bracket_Cmax_hi": bracket[3],
                  "C_bar_at_t_dry": float(c_bar[-1]), "C_max_at_t_dry": float(c_max[-1]),
                  "C_center_at_t_dry": round_out(prof_dry[0]),
                  "C_surface_at_t_dry": round_out(prof_dry[-1])},
        "result_grid": {"times_s": [float(v) for v in grid_t], "radii_cm": list(OUT_RADII_CM),
                        "C": [[round_out(v) for v in row] for row in c_cols[keep]]},
        "table5_C": {"times_s": row_times,
                     "times_h": [float(v / SEC_PER_HOUR) for v in row_times],
                     "radii_cm": [OUT_RADII_CM[j] for j in col_tab], "values": rows,
                     "last_row_is_t_dry": True},
        "criterion_series": {"times_s": [float(grid_t[i]) for i in dsi],
                             "times_h": [float(grid_t[i] / SEC_PER_HOUR) for i in dsi],
                             "C_max": [float(c_max[keep][i]) for i in dsi],
                             "C_bar": [float(c_bar[keep][i]) for i in dsi],
                             "C_surface": [float(c_cols[keep][i, -1]) for i in dsi],
                             "T_center": [float(t_cols[keep][i, 0]) for i in dsi],
                             "T_surface": [float(t_cols[keep][i, -1]) for i in dsi],
                             "Ta": [float(t_amb(grid_t[i])) for i in dsi],
                             "Ca": [float(c_amb(grid_t[i])) for i in dsi],
                             "D_center": [float(res["D_center"][keep][i]) for i in dsi],
                             "D_surface": [float(res["D_surface"][keep][i]) for i in dsi],
                             "Bi_m_surface": [
                                 float(props.bi_mass(np.array([c_cols[keep][i, -1]]),
                                                     np.array([t_cols[keep][i, -1]]))[0])
                                 for i in dsi]},
        "mass_series": {"times_s": [float(grid_t[i]) for i in dsi],
                        "flux_cum": [float(res["flux_cum"][keep][i]) for i in dsi]},
        "final_state": {"t_s": float(t_out[-1]), "t_h": float(t_out[-1] / SEC_PER_HOUR),
                        "C_center": float(c_field[-1, 0]), "C_surface": float(c_field[-1, -1]),
                        "T_center": float(t_field[-1, 0]), "T_surface": float(t_field[-1, -1]),
                        "C_center_r4": round_out(c_field[-1, 0]),
                        "C_surface_r4": round_out(c_field[-1, -1]),
                        "T_center_r4": round_out(t_field[-1, 0]),
                        "T_surface_r4": round_out(t_field[-1, -1])},
        "diagnostics": {"eps_mass": float(res["eps_mass"]), "m0": float(res["m0"]),
                        "m_end": float(res["m_end"]), "flux_total": float(res["flux_total"]),
                        "picard_mean": float(res["picard_mean"]),
                        "picard_max": int(res["picard_max"]), "wall_s": float(res["wall_s"]),
                        "C_range": [float(c_field.min()), float(c_field.max())],
                        "T_range": [float(t_field.min()), float(t_field.max())],
                        "Ta_max_used": float(max(t_amb(v) for v in grid_t)),
                        "D_surface_range": [float(res["D_surface"].min()),
                                            float(res["D_surface"].max())]},
    }
    bi = out["criterion_series"]["Bi_m_surface"]
    out["diagnostics"]["Bi_m_range"] = [float(min(bi)), float(max(bi))]
    return out


def validate_constraints(data: dict) -> bool:
    """K01—K05、K10—K13：从写出的结果复核。"""
    c_mat = np.asarray(data["result_grid"]["C"], dtype=float)
    times = np.asarray(data["result_grid"]["times_s"], dtype=float)
    diag, td = data["diagnostics"], data["t_dry"]

    assert c_mat.shape[1] == 21, f"K13 列数 {c_mat.shape[1]} 不是 21"
    assert np.allclose(np.diff(times), P3_OUT_DT_S), "时间列步长非 60 s"
    assert abs(times[0] - P3_OUT_DT_S) < 1e-9, "时间列起点非 60 s"
    assert np.max(np.abs(c_mat - np.round(c_mat, DECIMALS))) < 1e-12, "K13 未量化到四位小数"
    assert not np.isnan(c_mat).any(), "K13 存在空缺"

    assert c_mat.min() >= -1e-9 and c_mat.max() <= C_INIT + 1e-9, \
        f"K01 含水率越界 [{c_mat.min()},{c_mat.max()}]"
    lo_t, hi_t = diag["T_range"]
    assert lo_t >= T_INIT - 1e-9 and hi_t <= diag["Ta_max_used"] + 0.5, \
        f"K02 温度越界 [{lo_t},{hi_t}]"

    dc_r = float(np.max(np.diff(c_mat, axis=1)))
    assert dc_r <= 1e-9, f"K03 含水率沿半径出现上升 max差={dc_r:.3e}"
    assert np.allclose(c_mat.max(axis=1), c_mat[:, 0], atol=1e-12), "K03 逐时刻最大值不在中心"

    c_bar = np.asarray(data["criterion_series"]["C_bar"], dtype=float)
    assert float(np.max(np.diff(c_bar))) <= 1e-12, "K04 平均含水率出现回升"
    assert diag["eps_mass"] < 1e-2, f"K05 质量守恒误差 {diag['eps_mass']:.3e} 超 1%"

    # K11 跨阈值区间：两端必须真正跨过 0.15，且 t_dry 落在区间内。
    assert td["bracket_Cmax_lo"] >= C_THRESHOLD > td["bracket_Cmax_hi"], \
        f"K11 区间未跨阈值 {td['bracket_Cmax_lo']} / {td['bracket_Cmax_hi']}"
    assert td["bracket_t_lo_s"] <= td["seconds"] <= td["bracket_t_hi_s"], "K11 t_dry 不在区间内"
    assert td["bracket_t_hi_s"] - td["bracket_t_lo_s"] <= P3_OUT_DT_S + 1e-9, "K11 区间宽于一个输出步"

    # K10 达标时刻须落在附件 2 覆盖窗内（问题 4 需要，问题 3 同步登记）。
    assert td["seconds"] < ATT2_T_END_S, f"K10 t_dry={td['seconds']} 超出附件 2 覆盖"

    # 逐点判据与均值判据必须可区分：t_dry 处 max_r C 刚跌破 0.15 而 C̄ 明显更低。
    assert td["C_max_at_t_dry"] < C_THRESHOLD, "t_dry 处 max_r C 未跌破阈值"
    assert td["C_bar_at_t_dry"] < C_THRESHOLD - 1e-3, "t_dry 处 C̄ 未明显低于阈值"

    # 表 5 末行须为 t_dry，且各行与全场同源。
    tab = data["table5_C"]
    assert tab["last_row_is_t_dry"] and abs(tab["times_s"][-1] - td["seconds"]) < 1e-9, \
        "表 5 末行不是烘干结束时间"
    step_s = P3_TABLE_DT_H * SEC_PER_HOUR
    for k, tt in enumerate(tab["times_s"][:-1]):
        assert abs(tt - (k + 1) * step_s) < 1e-9, "表 5 时刻步长非 6 h"
        row = int(round(tt / P3_OUT_DT_S)) - 1
        for m, rc in enumerate(tab["radii_cm"]):
            col = int(round(rc / OUT_DR_CM))
            assert abs(tab["values"][k][m] - c_mat[row, col]) < 1e-12, "表 5 与全场不一致"
    return True


def validate_capability(data: dict) -> dict:
    """P3-C1 / P3-C2 的可反驳检验：判据作用于径向最大值，且 t_dry 严格落在跨阈值区间。"""
    verdict = {}
    body = "\n".join(ln for ln in (Path(__file__).resolve().parent / "utils.py")
                     .read_text(encoding="utf-8").splitlines() if not ln.strip().startswith("#"))

    # ---- P3-C1：终止判断必须是 max_r C，不是均值、不是任意单点。
    assert "c_max_now = float(C.max())" in body, "P3-C1 终止量未取径向最大值"
    assert "if c_th is not None and c_max_now < c_th" in body, "P3-C1 终止判断未用 max_r C"
    stop_block = body.split("if c_th is not None and c_max_now < c_th")[1][:400]
    assert "rec_max" in stop_block and "rec_bar" not in stop_block, "P3-C1 终止定位混入均值序列"

    c_mat = np.asarray(data["result_grid"]["C"], dtype=float)
    td = data["t_dry"]
    cs = data["criterion_series"]
    # max_r C 等价于中心值：逐时刻核对（若判据错取表面值会立刻被证伪）。
    same_as_center = bool(np.allclose(c_mat.max(axis=1), c_mat[:, 0], atol=1e-12))
    surf_last = float(c_mat[-1, -1])
    assert same_as_center, "P3-C1 max_r C 不等于中心值"
    assert surf_last < td["C_max_at_t_dry"], "P3-C1 表面值应严格低于中心值"
    gap_bar = C_THRESHOLD - td["C_bar_at_t_dry"]
    assert gap_bar > 1e-3, "P3-C1 C̄ 与阈值差距过小，判据无法区分逐点与均值"
    verdict["P3-C1"] = {
        "status": "PASS",
        "evidence": (f"utils.run_case 的终止量为 c_max_now = float(C.max())，判断式 "
                     f"`if c_th is not None and c_max_now < c_th`，作用于径向数组最大值；"
                     f"逐时刻核对 max_r C 恒等于 C(0,t)（allclose atol=1e-12）。"
                     f"在 t_dry={td['hours']:.4f} h 处 max_r C={td['C_max_at_t_dry']:.6f} 刚跌破 0.15，"
                     f"而 C̄={td['C_bar_at_t_dry']:.6f}（低于阈值 {gap_bar:.6f}），"
                     f"表面 C={surf_last:.6f}；两者并存，说明若改用均值判据会提前 "
                     f"{_mean_criterion_lead(cs):.4f} h 误判达标。"),
        "C_max_at_t_dry": td["C_max_at_t_dry"], "C_bar_at_t_dry": td["C_bar_at_t_dry"],
    }

    # ---- P3-C2：跨阈值区间存在且 t_dry 落在其中（不早于首次跌破、不晚于其后一个输出步）。
    lo_s, hi_s = td["bracket_t_lo_s"], td["bracket_t_hi_s"]
    assert td["bracket_Cmax_lo"] >= C_THRESHOLD > td["bracket_Cmax_hi"], "P3-C2 区间未跨阈值"
    assert lo_s < td["seconds"] < hi_s, "P3-C2 t_dry 未严格落在区间内部"
    assert abs((hi_s - lo_s) - P3_OUT_DT_S) < 1e-9, "P3-C2 区间宽度不等于一个输出步"
    verdict["P3-C2"] = {
        "status": "PASS",
        "evidence": (f"存在输出时刻对 (t_k, t_{{k+1}}) = ({lo_s:.0f}, {hi_s:.0f}) s，"
                     f"max_r C(t_k)={td['bracket_Cmax_lo']:.6f} ≥ 0.15 > "
                     f"max_r C(t_{{k+1}})={td['bracket_Cmax_hi']:.6f}；线性插值得 "
                     f"t_dry={td['seconds']:.3f} s = {td['hours']:.4f} h，严格落在该区间内部，"
                     f"区间宽度恰为一个输出步 {P3_OUT_DT_S:.0f} s。"),
        "bracket_s": [lo_s, hi_s],
    }

    verdict["P3-C4"] = {
        "status": "PASS",
        "evidence": (f"result_grid 为 {c_mat.shape[0]}×{c_mat.shape[1]}（每隔 60 s × 每隔 0.1 cm，"
                     f"直到烘干结束时刻 {hi_s:.0f} s），全部四位小数无空缺；表 5 为 "
                     f"{len(data['table5_C']['times_s'])} 行（每 6 h 一行，末行为 t_dry="
                     f"{td['hours']:.4f} h）× {len(data['table5_C']['radii_cm'])} 列（每 0.5 cm）。"),
    }
    return verdict


def _mean_criterion_lead(cs: dict) -> float:
    """若误用均值判据，会提前多少小时判达标（供 P3-C1 的可区分性证据）。"""
    bar = np.asarray(cs["C_bar"], dtype=float)
    hrs = np.asarray(cs["times_h"], dtype=float)
    hit = np.nonzero(bar < C_THRESHOLD)[0]
    if hit.size == 0:
        return 0.0
    return float(hrs[-1] - hrs[hit[0]])


def main() -> int:
    ensure_dirs()
    ckpt = CaseCheckpoint(
        ROOT, CASE_KEY,
        dependencies=["code/utils.py", "code/params.py", "code/problem3.py",
                      "user_data/附件1.xlsx", "user_data/附件2.xlsx"],
        parameters={"cfg": CFG, "picard_tol": PICARD_TOL, "gauss_n": GAUSS_N,
                    "numpy": np.__version__},
    )
    cached = ckpt.load(validate_constraints)
    if cached is None:
        print(f"[P3] 无可用检查点，开始求解：附录 {CFG['appendix']} N={N_GRID} "
              f"阈值={C_THRESHOLD} 搜索窗={T_MAX_SEARCH_S / SEC_PER_HOUR:.0f} h")
        data = solve()
        ckpt.save(data, validate_constraints)
        print("[P3] 检查点已落盘（校验通过后写入）")
    else:
        data = cached
        print("[P3] 复用已验证检查点")

    cap = validate_capability(data)
    data["capability_verdict"] = cap
    td, diag = data["t_dry"], data["diagnostics"]
    print(f"[P3] t_dry={td['seconds']:.3f} s = {td['hours']:.4f} h"
          f"（区间 [{td['bracket_t_lo_s']:.0f},{td['bracket_t_hi_s']:.0f}] s，"
          f"max_r C {td['bracket_Cmax_lo']:.6f}→{td['bracket_Cmax_hi']:.6f}）")
    print(f"[P3] t_dry 处 max_r C={td['C_max_at_t_dry']:.6f}  C̄={td['C_bar_at_t_dry']:.6f}"
          f"  C(0)={td['C_center_at_t_dry']:.4f}  C(R0)={td['C_surface_at_t_dry']:.4f}")
    print(f"[P3] ε_mass={diag['eps_mass']:.3e}  Picard均值={diag['picard_mean']:.2f}"
          f"  峰值={diag['picard_max']}  步数={data['schedule']['n_steps']}"
          f"  墙钟={diag['wall_s']:.1f} s")
    print(f"[P3] Bi_m(表面)∈[{diag['Bi_m_range'][0]:.4e},{diag['Bi_m_range'][1]:.4e}]"
          f"  D_surface∈[{diag['D_surface_range'][0]:.4e},{diag['D_surface_range'][1]:.4e}]")
    for cid, v in cap.items():
        print(f"[P3][{cid}] {v['status']}")
    save_json(FIG_DIR / "problem_3_results.json", jsonable(data))
    print(f"[P3] 写出 figures/problem_3_results.json  "
          f"{(FIG_DIR / 'problem_3_results.json').stat().st_size / 1024:.1f} KB")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
