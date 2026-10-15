# -*- coding: utf-8 -*-
"""数据读取校验（CP1）+ 求解器冒烟测试。

行数基线来自 DATA_PROFILE.json，不硬编码；附件 2 半径单调不增（规则 R7、约束 K09）。
运行：python _utils/run_compute.py code/data_check.py
"""
from __future__ import annotations

import numpy as np

from params import *  # noqa: F401,F403
from utils import (Props, ensure_dirs, load_attachment1, load_attachment2, make_ambient,
                   make_grid, make_radius, run_case, save_json, _THOMAS_ERR)


def check_attachments() -> dict:
    t_air, temp_air, conc_air = load_attachment1()
    t_rad, r_cm = load_attachment2()
    print(f"[附件1] 行数={len(t_air)}（基线 {ATT1_ROWS}）  dt={np.unique(np.diff(t_air))} s"
          f"  t∈[{t_air[0]:.0f},{t_air[-1]:.0f}]")
    print(f"        T_a∈[{temp_air.min():.3f},{temp_air.max():.3f}] °C"
          f"  C_a∈[{conc_air.min():.5f},{conc_air.max():.5f}] kg/kg")
    print(f"[附件2] 行数={len(t_rad)}（基线 {ATT2_ROWS}）  dt={np.unique(np.diff(t_rad))} s"
          f"  R∈[{r_cm.min():.3f},{r_cm.max():.3f}] cm")

    tail = t_air >= PLATEAU_WINDOW_S
    ta_tail = float(temp_air[tail].mean())
    ca_tail = float(conc_air[tail].mean())
    assert int(tail.sum()) == PLATEAU_N, f"平台段点数 {int(tail.sum())} 偏离基线 {PLATEAU_N}"
    assert abs(ta_tail - PLATEAU_TA) < 5e-4, f"平台温度均值实测 {ta_tail:.4f} 偏离登记值"
    assert abs(ca_tail - PLATEAU_CA) < 5e-6, f"平台空气含水率均值实测 {ca_tail:.6f} 偏离登记值"
    print(f"[平台段] t>={PLATEAU_WINDOW_S:.0f}s  n={int(tail.sum())}"
          f"  T_a均值={ta_tail:.4f} °C  C_a均值={ca_tail:.6f} kg/kg（延拓用值）")

    dr = np.diff(r_cm)
    assert dr.max() <= 1e-12, "附件 2 半径出现上升"
    radius, radius_dot, rdot_max = make_radius()
    print(f"[PCHIP] dR/dt 最大值={rdot_max:.3e} m/s（须≤0）"
          f"  R(0)={radius(0.0) * CM_PER_M:.4f} cm  R(末)={radius(t_rad[-1]) * CM_PER_M:.4f} cm")

    t_amb, c_amb = make_ambient()
    print(f"[延拓] T_a(14400)={t_amb(ATT1_T_END_S):.4f}  T_a(200000)={t_amb(2e5):.4f}"
          f"  C_a(200000)={c_amb(2e5):.6f}")
    return {
        "att1_rows": int(len(t_air)), "att1_rows_baseline": ATT1_ROWS,
        "att2_rows": int(len(t_rad)), "att2_rows_baseline": ATT2_ROWS,
        "Ta_range": [float(temp_air.min()), float(temp_air.max())],
        "Ca_range": [float(conc_air.min()), float(conc_air.max())],
        "R_range_cm": [float(r_cm.min()), float(r_cm.max())],
        "plateau_n": int(tail.sum()), "plateau_Ta_mean": ta_tail, "plateau_Ca_mean": ca_tail,
        "R_monotone_max_diff_cm": float(dr.max()), "rdot_max_m_per_s": float(rdot_max),
        "R_ratio_end": float(r_cm[-1] / r_cm[0]),
    }


def check_properties() -> dict:
    """物性口径核对：R1（开尔文）、R2（D 随变干减小）、R3（湿基自变量）。"""
    out = {}
    for ap in (2, 3, 4):
        pr = Props(ap)
        c_hi, c_lo = np.array([C_INIT]), np.array([C_THRESHOLD])
        t_hot = np.array([PLATEAU_TA])
        d_hi = float(pr.D(c_hi, t_hot)[0])
        d_lo = float(pr.D(c_lo, t_hot)[0])
        assert d_hi > d_lo, f"附录 {ap} 的 D 未随变干减小（规则 R2）"
        out[f"appendix{ap}"] = {
            "rho_at_C0": float(pr.rho(c_hi)[0]), "cp_at_C0": float(pr.cp(c_hi)[0]),
            "lam_at_C0": float(pr.lam(c_hi)[0]),
            "D_at_C0_50C": d_hi, "D_at_Cth_50C": d_lo,
            "D_ratio": d_hi / d_lo,
            "Bi_h_at_C0": float(pr.bi_heat(c_hi)[0]),
            "Bi_m_at_C0": float(pr.bi_mass(c_hi, t_hot)[0]),
            "Bi_m_at_Cth": float(pr.bi_mass(c_lo, t_hot)[0]),
        }
        print(f"[附录{ap}] ρ={out[f'appendix{ap}']['rho_at_C0']:.2f}"
              f" c_p={out[f'appendix{ap}']['cp_at_C0']:.1f}"
              f" λ={out[f'appendix{ap}']['lam_at_C0']:.4f}"
              f" D(2.55,50℃)={d_hi:.4e} D(0.15,50℃)={d_lo:.4e}"
              f" Bi_h={out[f'appendix{ap}']['Bi_h_at_C0']:.4f}")
    pr3 = Props(3)
    c_surf = np.array([PLATEAU_CA])
    out["appendix3"]["D_at_Ca_50C"] = float(pr3.D(c_surf, np.array([PLATEAU_TA]))[0])
    out["appendix3"]["Bi_m_at_Ca"] = float(pr3.bi_mass(c_surf, np.array([PLATEAU_TA]))[0])
    print(f"[附录3] 表面平衡态 C={PLATEAU_CA:.5f} 时 D={out['appendix3']['D_at_Ca_50C']:.4e}"
          f"  Bi_m={out['appendix3']['Bi_m_at_Ca']:.4e}")
    return out


def smoke_solver() -> dict:
    """冒烟测试：附录 2、N=20、Δt=1 s 推进到 1800 s，核对与报告 10.2 预检量级一致。"""
    res = run_case({"appendix": 2, "N": 20, "out_dt": 100.0, "t_end": P1_WINDOW_S,
                    "dt1": 1.0, "record": "full", "tag": "smoke-P1-N20"})
    t_f, c_f = res["T_final"], res["C_final"]
    print(f"[冒烟] 1800 s: T(0)={t_f[0]:.4f} T(R0)={t_f[-1]:.4f}"
          f"  C(0)={c_f[0]:.4f} C(R0)={c_f[-1]:.4f}")
    print(f"[冒烟] ε_mass={res['eps_mass']:.3e}  Picard均值={res['picard_mean']:.2f}"
          f"  峰值={res['picard_max']}  步数={res['n_steps']}")
    assert 28.0 <= t_f.min() and t_f.max() < 50.8, "温度越界（约束 K02）"
    assert -1e-12 <= c_f.min() and c_f.max() <= C_INIT + 1e-9, "含水率越界（约束 K01）"
    assert res["eps_mass"] < 1e-2, "质量守恒误差超 1%（约束 K05）"
    assert np.all(np.diff(c_f) <= 1e-9), "含水率径向应单调不增（约束 K03）"
    return {"T0": float(t_f[0]), "TR": float(t_f[-1]), "C0": float(c_f[0]),
            "CR": float(c_f[-1]), "eps_mass": float(res["eps_mass"]),
            "picard_mean": float(res["picard_mean"]), "picard_max": int(res["picard_max"])}


def main() -> int:
    ensure_dirs()
    print(f"[口径] CALIBER_OK={CALIBER_OK}  Thomas/带状解偏差={_THOMAS_ERR:.3e}")
    grid = make_grid(N_GRID)
    print(f"[网格] N={N_GRID} Δs={grid[1]:.5f} ΣV={grid[2].sum():.15f}（须=0.5）"
          f"  输出步长对齐：stride={OUT_STRIDE}")
    payload = {
        "caliber_ok": bool(CALIBER_OK), "thomas_vs_banded_max_abs_diff": float(_THOMAS_ERR),
        "grid": {"N": N_GRID, "ds": float(grid[1]), "volume_sum": float(grid[2].sum()),
                 "out_stride": OUT_STRIDE},
        "attachments": check_attachments(),
        "properties": check_properties(),
        "smoke": smoke_solver(),
    }
    save_json(FIG_DIR / "data_check.json", payload)
    print("[data_check] 全部校验通过")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
