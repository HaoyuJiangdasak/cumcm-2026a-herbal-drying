# -*- coding: utf-8 -*-
"""一次一因子灵敏度与鲁棒性（MODELING_REPORT 第八节表；能力 P3-C3、检查点 CP6）。

13 个扰动算例 + 基准，全部用问题 3 的同一求解器（附录 3、N=40、判据 max_r C < 0.15），
每次只改一个量，其余保持不变：

  D ±10%          扩散系数整体缩放（scale_D）
  h_m ±10%        表面传质系数缩放（scale_hm）
  h ±10%          对流换热系数缩放（scale_h）
  C_a ±10%        烘房空气含水率缩放（scale_Ca，含尾段平台值同步缩放）
  C_th=0.165/0.135 判据阈值 ±10%（同一条解轨迹上的不同停时）
  驱动力解释 B     表面通量以体积浓度作差，等价于空气侧参考浓度取 0（drive="B"）
  冻结温度         T ≡ T_a(t)，切断温度对 D 的耦合（freeze_T）
  含蒸发潜热       表面能量方程加 -L·h_m·ρ_d·(C_R-C_a) 汇项（latent_heat）

每个算例独立落检查点，可分批运行：
  python _utils/run_compute.py code/sensitivity.py -- --batch 1/3
全部检查点齐备后自动写出 figures/sensitivity_results.json。
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_utils"))
from compute_checkpoint import CaseCheckpoint  # noqa: E402

from params import *  # noqa: F401,F403,E402
from utils import (downsample_index, ensure_dirs, jsonable, run_case, save_json)  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
DEPS = ["code/utils.py", "code/params.py", "code/sensitivity.py",
        "user_data/附件1.xlsx", "user_data/附件2.xlsx"]
BASE_CFG = {"appendix": APPENDIX_OF_PROBLEM[3], "N": N_GRID, "out_dt": P3_OUT_DT_S,
            "t_end": T_MAX_SEARCH_S, "stop_threshold": C_THRESHOLD, "dt1": DT_STAGE1,
            "dt2": DT_STAGE2, "dt3": DT_STAGE3, "record": "scalar"}
UP = 1.0 + SENS_FRAC
DOWN = 1.0 - SENS_FRAC

# (id, 报告表中的行名, 相对基准改动的 cfg 覆盖项, 被扰动的量)
CASES = [
    ("base", "基准", {}, "—"),
    ("D_up", f"D +{SENS_FRAC:.0%}", {"scale_D": UP}, "D"),
    ("D_dn", f"D -{SENS_FRAC:.0%}", {"scale_D": DOWN}, "D"),
    ("hm_up", f"h_m +{SENS_FRAC:.0%}", {"scale_hm": UP}, "h_m"),
    ("hm_dn", f"h_m -{SENS_FRAC:.0%}", {"scale_hm": DOWN}, "h_m"),
    ("h_up", f"h +{SENS_FRAC:.0%}", {"scale_h": UP}, "h"),
    ("h_dn", f"h -{SENS_FRAC:.0%}", {"scale_h": DOWN}, "h"),
    ("Ca_up", f"C_a +{SENS_FRAC:.0%}", {"scale_Ca": UP}, "C_a"),
    ("Ca_dn", f"C_a -{SENS_FRAC:.0%}", {"scale_Ca": DOWN}, "C_a"),
    ("Cth_up", "C_th=0.165", {"stop_threshold": round(C_THRESHOLD * UP, 10)}, "C_th"),
    ("Cth_dn", "C_th=0.135", {"stop_threshold": round(C_THRESHOLD * DOWN, 10)}, "C_th"),
    ("drive_B", "驱动力解释 B", {"drive": "B"}, "驱动力参考基"),
    ("freeze_T", "冻结温度 T≡T_a(t)", {"freeze_T": True}, "温度耦合"),
    ("latent", f"含蒸发潜热 {LATENT_HEAT_J_PER_KG:.1e} J/kg", {"latent_heat": True}, "潜热汇项"),
]
CASE_IDS = [c[0] for c in CASES]
assert len(CASES) == 14 and len(set(CASE_IDS)) == 14, "灵敏度算例应为基准 + 13 扰动"


def case_cfg(case_id: str) -> dict:
    override = dict(next(c[2] for c in CASES if c[0] == case_id))
    cfg = dict(BASE_CFG)
    cfg.update(override)
    cfg["tag"] = f"SENS-{case_id}"
    return cfg


def make_validator(case_id: str):
    """从写出的结果复核 K01—K05、K11（不重跑求解器）。"""

    def validate(data: dict) -> bool:
        assert data.get("case_id") == case_id, "检查点算例号不符"
        c_th = float(data["stop_threshold"])
        assert data["t_dry_s"] is not None and data["t_dry_s"] > 0.0, "未在搜索窗内达标"
        lo_t, hi_t, lo_c, hi_c = data["bracket"]
        assert lo_c >= c_th > hi_c, "K11 跨阈值区间未真正跨过阈值"
        assert lo_t <= data["t_dry_s"] <= hi_t, "K11 t_dry 不在跨阈值区间内"
        assert abs((hi_t - lo_t) - P3_OUT_DT_S) < 1e-9, "K11 区间宽度不等于一个输出步"
        assert 0.0 - 1e-9 <= data["C_center_end"] <= C_INIT + 1e-9, "K01 中心含水率越界"
        assert 0.0 - 1e-9 <= data["C_surface_end"] <= C_INIT + 1e-9, "K01 表面含水率越界"
        assert data["C_center_end"] >= data["C_surface_end"] - 1e-12, "K03 中心不高于表面"
        c_max = np.asarray(data["series"]["C_max"], dtype=float)
        c_bar = np.asarray(data["series"]["C_bar"], dtype=float)
        assert float(np.max(np.diff(c_max))) <= 1e-12, "K04 max_r C 出现回升"
        assert float(np.max(np.diff(c_bar))) <= 1e-12, "K04 平均含水率出现回升"
        assert data["eps_mass"] < 1e-2, "K05 质量守恒误差超 1%"
        assert T_INIT - 1e-9 <= data["T_center_end"] <= data["Ta_max"] + 0.5, "K02 温度越界"
        return True

    return validate


def solve_case(case_id: str) -> dict:
    cfg = case_cfg(case_id)
    res = run_case(cfg)
    assert res["t_dry"] is not None, f"{case_id} 未在 {T_MAX_SEARCH_S / SEC_PER_HOUR:.0f} h 内达标"
    t_arr = np.asarray(res["t"], dtype=float)
    dsi = downsample_index(len(t_arr), 400)
    c_fin = np.asarray(res["C_final"], dtype=float)
    t_fin = np.asarray(res["T_final"], dtype=float)
    from utils import make_ambient
    t_amb, _ = make_ambient(scale_Ca=float(cfg.get("scale_Ca", 1.0)))
    return {
        "case_id": case_id,
        "label": next(c[1] for c in CASES if c[0] == case_id),
        "perturbed": next(c[3] for c in CASES if c[0] == case_id),
        "override": {k: v for k, v in next(c[2] for c in CASES if c[0] == case_id).items()},
        "stop_threshold": float(cfg["stop_threshold"]),
        "t_dry_s": float(res["t_dry"]), "t_dry_h": float(res["t_dry"] / SEC_PER_HOUR),
        "bracket": [float(v) for v in res["bracket"]],
        "C_center_end": float(c_fin[0]), "C_surface_end": float(c_fin[-1]),
        "C_bar_end": float(res["C_bar"][-1]), "C_max_end": float(res["C_max"][-1]),
        "T_center_end": float(t_fin[0]), "T_surface_end": float(t_fin[-1]),
        "Ta_max": float(max(t_amb(v) for v in t_arr)),
        "eps_mass": float(res["eps_mass"]), "n_steps": int(res["n_steps"]),
        "picard_max": int(res["picard_max"]), "picard_mean": float(res["picard_mean"]),
        "wall_s": float(res["wall_s"]),
        "D_center_range": [float(res["D_center"].min()), float(res["D_center"].max())],
        "D_surface_range": [float(res["D_surface"].min()), float(res["D_surface"].max())],
        "series": {"times_h": [float(t_arr[i] / SEC_PER_HOUR) for i in dsi],
                   "C_max": [float(res["C_max"][i]) for i in dsi],
                   "C_bar": [float(res["C_bar"][i]) for i in dsi]},
    }


def checkpoint(case_id: str) -> CaseCheckpoint:
    return CaseCheckpoint(ROOT, f"2026A/sensitivity/{case_id}", dependencies=DEPS,
                          parameters={"cfg": case_cfg(case_id), "picard_tol": PICARD_TOL,
                                      "gauss_n": GAUSS_N, "numpy": np.__version__})


def load_all() -> dict | None:
    """收集全部算例的已验证检查点；任一缺失返回 None。"""
    out = {}
    for cid in CASE_IDS:
        data = checkpoint(cid).load(make_validator(cid))
        if data is None:
            return None
        out[cid] = data
    return out


def assemble(cases: dict) -> dict:
    """组装报告第八节的灵敏度表：相对变化以基准 t_dry 为分母。"""
    base_h = cases["base"]["t_dry_h"]
    rows = []
    for cid in CASE_IDS:
        d = cases[cid]
        rel = (d["t_dry_h"] - base_h) / base_h * 100.0
        elast = None
        if cid.endswith(("_up", "_dn")) and cid != "drive_B":
            frac = SENS_FRAC if cid.endswith("_up") else -SENS_FRAC
            elast = (rel / 100.0) / frac
        rows.append({"case_id": cid, "label": d["label"], "perturbed": d["perturbed"],
                     "t_dry_h": d["t_dry_h"], "t_dry_h_r4": round(d["t_dry_h"], DECIMALS),
                     "t_dry_s": d["t_dry_s"],
                     "rel_change_pct": None if cid == "base" else rel,
                     "rel_change_pct_r4": None if cid == "base" else round(rel, DECIMALS),
                     "elasticity": None if elast is None else round(elast, DECIMALS),
                     "stop_threshold": d["stop_threshold"],
                     "C_bar_end": d["C_bar_end"], "eps_mass": d["eps_mass"],
                     "n_steps": d["n_steps"], "wall_s": d["wall_s"]})
    by_id = {r["case_id"]: r for r in rows}
    ranked = sorted((r for r in rows if r["case_id"] != "base"),
                    key=lambda r: abs(r["rel_change_pct"]), reverse=True)

    # 阈值敏感性的不对称：阈值降 10% 的增幅远大于升 10% 的降幅。
    asym = {"Cth_dn_pct": by_id["Cth_dn"]["rel_change_pct"],
            "Cth_up_pct": by_id["Cth_up"]["rel_change_pct"],
            "asymmetry_ratio": abs(by_id["Cth_dn"]["rel_change_pct"]
                                   / by_id["Cth_up"]["rel_change_pct"])}
    return {
        "problem": 3, "appendix": BASE_CFG["appendix"], "n_cases": len(rows),
        "sens_frac": SENS_FRAC, "baseline_t_dry_h": base_h,
        "baseline_t_dry_h_r4": round(base_h, DECIMALS),
        "baseline_case_id": "base",
        "criterion": {"expr": "t_dry = inf{t : max_r C(r,t) < 0.15}",
                      "threshold": C_THRESHOLD, "scope": CRITERION_SCOPE},
        "solver": {"N": N_GRID, "out_dt_s": P3_OUT_DT_S, "dt_stage_s": [DT_STAGE1, DT_STAGE2,
                                                                       DT_STAGE3],
                   "picard_tol": PICARD_TOL, "gauss_n": GAUSS_N,
                   "note": "全部算例共用问题 3 的求解器，一次只改一个量"},
        "table": rows,
        "ranking_by_abs_effect": [{"case_id": r["case_id"], "label": r["label"],
                                   "abs_rel_pct": abs(r["rel_change_pct"])} for r in ranked],
        "dominant_factor": ranked[0]["case_id"],
        "threshold_asymmetry": asym,
        "drive_B": {"t_dry_h": by_id["drive_B"]["t_dry_h"],
                    "rel_change_pct": by_id["drive_B"]["rel_change_pct"],
                    "interpretation": "表面通量以体积浓度 ρ_d C 作差，等价于空气侧参考浓度取 0",
                    "capability": "P3-C3"},
        "freeze_T": {"t_dry_h": by_id["freeze_T"]["t_dry_h"],
                     "rel_change_pct": by_id["freeze_T"]["rel_change_pct"],
                     "interpretation": "恒温段温度对传质的耦合可视为参数；预热段不可冻结"},
        "latent_heat": {"t_dry_h": by_id["latent"]["t_dry_h"],
                        "rel_change_pct": by_id["latent"]["rel_change_pct"],
                        "L_J_per_kg": LATENT_HEAT_J_PER_KG,
                        "interpretation": "假设(3) 的模型选择，影响量级最大，如实登记不宣称可忽略"},
        "robustness": {
            "all_within_S9_window": all(S9_DAYS_MIN * HOURS_PER_DAY <= r["t_dry_h"]
                                        <= S9_DAYS_MAX * HOURS_PER_DAY for r in rows),
            "S9_window_h": [S9_DAYS_MIN * HOURS_PER_DAY, S9_DAYS_MAX * HOURS_PER_DAY],
            "t_dry_h_range": [min(r["t_dry_h"] for r in rows),
                              max(r["t_dry_h"] for r in rows)],
            "max_eps_mass": max(r["eps_mass"] for r in rows),
            "note": "阈值扰动算例改的是判据本身，落在 S9 区间外不代表主模型失效"},
        "series_by_case": {cid: cases[cid]["series"] for cid in CASE_IDS},
        "detail_by_case": {cid: {k: v for k, v in cases[cid].items() if k != "series"}
                           for cid in CASE_IDS},
    }


def logic_probes(summary: dict) -> dict:
    """把合同 monotonic 的方向声明变成实算探针（logic_audit 的 audit_direction 用）。"""
    by_id = {r["case_id"]: r for r in summary["table"]}

    def sign(cid: str) -> int:
        v = by_id[cid]["rel_change_pct"]
        return 0 if v == 0 else (1 if v > 0 else -1)

    return {"monotonic": [
        {"more": "D_scale", "then": "t_dry", "observed_sign": sign("D_up"),
         "expect_sign": -1, "expect_dir": "better",
         "basis": "扩散越快越早达标（合同 monotonic 第 5 条）"},
        {"more": "C_threshold", "then": "t_dry", "observed_sign": sign("Cth_up"),
         "expect_sign": -1, "expect_dir": "better",
         "basis": "阈值越宽越早达标（合同 monotonic 第 6 条）"},
        {"more": "h_m_scale", "then": "t_dry", "observed_sign": sign("hm_up"),
         "expect_sign": -1, "expect_dir": "better",
         "basis": "表面传质越强越早达标"},
        {"more": "C_a_scale", "then": "t_dry", "observed_sign": sign("Ca_up"),
         "expect_sign": 1, "expect_dir": "worse",
         "basis": "空气越湿则驱动力越小、达标越晚"},
    ]}


def validate_capability(summary: dict) -> dict:
    """P3-C3：驱动力解释 B 必须是实算算例，并给出与主模型的相对偏差数值。"""
    b = summary["drive_B"]
    detail = summary["detail_by_case"]["drive_B"]
    assert detail["override"] == {"drive": "B"}, "P3-C3 解释 B 算例的覆盖项不是 drive=B"
    assert detail["t_dry_s"] > 0.0 and detail["n_steps"] > 0, "P3-C3 解释 B 没有实算步数"
    assert detail["eps_mass"] < 1e-2, "P3-C3 解释 B 算例质量守恒不达标"
    assert b["rel_change_pct"] is not None, "P3-C3 缺相对偏差数值"
    th = summary["threshold_asymmetry"]
    assert abs(b["rel_change_pct"]) < abs(th["Cth_dn_pct"]), \
        "P3-C3 需证明参考基模糊性的影响小于阈值本身的敏感度"
    return {"P3-C3": {
        "status": "PASS",
        "evidence": (f"灵敏度表含 drive=\"B\" 的实算算例（体积浓度作差，等价于空气侧参考浓度取 0）："
                     f"t_dry={b['t_dry_h']:.4f} h，相对主模型 {summary['baseline_t_dry_h']:.4f} h "
                     f"偏差 {b['rel_change_pct']:+.4f}%（实算步数 {detail['n_steps']}，"
                     f"ε_mass={detail['eps_mass']:.3e}，墙钟 {detail['wall_s']:.1f} s）；"
                     f"该偏差量级远小于阈值本身的敏感度（C_th-10% 为 "
                     f"{th['Cth_dn_pct']:+.4f}%），故参考基模糊性不影响结论的定性与工程量级。"),
        "t_dry_h": b["t_dry_h"], "rel_change_pct": b["rel_change_pct"]}}


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--batch", default=None, help="形如 2/3：只算第 2 批（共 3 批）")
    ap.add_argument("--case", default=None, help="只算指定算例号")
    args = ap.parse_args()

    ensure_dirs()
    todo = list(CASE_IDS)
    if args.case:
        assert args.case in CASE_IDS, f"未知算例 {args.case}"
        todo = [args.case]
    elif args.batch:
        k, m = (int(v) for v in args.batch.split("/"))
        assert 1 <= k <= m, "批次号非法"
        todo = [c for i, c in enumerate(CASE_IDS) if i % m == k - 1]

    for cid in todo:
        ck = checkpoint(cid)
        val = make_validator(cid)
        if ck.load(val) is not None:
            print(f"[SENS:{cid}] 复用已验证检查点")
            continue
        cfg = case_cfg(cid)
        print(f"[SENS:{cid}] 开始求解 {next(c[1] for c in CASES if c[0] == cid)}"
              f"  覆盖={cfg_override_str(cid)}")
        data = solve_case(cid)
        ck.save(data, val)
        print(f"[SENS:{cid}] t_dry={data['t_dry_h']:.4f} h  ε_mass={data['eps_mass']:.3e}"
              f"  步数={data['n_steps']}  墙钟={data['wall_s']:.1f} s")

    cases = load_all()
    if cases is None:
        missing = [c for c in CASE_IDS if checkpoint(c).load(make_validator(c)) is None]
        print(f"[SENS] 尚缺 {len(missing)} 个算例：{' '.join(missing)}；"
              f"补齐后自动写出 figures/sensitivity_results.json")
        return 0

    summary = assemble(cases)
    summary["logic_probes"] = logic_probes(summary)
    summary["capability_verdict"] = validate_capability(summary)
    print(f"\n[SENS] 基准 t_dry={summary['baseline_t_dry_h']:.4f} h"
          f"  共 {summary['n_cases']} 个算例（含基准）")
    print(f"{'算例':<26}{'t_dry (h)':>12}{'相对变化':>12}{'弹性':>8}")
    for r in summary["table"]:
        rel = "—" if r["rel_change_pct"] is None else f"{r['rel_change_pct']:+.2f}%"
        el = "" if r["elasticity"] is None else f"{r['elasticity']:+.2f}"
        print(f"{r['label']:<26}{r['t_dry_h']:>12.4f}{rel:>12}{el:>8}")
    th = summary["threshold_asymmetry"]
    print(f"[SENS] 主导因素={summary['dominant_factor']}  阈值不对称比="
          f"{th['asymmetry_ratio']:.4f}（-10%→{th['Cth_dn_pct']:+.2f}% vs "
          f"+10%→{th['Cth_up_pct']:+.2f}%）")
    print(f"[SENS] t_dry 全域 [{summary['robustness']['t_dry_h_range'][0]:.4f},"
          f"{summary['robustness']['t_dry_h_range'][1]:.4f}] h"
          f"  最大 ε_mass={summary['robustness']['max_eps_mass']:.3e}")
    for cid, v in summary["capability_verdict"].items():
        print(f"[SENS][{cid}] {v['status']}")
    save_json(FIG_DIR / "sensitivity_results.json", jsonable(summary))
    return 0


def cfg_override_str(case_id: str) -> str:
    ov = next(c[2] for c in CASES if c[0] == case_id)
    return "（基准）" if not ov else " ".join(f"{k}={v}" for k, v in ov.items())


if __name__ == "__main__":
    raise SystemExit(main())
