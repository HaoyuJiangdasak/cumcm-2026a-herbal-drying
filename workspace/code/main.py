# -*- coding: utf-8 -*-
"""全流程编排入口：按依赖顺序跑齐八个计算脚本，再汇总 figures/all_results.json。

本文件自己不算任何物理量——它只做两件事：
  1. 调度：data_check → problem1 → problem2 → problem3 → problem4 → verification
     → sensitivity → model_compare → export_excel → constraint_audit。
     每个脚本内部都有 CaseCheckpoint（sha256 指纹依赖），已算过的算例会复用，
     所以重复运行 main.py 是幂等的、不会重算已成立的结果。
  2. 汇总：把各产物的顶层标量收拢到 figures/all_results.json，并把方向探针
     （logic_probes）合并成一份。all_results.json 只放标量与探针，不放场数据
     ——场数据留在各 problem_N_results.json 里，避免同一份数据两处存放。

用法：
  python _utils/run_compute.py code/main.py                    # 全流程（复用检查点）
  python _utils/run_compute.py code/main.py -- --only problem4 # 只跑一个阶段
  python _utils/run_compute.py code/main.py -- --summary-only  # 只重算汇总
"""
from __future__ import annotations

import argparse
import json
import subprocess
import sys
import time
from pathlib import Path

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "_utils"))

from params import *  # noqa: F401,F403,E402
from utils import ensure_dirs, jsonable, save_json  # noqa: E402

ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "_utils" / "run_compute.py"

# (阶段名, 脚本, 该阶段应当产出的文件)
STAGES = [
    ("data_check", "code/data_check.py", "figures/data_check.json"),
    ("problem1", "code/problem1.py", "figures/problem_1_results.json"),
    ("problem2", "code/problem2.py", "figures/problem_2_results.json"),
    ("problem3", "code/problem3.py", "figures/problem_3_results.json"),
    ("problem4", "code/problem4.py", "figures/problem_4_results.json"),
    ("verification", "code/verification.py", "figures/verification_results.json"),
    ("sensitivity", "code/sensitivity.py", "figures/sensitivity_results.json"),
    ("model_compare", "code/model_compare.py", "figures/model_compare_results.json"),
    ("export_excel", "code/export_excel.py", "result4.xlsx"),
    ("constraint_audit", "code/constraint_audit.py", "figures/constraint_audit_results.json"),
]
STAGE_NAMES = [s[0] for s in STAGES]
AUDIT_STAGE = "constraint_audit"      # 该阶段要复核 all_results.json，故汇总在它之前先写一遍


def run_stage(name: str, script: str, product: str) -> dict:
    """子进程跑一个阶段（隔离各脚本的 params 全局态），产物必须落盘才算成功。"""
    t0 = time.time()
    proc = subprocess.run([sys.executable, str(RUNNER), script],
                          cwd=str(ROOT), text=True, encoding="utf-8", errors="replace",
                          capture_output=True)
    wall = time.time() - t0
    tail = [ln for ln in (proc.stdout or "").splitlines() if ln.strip()][-3:]
    path = ROOT / product
    if not path.exists():
        why = f"产物 {product} 缺失"
    elif path.stat().st_size == 0:
        why = f"产物 {product} 为空"
    elif proc.returncode != 0:
        why = f"产物 {product} 已落盘（{path.stat().st_size / 1024.0:.1f} KB）但脚本自报失败"
    else:
        why = ""
    if why:
        sys.stdout.write(proc.stdout or "")
        sys.stderr.write(proc.stderr or "")
    assert not why, f"阶段 {name} 失败：returncode={proc.returncode}，{why}"
    print(f"[MAIN] {name:<17} {wall:7.1f} s  → {product}")
    for ln in tail:
        print(f"        {ln}")
    return {"stage": name, "script": script, "product": product,
            "wall_s": round(wall, 2), "size_kb": round(path.stat().st_size / 1024.0, 1),
            "returncode": proc.returncode, "log_tail": tail}


def _load(name: str) -> dict:
    path = FIG_DIR / name
    assert path.exists(), f"缺 {name}，先跑对应阶段"
    return json.loads(path.read_text(encoding="utf-8"))


def collect_scalars(with_audit: bool = True) -> dict:
    """顶层标量：键名带问题号，避免 logic_audit 的同名歧义判定。

    with_audit=False 时不含约束审计的计数（第一遍汇总在审计之前写，见 main 的两遍顺序）。
    """
    p1, p2, p3, p4 = (_load(f"problem_{i}_results.json") for i in (1, 2, 3, 4))
    ver, sens, cmp_ = (_load(f"{n}_results.json") for n in
                       ("verification", "sensitivity", "model_compare"))
    td3, td4 = p3["t_dry"], p4["main"]["t_dry"]
    eps = [p1["diagnostics"]["eps_mass"], p2["diagnostics"]["eps_mass"],
           p3["diagnostics"]["eps_mass"], p4["main"]["diagnostics"]["eps_mass"]]
    c_all = [np.asarray([[np.nan if v is None else v for v in row]
                         for row in _load(f"problem_{i}_results.json")["result_grid"]["C"]],
                        dtype=float) for i in (1, 2, 3, 4)]
    t_vals = [np.asarray(p["result_grid"]["T"], dtype=float).ravel() for p in (p1, p2)]
    for src in (p3["criterion_series"], p4["main"]["criterion_series"]):
        t_vals += [np.asarray(src["T_center"], dtype=float),
                   np.asarray(src["T_surface"], dtype=float)]
    gc = ver["grid_convergence"]
    worst_grid = max(abs(c["cases"][2][c["key_output"]] - c["cases"][1][c["key_output"]])
                     / abs(c["cases"][1][c["key_output"]]) * 100.0
                     for c in gc.values() if isinstance(c, dict) and "cases" in c)
    out = {
        "C_min_all": float(min(np.nanmin(c) for c in c_all)),
        "C_max_all": float(max(np.nanmax(c) for c in c_all)),
        "T_min_all": float(min(v.min() for v in t_vals)),
        "T_max_all": float(max(v.max() for v in t_vals)),
        "t_dry_Q3_s": float(td3["seconds"]), "t_dry_Q3_h": float(td3["hours"]),
        "t_dry_Q4_s": float(td4["seconds"]), "t_dry_Q4_h": float(td4["hours"]),
        "t_dry_Q3_days": float(td3["seconds"] / SEC_PER_HOUR / HOURS_PER_DAY),
        "C_bar_at_t_dry_Q3": float(td3["C_bar_at_t_dry"]),
        "C_bar_at_t_dry_Q4": float(td4["C_bar_at_t_dry"]),
        "R_ratio_Q4": float(td4["R_ratio"]),
        "mass_balance_error": float(max(eps)),
        "mass_balance_error_lab_frame_Q4": float(p4["lab"]["diagnostics"]["eps_mass"]),
        "grid_doubling_max_change_pct": float(worst_grid),
        "mass_conservation_eps_max_all_cases": float(ver["mass_conservation"]["eps_max"]),
        "t_dry_Q4_fixed_domain_h": float(p4["fixed"]["t_dry"]["hours"]),
        "frame_rel_diff_Q4": float(p4["frame_comparison"]["rel_diff"]),
        "sensitivity_baseline_t_dry_h": float(sens["baseline_t_dry_h"]),
        "sensitivity_n_cases": int(sens["n_cases"]),
        "t_dry_M1_lumped_h": float(cmp_["lumped_detail"]["t_dry_h"]),
        "t_dry_M2_frozen_h": float(cmp_["frozen_detail"]["t_dry_h"]),
    }
    if with_audit:
        ca = _load("constraint_audit_results.json")
        out["n_constraints_checked"] = int(ca["n_constraints"])
        out["n_constraints_failed"] = int(ca["n_failed"])
    return out


def _probe_d_vs_c(cs3: dict, sign) -> dict:
    """D 随 C 增（合同 monotonic 第 7 条）：只能在近等温段上测，否则被升温段污染。

    D = D0 exp(-a/C) exp(-b/T_K) 里 C 与 T 同时变。全程 0—28 min 是升温段，T 由 28 升到
    50 ℃ 抬高 D，同时 C 在降——两个反向效应叠在一起，全程相关系数会给出错的符号（实测
    -0.079）。故把探针限制在 T_center > 49.9 ℃ 的近等温平台段：那里 T 只在 0.1 ℃ 内变动，
    ΔD 的符号由 ΔC 支配。
    """
    c = np.asarray(cs3["C_max"], dtype=float)
    d = np.asarray(cs3["D_center"], dtype=float)
    t = np.asarray(cs3["T_center"], dtype=float)
    mask = t > 49.9
    cp, dp, tp = c[mask], d[mask], t[mask]
    dc, dd = np.diff(cp), np.diff(dp)
    n_agree = int(np.count_nonzero(np.sign(dc) == np.sign(dd)))
    return {"more": "C", "then": "D",
            "observed_sign": sign(float((dp[-1] - dp[0]) / (cp[-1] - cp[0]))),
            "expect_sign": 1, "expect_dir": "worse",
            "n_points": int(cp.size), "n_agreeing_pairs": n_agree,
            "n_pairs": int(dc.size), "T_span_degC": float(tp.max() - tp.min()),
            "basis": ("exp(-a/C) 的 C 在分母，故 D 随 C 增（第 7 条）。限定在 T>49.9 ℃ 的"
                      "近等温平台段测量，避免升温段的 T 效应反号；平台段内 ΔC 与 ΔD 同号"),
            "source": "problem_3_results.json 的 criterion_series 的 C_max/D_center/T_center"}


def build_probes() -> dict:
    """方向探针：合同 monotonic 七条逐条给实测符号，全部从已落盘产物现算。

    合同的 bounds 五项是物理界（极大值原理给出的 0≤C≤2.55 一类），不是「由删失值反解
    的上/下界」，没有「把输入朝真值方向推一点」的操作定义，故 bounds 留空不臆造探针。
    """
    p3, p4 = _load("problem_3_results.json"), _load("problem_4_results.json")
    sens = _load("sensitivity_results.json")
    cmp_ = _load("model_compare_results.json")
    cs3 = p3["criterion_series"]
    grid4 = p4["result_grid"]

    def sign(x: float) -> int:
        return 0 if x == 0.0 else (1 if x > 0 else -1)

    def trend(arr: np.ndarray) -> dict:
        """趋势符号取首末差（对含平段的序列稳健），非增性另记最大上升量为证据。"""
        return {"observed_sign": sign(float(arr[-1] - arr[0])),
                "max_rise": float(np.max(np.diff(arr))), "n_points": int(arr.size)}

    cmax = np.asarray(cs3["C_max"], dtype=float)
    cbar = np.asarray(cs3["C_bar"], dtype=float)
    crow = np.asarray(p3["result_grid"]["C"][-1], dtype=float)
    r_cm = np.asarray(grid4["R_cm"], dtype=float)
    probes = [
        {"more": "t", "then": "max_r C", **trend(cmax),
         "expect_sign": -1, "expect_dir": "better",
         "basis": "无源扩散 + 表面外向通量的极大值原理（合同 monotonic 第 1 条）",
         "source": "problem_3_results.json 的 criterion_series.C_max 全序列"},
        {"more": "t", "then": "volume_mean_C", **trend(cbar),
         "expect_sign": -1, "expect_dir": "better",
         "basis": "全局质量守恒 + 表面通量恒向外（第 2 条）",
         "source": "problem_3_results.json 的 criterion_series.C_bar"},
        {"more": "r", "then": "C", **trend(crow),
         "expect_sign": -1, "expect_dir": "better",
         "basis": "轴线零通量 + 表面排湿 + 均匀初值（第 3 条）",
         "source": "problem_3_results.json 的 result_grid.C 末行沿半径"},
        {"more": "t", "then": "R", **trend(r_cm),
         "expect_sign": -1, "expect_dir": "better",
         "basis": "附件 2 实测 + PCHIP 保形（第 4 条）；R 有平段，故趋势取首末差、"
                  "非增性由 max_rise=0 单独佐证",
         "source": "problem_4_results.json 的 result_grid.R_cm 全序列"},
        _probe_d_vs_c(cs3, sign),
    ]
    probes += list(sens["logic_probes"]["monotonic"])      # D/C_th/h_m/C_a 四条扰动实算
    # 收缩使达标更早：问题 4 的 shrink=1 主算例 vs shrink=0 固定域退化算例
    probes.append({"more": "shrinkage_enabled", "then": "t_dry",
                   "observed_sign": sign(p4["main"]["t_dry"]["hours"]
                                         - p4["fixed"]["t_dry"]["hours"]),
                   "expect_sign": -1, "expect_dir": "better",
                   "basis": "收缩缩短扩散路径，达标更早（式 eq:chain 的动边界项）",
                   "source": "problem_4_results.json 的 main 与 fixed 两分支 t_dry"})
    # 简化程度越高越早（错误地）判达标：M1 < M2 < M3
    by_model = {r["model"]: r["t_dry_h"] for r in cmp_["table"]}
    probes.append({"more": "model_fidelity", "then": "t_dry",
                   "observed_sign": sign(by_model["M3"] - by_model["M1"]),
                   "expect_sign": 1, "expect_dir": "worse",
                   "basis": "M1 集中参数忽略内部扩散阻力，会低估时长；保真度越高时长越长",
                   "source": "model_compare_results.json 的 table（M1 11.17 h → M3 57.53 h）"})
    return {"bounds": [], "monotonic": probes,
            "note": ("bounds 留空：合同的五项界是极大值原理给出的物理界，不是由删失/封顶值"
                     "反解的上下界，没有「把输入朝真值方向推一点再重算」的操作定义，"
                     "臆造探针会把物理界误判成方向错。monotonic 十一条全部实算。")}


def summarize(stage_log: list, with_audit: bool = True) -> dict:
    scal = collect_scalars(with_audit)
    probes = build_probes()
    p4 = _load("problem_4_results.json")
    payload = {
        "problem": "2026A", "title": "药材的烘干问题", "n_subproblems": 4,
        "generated_by": "code/main.py",
        **scal,
        "logic_probes": probes,
        "stages": stage_log,
        "per_problem": {
            "1": {"appendix": 2, "window_s": P1_WINDOW_S, "out_dt_s": P1_OUT_DT_S,
                  "deliverable": RESULT_FILES[1], "eps_mass": scal["mass_balance_error"],
                  "task": "预热阶段（常物性附录 2）逐秒温度场与水分场"},
            "2": {"appendix": 3, "window_s": P2_WINDOW_S, "out_dt_s": P2_OUT_DT_S,
                  "deliverable": RESULT_FILES[2],
                  "task": "变物性（附录 3）3 h 逐秒温度场与水分场"},
            "3": {"appendix": 3, "out_dt_s": P3_OUT_DT_S, "deliverable": RESULT_FILES[3],
                  "t_dry_h": scal["t_dry_Q3_h"], "t_dry_days": scal["t_dry_Q3_days"],
                  "task": "逐点判据 max_r C < 0.15 的最早达标时刻"},
            "4": {"appendix": 4, "out_dt_s": P4_OUT_DT_S, "deliverable": RESULT_FILES[4],
                  "t_dry_h": scal["t_dry_Q4_h"], "R_ratio": scal["R_ratio_Q4"],
                  "frame": p4["frame_main"],
                  "task": "动边界收缩域（物质坐标）下的达标时刻与移动表面含水率"},
        },
        "criterion": {"expr": "t_dry = inf{t : max_r C(r,t) < 0.15}",
                      "threshold": C_THRESHOLD, "scope": CRITERION_SCOPE},
        "provenance": ("标量与探针均从 figures/*_results.json 现算；场数据不在本文件重复存放。"
                       "各阶段由 code/main.py 子进程调度，检查点保证幂等。"),
    }
    return payload


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--only", action="append", default=[],
                    help=f"只跑指定阶段（可重复）：{' '.join(STAGE_NAMES)}")
    ap.add_argument("--summary-only", action="store_true", help="跳过计算，只重算汇总")
    args = ap.parse_args()

    ensure_dirs()
    for bad in [n for n in args.only if n not in STAGE_NAMES]:
        raise AssertionError(f"未知阶段 {bad}；可选：{' '.join(STAGE_NAMES)}")
    todo = [s for s in STAGES if not args.only or s[0] in args.only]

    log = []
    if not args.summary_only:
        print(f"[MAIN] 依赖顺序调度 {len(todo)} 个阶段（各脚本内部按检查点复用，重复运行幂等）")
        t0 = time.time()
        # 约束审计要复核 all_results.json，故汇总分两遍：审计前先按已有产物写一遍
        # （不含审计计数），审计跑完再写终版（补上计数与完整阶段日志）。
        for name, script, product in todo:
            if name == AUDIT_STAGE:
                save_json(FIG_DIR / "all_results.json",
                          jsonable(summarize(log, with_audit=False)))
                print("[MAIN] 汇总第一遍已落盘（不含审计计数），供约束审计复核")
            log.append(run_stage(name, script, product))
        print(f"[MAIN] 计算阶段合计 {time.time() - t0:.1f} s")
    else:
        print("[MAIN] --summary-only：跳过计算，仅从已落盘产物重算汇总")

    with_audit = (FIG_DIR / "constraint_audit_results.json").exists()
    payload = summarize(log, with_audit=with_audit)
    save_json(FIG_DIR / "all_results.json", jsonable(payload))
    print(f"[MAIN] 问题 3 t_dry={payload['t_dry_Q3_h']:.4f} h "
          f"（{payload['t_dry_Q3_days']:.2f} d，题面 S9 的 2—3 天内）"
          f"  C̄={payload['C_bar_at_t_dry_Q3']:.4f}")
    print(f"[MAIN] 问题 4 t_dry={payload['t_dry_Q4_h']:.4f} h  R/R0={payload['R_ratio_Q4']:.4f}"
          f"  C̄={payload['C_bar_at_t_dry_Q4']:.4f}")
    audit_txt = (f"；约束 {payload['n_constraints_checked']} 条违规 "
                 f"{payload['n_constraints_failed']} 条" if with_audit else "")
    print(f"[MAIN] 守恒最差 ε_mass={payload['mass_balance_error']:.3e}；"
          f"网格加倍最大变化 {payload['grid_doubling_max_change_pct']:.5f}%{audit_txt}")
    mono = payload["logic_probes"]["monotonic"]
    n_bad = sum(1 for m in mono if m.get("expect_sign") and m["observed_sign"]
                and (m["observed_sign"] > 0) != (m["expect_sign"] > 0))
    print(f"[MAIN] 方向探针 {len(mono)} 条，与合同声明相反 {n_bad} 条")
    assert n_bad == 0, "存在与合同声明方向相反的探针，先查模型符号"
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
