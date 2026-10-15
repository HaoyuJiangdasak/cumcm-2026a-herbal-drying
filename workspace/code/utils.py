# -*- coding: utf-8 -*-
"""四问共用的顶点式有限体积求解器（MODELING_REPORT 第五节）。

要点：控制体体积由 V_i=∫s ds 精确积分；界面扩散系数取式 (eq:iface) 的
Gauss-Legendre 积分平均（不是调和平均，见报告 1.1 节）；全隐式后向 Euler +
Picard 内迭代 + 三对角直接解；问题 4 用物质坐标 s=r/R(t)，1/R^2 进入算子、
1/R 进入 Robin 项，对流项被坐标变换吸收（报告 1.2 节），实验室坐标作为对照分支。
"""
from __future__ import annotations

import json
import math
import time

import numpy as np
import pandas as pd
from scipy.interpolate import PchipInterpolator
from scipy.linalg import solve_banded
from scipy.special import j0, j1

from params import *  # noqa: F401,F403

GAUSS_X, GAUSS_W = np.polynomial.legendre.leggauss(GAUSS_N)
C_FLOOR = 1e-6          # 物性求值前的含水率下截断（报告 5.3，远低于阈值 0.15）
DEGEN_TOL = 1e-12       # 界面积分平均退化为中点取值的判定


def make_grid(n_cells: int):
    """顶点式网格与无量纲控制体体积（式 eq:vol）。"""
    ds = 1.0 / n_cells
    s = np.arange(n_cells + 1, dtype=float) * ds
    V = np.empty(n_cells + 1)
    V[0] = ds * ds / 8.0
    V[1:n_cells] = np.arange(1, n_cells, dtype=float) * ds * ds
    V[n_cells] = (ds * ds / 2.0) * (n_cells - 0.25)
    face = (np.arange(n_cells, dtype=float) + 0.5) * ds
    assert abs(V.sum() - 0.5) < 1e-14, "离散域体积必须严格等于 ∫s ds=1/2"
    return s, ds, V, face


def thomas(a, b, c, d):
    """三对角直接解（参考实现，O(N) 前向消元 + 回代）。"""
    n = len(b)
    cp = np.zeros(n)
    dp = np.zeros(n)
    cp[0] = c[0] / b[0]
    dp[0] = d[0] / b[0]
    for i in range(1, n):
        m = b[i] - a[i] * cp[i - 1]
        if i < n - 1:
            cp[i] = c[i] / m
        dp[i] = (d[i] - a[i] * dp[i - 1]) / m
    x = np.empty(n)
    x[n - 1] = dp[n - 1]
    for i in range(n - 2, -1, -1):
        x[i] = dp[i] - cp[i] * x[i + 1]
    return x


def solve_tridiag(a, b, c, d, ab=None):
    """生产路径的三对角解：LAPACK 带状直接解，与 thomas() 数学等价（下方启动时断言）。"""
    n = len(b)
    if ab is None or ab.shape[1] != n:
        ab = np.zeros((3, n))
    ab[0, 0] = 0.0
    ab[0, 1:] = c[:-1]
    ab[1, :] = b
    ab[2, :-1] = a[1:]
    ab[2, -1] = 0.0
    return solve_banded((1, 1), ab, d, overwrite_ab=False, check_finite=False)


def assert_thomas_equivalence(seed: int = 7, n: int = 41, tol: float = 1e-12) -> float:
    """证明生产带状解与参考 Thomas 消元逐位一致（供 MC-THOMAS 与能力验收引用）。"""
    rng = np.random.default_rng(seed)
    off_l = -rng.uniform(0.1, 1.0, n)
    off_u = -rng.uniform(0.1, 1.0, n)
    off_l[0] = 0.0
    off_u[-1] = 0.0
    diag = np.abs(off_l) + np.abs(off_u) + rng.uniform(0.5, 1.5, n)
    rhs = rng.normal(size=n)
    x1 = thomas(off_l, diag, off_u, rhs)
    x2 = solve_tridiag(off_l, diag, off_u, rhs)
    err = float(np.max(np.abs(x1 - x2)))
    assert err < tol, f"带状解与 Thomas 消元不一致：{err:.3e}"
    return err


def d_iface(C, T, d_fun):
    """界面扩散系数 = 式 (eq:iface) 的 Kirchhoff 积分平均，8 点 Gauss-Legendre。

    D_{i+1/2} = 1/(C_{i+1}-C_i) ∫_{C_i}^{C_{i+1}} D(c,(T_i+T_{i+1})/2) dc
    """
    c_l = C[:-1]
    c_r = C[1:]
    t_mid = 0.5 * (T[:-1] + T[1:])
    d_c = c_r - c_l
    mid = 0.5 * (c_l + c_r)
    half = 0.5 * d_c
    c_q = mid[:, None] + half[:, None] * GAUSS_X[None, :]
    val = 0.5 * (d_fun(c_q, t_mid[:, None]) * GAUSS_W[None, :]).sum(axis=1)
    degen = np.abs(d_c) < DEGEN_TOL
    if degen.any():
        val[degen] = d_fun(c_l[degen], t_mid[degen])
    return val


def lam_iface(lam):
    """热传导系数界面值取调和平均（报告 5.2：λ 变化不足一倍，平均方式无实质差别）。"""
    left_v = lam[:-1]
    right_v = lam[1:]
    return 2.0 * left_v * right_v / (left_v + right_v)


def to_kelvin(t_degc):
    """Arrhenius 项必须用绝对温度（规则 R1）。摄氏度直接代入会使 D 相差约 1e54 倍。"""
    t_k = np.asarray(t_degc, dtype=float) + KELVIN
    assert np.all(t_k > 273.0), "D 求值入参不是绝对温度（K），检查 273.15 换算"
    return t_k


def wet_basis(C):
    """湿基含水率 C/(C+1)，c_p 与 λ 的自变量（规则 R3），不是干基 C。"""
    c = np.maximum(np.asarray(C, dtype=float), C_FLOOR)
    return c / (c + 1.0)


class Props:
    """按附录号装载物性；附录归属互斥（规则 R5），构造时断言不含跨附录常数。"""

    def __init__(self, appendix: int, *, scale_D: float = 1.0, frozen_state=None):
        assert appendix in (2, 3, 4), f"未知附录号 {appendix}"
        self.appendix = int(appendix)
        self.scale_D = float(scale_D)
        self.frozen_state = frozen_state   # (C,T) 时冻结全部物性 → 模型 M2
        self._check_isolation()

    def rho(self, C):
        if self.frozen_state is not None:
            C = np.full_like(np.asarray(C, dtype=float), self.frozen_state[0])
        c = np.maximum(np.asarray(C, dtype=float), C_FLOOR)
        if self.appendix == 2:
            return np.full_like(c, RHO_A2)
        if self.appendix == 3:
            return RHO3_A + RHO3_B * c
        return RHO4_A + RHO4_B * c

    def cp(self, C):
        if self.frozen_state is not None:
            C = np.full_like(np.asarray(C, dtype=float), self.frozen_state[0])
        w = wet_basis(C)
        if self.appendix == 2:
            return np.full_like(w, CP_A2)
        if self.appendix == 3:
            return CP3_A + CP3_B * w
        return CP4_A + CP4_B * w

    def lam(self, C):
        if self.frozen_state is not None:
            C = np.full_like(np.asarray(C, dtype=float), self.frozen_state[0])
        w = wet_basis(C)
        if self.appendix == 2:
            return np.full_like(w, LAM_A2)
        if self.appendix == 3:
            return LAM3_A + LAM3_B * w
        return LAM4_A + LAM4_B * w

    def D(self, C, T):
        if self.frozen_state is not None:
            C = np.full_like(np.asarray(C, dtype=float), self.frozen_state[0])
            T = np.full_like(np.asarray(C, dtype=float), self.frozen_state[1])
        c = np.maximum(np.asarray(C, dtype=float), C_FLOOR)
        if self.appendix == 2:
            return self.scale_D * D0_A2 * np.exp(-DC_A2 / c)
        t_k = to_kelvin(T)
        if self.appendix == 3:
            return self.scale_D * D0_A3 * np.exp(-DC_A3 / c) * np.exp(-DT_ARR_A3 / t_k)
        return self.scale_D * D0_A4 * np.exp(-DC_A4 / c) * np.exp(-DT_ARR_A4 / t_k)

    def _check_isolation(self):
        """运行时核对附录互斥（K07、K08）：本附录的常数必须命中，别的附录常数必须不命中。

        冻结物性（模型 M2、退化算例 V5）只改变入参取值、不改变附录归属，
        故 C→0 探针须临时绕过冻结，否则读到的是冻结态而非绝干极限。
        """
        frozen, self.frozen_state = self.frozen_state, None
        try:
            zero = np.zeros(1)
            dry_rho = float(self.rho(zero)[0])
            dry_lam = float(self.lam(zero)[0])
        finally:
            self.frozen_state = frozen
        want = {2: (RHO_A2, LAM_A2), 3: (RHO3_A, LAM3_A), 4: (RHO4_A, LAM4_A)}
        other = {2: (RHO3_A, RHO4_A), 3: (RHO_A2, RHO4_A), 4: (RHO_A2, RHO3_A)}
        rho_ref, lam_ref = want[self.appendix]
        # 容差取 1e-3：C_FLOOR=1e-6 使 C→0 极限有 O(1e-4) 残差，但三套附录常数相差 ≥90，可分辨
        assert abs(dry_rho - rho_ref) < 1e-3, f"附录 {self.appendix} 密度常数不符"
        assert abs(dry_lam - lam_ref) < 1e-5, f"附录 {self.appendix} 热传导常数不符"
        for wrong in other[self.appendix]:
            assert abs(dry_rho - wrong) > 1.0, f"附录 {self.appendix} 命中了别的附录常数"

    def bi_heat(self, C):
        return H_HEAT * R0 / self.lam(C)

    def bi_mass(self, C, T, radius=R0):
        return H_MASS * radius / self.D(C, T)


_AMBIENT_CACHE = {}
_RADIUS_CACHE = {}


def load_attachment1():
    """附件 1：烘房空气温度与水分浓度实测序列（药材侧 C 与空气侧 C_a 独立命名，K16）。"""
    if "att1" in _AMBIENT_CACHE:
        return _AMBIENT_CACHE["att1"]
    book = pd.read_excel(DATA_DIR / ATT1_FILE, sheet_name=None)
    assert SHEET_SINGLE in book, f"{ATT1_FILE} 缺 {SHEET_SINGLE} 工作表"
    df = book[SHEET_SINGLE]
    assert len(df) == ATT1_ROWS, f"{ATT1_FILE} 行数 {len(df)} 偏离 DATA_PROFILE 基线 {ATT1_ROWS}"
    assert list(df.columns) == ATT1_COLS, f"{ATT1_FILE} 列名偏离基线：{list(df.columns)}"
    t_air = df[ATT1_COLS[0]].to_numpy(dtype=float)
    temp_air = df[ATT1_COLS[1]].to_numpy(dtype=float)
    conc_air = df[ATT1_COLS[2]].to_numpy(dtype=float)
    assert np.all(np.diff(t_air) > 0.0), "附件 1 时间列必须严格递增"
    assert abs(t_air[-1] - ATT1_T_END_S) < 1e-9, "附件 1 时间跨度偏离基线"
    out = (t_air, temp_air, conc_air)
    _AMBIENT_CACHE["att1"] = out
    return out


def make_ambient(*, scale_Ca: float = 1.0, const_Ta=None, const_Ca=None):
    """建 T_a(t)、C_a(t) 线性插值器；越界按尾段实测平台均值延拓（假设 A3）。"""
    t_air, temp_air, conc_air = load_attachment1()
    hold_t = PLATEAU_TA
    hold_c = PLATEAU_CA * scale_Ca
    conc_scaled = conc_air * scale_Ca

    def t_amb(t):
        if const_Ta is not None:
            return float(const_Ta)
        return float(np.interp(t, t_air, temp_air, left=temp_air[0], right=hold_t))

    def c_amb(t):
        if const_Ca is not None:
            return float(const_Ca)
        return float(np.interp(t, t_air, conc_scaled, left=conc_scaled[0], right=hold_c))

    return t_amb, c_amb


def load_attachment2():
    """附件 2：药材半径实测序列（仅问题 4 使用，规则 R7 禁止自建收缩律）。"""
    if "att2" in _RADIUS_CACHE:
        return _RADIUS_CACHE["att2"]
    book = pd.read_excel(DATA_DIR / ATT2_FILE, sheet_name=None)
    assert SHEET_SINGLE in book, f"{ATT2_FILE} 缺 {SHEET_SINGLE} 工作表"
    df = book[SHEET_SINGLE]
    assert len(df) == ATT2_ROWS, f"{ATT2_FILE} 行数 {len(df)} 偏离 DATA_PROFILE 基线 {ATT2_ROWS}"
    assert list(df.columns) == ATT2_COLS, f"{ATT2_FILE} 列名偏离基线：{list(df.columns)}"
    t_rad = df[ATT2_COLS[0]].to_numpy(dtype=float)
    r_cm = df[ATT2_COLS[1]].to_numpy(dtype=float)
    assert np.all(np.diff(t_rad) > 0.0), "附件 2 时间列必须严格递增"
    assert np.all(np.diff(r_cm) <= 1e-12), "附件 2 半径必须单调不增（规则 R7、约束 K09）"
    assert abs(r_cm[0] - R0_CM) < 1e-9, "附件 2 首点半径应等于初始半径"
    out = (t_rad, r_cm)
    _RADIUS_CACHE["att2"] = out
    return out


def make_radius():
    """R(t) 的 PCHIP 保形插值（不用三次样条，避免拐点过冲造成 dR/dt>0）。"""
    t_rad, r_cm = load_attachment2()
    r_m = r_cm / CM_PER_M
    spline = PchipInterpolator(t_rad, r_m, extrapolate=True)
    deriv = spline.derivative()
    probe = np.linspace(t_rad[0], t_rad[-1], 20001)
    rdot_probe = deriv(probe)
    assert np.max(rdot_probe) <= 1e-14, f"PCHIP 出现 dR/dt>0：max={np.max(rdot_probe):.3e}"

    def radius(t):
        tc = min(max(float(t), t_rad[0]), t_rad[-1])
        return float(spline(tc))

    def radius_dot(t):
        tc = min(max(float(t), t_rad[0]), t_rad[-1])
        return float(deriv(tc))

    return radius, radius_dot, float(np.max(rdot_probe))


class Progress:
    """长求解的心跳输出：至少每 20 s 打印一次进度，避免无输出被误判为卡死。"""

    def __init__(self, tag: str, every: float = 20.0):
        self.tag = tag
        self.every = float(every)
        self.t0 = time.time()
        self.last = self.t0

    def tick(self, msg: str, force: bool = False):
        now = time.time()
        if force or now - self.last >= self.every:
            self.last = now
            print(f"  [{self.tag}] +{now - self.t0:7.1f}s  {msg}", flush=True)

    def done(self, msg: str = ""):
        print(f"  [{self.tag}] 完成，用时 {time.time() - self.t0:.1f}s  {msg}", flush=True)
        return time.time() - self.t0


def stage_dt(t: float, cfg: dict) -> float:
    """分段内部时间步（报告 10.4 的四位小数政策）。"""
    if t < T_SWITCH_1 - 1e-12:
        return float(cfg.get("dt1", DT_STAGE1))
    if t < T_SWITCH_2 - 1e-12:
        return float(cfg.get("dt2", DT_STAGE2))
    return float(cfg.get("dt3", DT_STAGE3))


def _const_radius(_t):
    return R0


def _zero_rate(_t):
    return 0.0


def run_case(cfg: dict) -> dict:
    """全隐式后向 Euler + Picard 内迭代推进耦合传热传质方程。

    cfg 关键字：appendix, N, out_dt, t_end, stop_threshold, extra_out_times,
    shrink, frame, dt1/dt2/dt3, scale_D/scale_hm/scale_h/scale_Ca, drive,
    freeze_T, latent_heat, frozen_state, const_Ta, const_Ca, record, tag。
    """
    n = int(cfg.get("N", N_GRID))
    s_node, ds, V, face = make_grid(n)
    props = Props(int(cfg["appendix"]), scale_D=float(cfg.get("scale_D", 1.0)),
                  frozen_state=cfg.get("frozen_state"))
    t_amb, c_amb = make_ambient(scale_Ca=float(cfg.get("scale_Ca", 1.0)),
                               const_Ta=cfg.get("const_Ta"), const_Ca=cfg.get("const_Ca"))
    h_heat = H_HEAT * float(cfg.get("scale_h", 1.0))
    h_mass = H_MASS * float(cfg.get("scale_hm", 1.0))
    drive = str(cfg.get("drive", DRIVE_REF))
    freeze_T = bool(cfg.get("freeze_T", False))
    latent = bool(cfg.get("latent_heat", LATENT_HEAT))
    frame = str(cfg.get("frame", "material"))
    record_full = str(cfg.get("record", "full")) == "full"
    assert drive in ("A", "B") and frame in ("material", "lab")

    if bool(cfg.get("shrink", False)):
        radius, radius_dot, rdot_max = make_radius()
    else:
        radius, radius_dot, rdot_max = _const_radius, _zero_rate, 0.0

    out_dt = float(cfg["out_dt"])
    c_th = cfg.get("stop_threshold")
    t_limit = float(cfg["t_end"]) if cfg.get("t_end") is not None else T_MAX_SEARCH_S
    targets = np.arange(1, int(round(t_limit / out_dt)) + 1, dtype=float) * out_dt
    if cfg.get("extra_out_times"):
        targets = np.union1d(targets, np.asarray(cfg["extra_out_times"], dtype=float))
        targets = targets[targets <= t_limit + 1e-9]

    T = np.full(n + 1, T_INIT)
    C = np.full(n + 1, C_INIT)
    if freeze_T:
        T[:] = t_amb(0.0)
    m0 = float(V @ C)
    flux_acc = 0.0
    ab_buf = np.zeros((3, n + 1))
    a_c = np.zeros(n + 1)
    b_c = np.zeros(n + 1)
    c_c = np.zeros(n + 1)
    d_c = np.zeros(n + 1)
    picard_hist = []

    def implicit_step(t_old, dt):
        """推进一步；返回 (T_new, C_new, 迭代次数, 表面通量项)。"""
        nonlocal T, C
        t_new = t_old + dt
        ta = t_amb(t_new)
        ca = c_amb(t_new)
        r_now = radius(t_new)
        inv_r2 = 1.0 / r_now**2
        rdot = radius_dot(t_new) if frame == "lab" else 0.0
        u_adv = s_node * rdot / r_now if frame == "lab" else None
        T_it = T.copy()
        C_it = C.copy()
        T_new = T_it
        C_new = C_it
        n_it = 0
        for _ in range(PICARD_MAX):
            n_it += 1
            if freeze_T:
                T_new = np.full(n + 1, ta)
            else:
                lam_f = lam_iface(props.lam(C_it)) * face * inv_r2 / ds
                cap = props.rho(C_it) * props.cp(C_it) * V / dt
                a_c[:] = 0.0
                c_c[:] = 0.0
                b_c[:] = cap
                d_c[:] = cap * T
                b_c[:n] += lam_f
                c_c[:n] = -lam_f
                b_c[1:] += lam_f
                a_c[1:] = -lam_f
                b_c[n] += h_heat / r_now
                d_c[n] += h_heat / r_now * ta
                if latent:
                    rho_dry = props.rho(C_it[n]) / (1.0 + max(C_it[n], C_FLOOR))
                    d_c[n] -= LATENT_HEAT_J_PER_KG * h_mass * rho_dry * (C_it[n] - ca) / r_now
                if u_adv is not None:
                    adv = props.rho(C_it) * props.cp(C_it) * V * u_adv / ds
                    b_c[1:] -= adv[1:]
                    a_c[1:] += adv[1:]
                T_new = solve_tridiag(a_c, b_c, c_c, d_c, ab_buf)

            d_f = d_iface(C_it, T_new, props.D) * face * inv_r2 / ds
            vcap = V / dt
            a_c[:] = 0.0
            c_c[:] = 0.0
            b_c[:] = vcap
            d_c[:] = vcap * C
            b_c[:n] += d_f
            c_c[:n] = -d_f
            b_c[1:] += d_f
            a_c[1:] = -d_f
            b_c[n] += h_mass / r_now
            if drive == "A":
                d_c[n] += h_mass / r_now * ca
            if u_adv is not None:
                adv_c = V * u_adv / ds
                b_c[1:] -= adv_c[1:]
                a_c[1:] += adv_c[1:]
            C_new = solve_tridiag(a_c, b_c, c_c, d_c, ab_buf)

            res = max(np.max(np.abs(C_new - C_it)), 1e-3 * np.max(np.abs(T_new - T_it)))
            C_it = C_new
            T_it = T_new
            if res < PICARD_TOL:
                break
        else:
            raise RuntimeError(f"Picard 未收敛（{PICARD_MAX} 次）于 t={t_new:.4f}s，残差 {res:.3e}")
        drive_gap = (C_new[n] - ca) if drive == "A" else C_new[n]
        return T_new, C_new, n_it, dt * h_mass / r_now * drive_gap

    prog = Progress(str(cfg.get("tag", "case")))
    t_cur = 0.0
    rec_t = [0.0]
    rec_T = [T.copy()]
    rec_C = [C.copy()]
    rec_flux = [0.0]
    rec_R = [radius(0.0)]
    rec_max = [float(C.max())]
    rec_bar = [2.0 * float(V @ C)]
    rec_dc = [float(props.D(np.array([C[0]]), np.array([T[0]]))[0])]
    rec_dsurf = [float(props.D(np.array([C[n]]), np.array([T[n]]))[0])]
    n_steps = 0
    it_max = 0
    it_sum = 0
    t_dry = None
    bracket = None
    for t_target in targets:
        dt_stage = stage_dt(t_cur, cfg)
        n_sub = max(1, int(math.ceil((t_target - t_cur) / dt_stage - 1e-9)))
        dt = (t_target - t_cur) / n_sub
        for _ in range(n_sub):
            T, C, n_it, d_flux = implicit_step(t_cur, dt)
            flux_acc += d_flux
            t_cur += dt
            n_steps += 1
            it_max = max(it_max, n_it)
            it_sum += n_it
        t_cur = float(t_target)
        c_max_now = float(C.max())
        rec_t.append(t_cur)
        rec_flux.append(flux_acc)
        rec_R.append(radius(t_cur))
        rec_max.append(c_max_now)
        rec_bar.append(2.0 * float(V @ C))
        rec_dc.append(float(props.D(np.array([C[0]]), np.array([T[0]]))[0]))
        rec_dsurf.append(float(props.D(np.array([C[n]]), np.array([T[n]]))[0]))
        if record_full:
            rec_T.append(T.copy())
            rec_C.append(C.copy())
        else:
            rec_T[0] = T.copy()
            rec_C[0] = C.copy()
        prog.tick(f"t={t_cur / SEC_PER_HOUR:8.3f} h  maxC={c_max_now:.5f}  steps={n_steps}")
        if c_th is not None and c_max_now < c_th:
            t_prev, c_prev = rec_t[-2], rec_max[-2]
            assert c_prev >= c_th > c_max_now, "跨阈值区间定位失败（约束 K11）"
            t_dry = t_prev + (c_prev - c_th) / (c_prev - c_max_now) * (t_cur - t_prev)
            bracket = (t_prev, t_cur, c_prev, c_max_now)
            break
    wall = prog.done(f"steps={n_steps} it_max={it_max}")

    m_end = float(V @ C)
    eps_mass = abs((m0 - m_end) - flux_acc) / m0
    out = {
        "cfg": {k: v for k, v in cfg.items() if not callable(v)},
        "N": n, "ds": ds, "s_nodes": s_node, "V": V,
        "t": np.asarray(rec_t), "C_max": np.asarray(rec_max), "C_bar": np.asarray(rec_bar),
        "flux_cum": np.asarray(rec_flux), "R": np.asarray(rec_R),
        "D_center": np.asarray(rec_dc), "D_surface": np.asarray(rec_dsurf),
        "T_final": T.copy(), "C_final": C.copy(),
        "t_end_reached": t_cur, "t_dry": t_dry, "bracket": bracket,
        "m0": m0, "m_end": m_end, "flux_total": flux_acc, "eps_mass": eps_mass,
        "n_steps": n_steps, "picard_max": it_max,
        "picard_mean": it_sum / max(1, n_steps), "wall_s": wall,
        "rdot_max": rdot_max, "record_full": record_full,
    }
    if record_full:
        out["T_field"] = np.asarray(rec_T)
        out["C_field"] = np.asarray(rec_C)
    return out


def bessel_roots(bi: float, n_terms: int) -> np.ndarray:
    """求 ζ J1(ζ)=Bi·J0(ζ) 的前 n_terms 个正根（无限长圆柱 Robin 特征值）。"""
    from scipy.optimize import brentq

    def char(z):
        return z * j1(z) - bi * j0(z)

    roots = []
    z = 1e-6
    step = 1e-3
    prev = char(z)
    while len(roots) < n_terms:
        z_next = z + step
        cur = char(z_next)
        if prev == 0.0:
            roots.append(z)
        elif prev * cur < 0.0:
            roots.append(brentq(char, z, z_next, xtol=1e-14, rtol=1e-15))
        z, prev = z_next, cur
    return np.asarray(roots)


def bessel_temperature(r_over_r0, t_s, *, bi: float, alpha: float, t_amb: float,
                       t_init: float, n_terms: int) -> np.ndarray:
    """无限长圆柱、常物性、恒定 T_a 与 Robin 边界的解析解（式 eq:bessel 前 n_terms 项）。"""
    zeta = bessel_roots(bi, n_terms)
    fo = alpha * t_s / R0**2
    coef = 2.0 * bi / ((zeta**2 + bi**2) * j0(zeta))
    x = np.atleast_1d(np.asarray(r_over_r0, dtype=float))
    series = (coef[None, :] * j0(zeta[None, :] * x[:, None])
              * np.exp(-zeta[None, :] ** 2 * fo)).sum(axis=1)
    return t_amb + (t_init - t_amb) * series


def ensure_dirs():
    for d in (FIG_DIR, OUT_DIR, TMP_DIR):
        d.mkdir(parents=True, exist_ok=True)


def node_columns(field, n_cells: int):
    """按节点抽取题面要求的 0.1 cm 输出列（偶数下标，无任何空间插值，MC-NODEOUT）。"""
    stride = n_cells // (N_OUT_RADII - 1)
    assert stride * (N_OUT_RADII - 1) == n_cells, "网格数必须是输出列数的整数倍"
    arr = np.asarray(field)
    return arr[..., ::stride]


def material_row_to_fixed_radii(row, s_nodes, r_now_cm):
    """问题 4：把物质坐标解映射到固定物理半径列；域外置空（约束 K14）。"""
    vals = []
    for r_cm in OUT_RADII_CM:
        if r_cm > r_now_cm + 1e-12:
            vals.append(None)
            continue
        s_query = min(r_cm / r_now_cm, 1.0)
        vals.append(float(np.interp(s_query, s_nodes, row)))
    return vals


def pick_time_index(t_arr, t_want: float) -> int:
    idx = int(np.argmin(np.abs(np.asarray(t_arr) - t_want)))
    assert abs(float(t_arr[idx]) - t_want) < 1e-6, f"输出时刻 {t_want} 未被记录"
    return idx


def downsample_index(n_total: int, max_rows: int) -> np.ndarray:
    if n_total <= max_rows:
        return np.arange(n_total)
    return np.unique(np.linspace(0, n_total - 1, max_rows).round().astype(int))


def round_out(x, nd: int = DECIMALS):
    if x is None:
        return None
    if isinstance(x, (list, tuple)):
        return [round_out(v, nd) for v in x]
    if isinstance(x, np.ndarray):
        return np.round(x, nd)
    return round(float(x), nd)


def jsonable(obj):
    if isinstance(obj, dict):
        return {str(k): jsonable(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [jsonable(v) for v in obj]
    if isinstance(obj, np.ndarray):
        return [jsonable(v) for v in obj.tolist()]
    if isinstance(obj, (np.floating, float)):
        val = float(obj)
        return None if (math.isnan(val) or math.isinf(val)) else val
    if isinstance(obj, (np.integer, int)):
        return int(obj)
    if isinstance(obj, (np.bool_, bool)):
        return bool(obj)
    return obj


def save_json(path, data) -> int:
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(jsonable(data), ensure_ascii=False, indent=1), encoding="utf-8")
    size = path.stat().st_size
    print(f"  [写出] {path.name}  {size / 1024:.1f} KB")
    return size


_THOMAS_ERR = assert_thomas_equivalence()
