# -*- coding: utf-8 -*-
"""本文件由 code/gen_params.py 自动生成，请勿手改。

物理常数来源 PROBLEM_FACTS.json（题面逐字登记），实测统计量来源 DATA_FACTS.json，
离散参数来源 MODELING_REPORT.md 第 5.3 与第 10 节的收敛性实测结论。
文件末尾的口径一致性断言在每次 import 时复核数值与量纲，防止长上下文记忆漂移。
"""
from __future__ import annotations

import json
import math
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SOURCE_FACTS = ROOT / "PROBLEM_FACTS.json"
SOURCE_DATA_FACTS = ROOT / "DATA_FACTS.json"
SOURCE_PROFILE = ROOT / "DATA_PROFILE.json"
DATA_DIR = ROOT / "user_data"
FIG_DIR = ROOT / "figures"
OUT_DIR = ROOT / "output"
TMP_DIR = ROOT / "_tmp"

PROBLEM_ID = '2026A'
PROBLEM_TITLE = '药材的烘干问题'
N_SUBPROBLEMS = 4

# ---------------- 几何（题面：长 25 cm，半径 2 cm） ----------------
CM_PER_M = 100.0
LENGTH_CM = 25.0
R0_CM = 2.0
LENGTH_M = LENGTH_CM / CM_PER_M
R0 = R0_CM / CM_PER_M
GEOM_SHAPE = '圆柱形'
GEOM_MODE = "radial_1d_infinite_cylinder"

# ---------------- 初始条件与判据 ----------------
T_INIT = 28.0
C_INIT = 2.55
MOISTURE_BASIS = '干基含水率'
C_THRESHOLD = 0.15
CRITERION_SCOPE = '药材各处（逐点，全域最大值）'
TARGET_UNIT = 'h'
SEC_PER_HOUR = 3600.0

# ---------------- 传递系数（全题通用，附录 2 给定） ----------------
H_HEAT = 25.0
H_MASS = 8e-7
KELVIN = 273.15

# ---------------- 附录 2：仅问题 1（常物性） ----------------
RHO_A2 = 820.0
CP_A2 = 2600.0
LAM_A2 = 0.36
D0_A2 = 7e-9
DC_A2 = 0.89
D_EXPR_A2 = '7e-9 * exp(-0.89/C)'

# ---------------- 附录 3：仅问题 2、3（变物性 + Arrhenius） ----------------
RHO3_A = 650.0
RHO3_B = 128.0
CP3_A = 1450.0
CP3_B = 2736.0
LAM3_A = 0.21
LAM3_B = 0.38
D0_A3 = 2.4e-3
DC_A3 = 0.45
DT_ARR_A3 = 3850.0
D_EXPR_A3 = '2.4e-3 * exp(-0.45/C) * exp(-3850/T_K)'

# ---------------- 附录 4：仅问题 4（变物性 + 收缩） ----------------
RHO4_A = 760.0
RHO4_B = 90.0
CP4_A = 1850.0
CP4_B = 2150.0
LAM4_A = 0.12
LAM4_B = 0.2
D0_A4 = 4.2e-4
DC_A4 = 0.3
DT_ARR_A4 = 3850.0
D_EXPR_A4 = '4.2e-4 * exp(-0.30/C) * exp(-3850/T_K)'
APPENDIX_LABEL_OF_PROBLEM = {1: '附录2', 2: '附录3', 3: '附录3', 4: '附录4'}
APPENDIX_OF_PROBLEM = {k: int(re.sub(r"\D", "", v))
                       for k, v in APPENDIX_LABEL_OF_PROBLEM.items()}

# ---------------- 输出规格（题面「所有结果保留四位小数」） ----------------
DECIMALS = 4
OUT_DR_CM = 0.1
P1_WINDOW_S = 1800.0
P1_OUT_DT_S = 1.0
P1_TABLE_TIMES_S = [100.0, 300.0, 600.0, 900.0, 1200.0, 1500.0, 1800.0]
P1_TABLE_RADII_CM = [0.0, 0.5, 1.0, 1.5, 2.0]
HOURS_PER_DAY = 24.0
S9_DAYS_MIN = 2.0
S9_DAYS_MAX = 3.0
S9_SENTENCE = '烘干过程一般持续 2-3 天'
P2_TABLE_TIMES_H = [0.5, 1.0, 1.5, 2.0, 2.5, 3.0]
P2_TABLE_RADII_CM = [0.0, 0.5, 1.0, 1.5, 2.0]
P2_OUT_DT_S = 1.0
P2_WINDOW_S = max(P2_TABLE_TIMES_H) * SEC_PER_HOUR
P3_OUT_DT_S = 60.0
P3_TABLE_DT_H = 6.0
P3_TABLE_DR_CM = 0.5
P4_OUT_DT_S = 60.0
P4_TABLE_DT_H = 6.0
P4_TABLE_DR_CM = 0.5
P4_TABLE_LAST_COL = '药材表面（跟随移动边界）'
RESULT_FILES = {1: 'result1.xlsx', 2: 'result2.xlsx', 3: 'result3.xlsx', 4: 'result4.xlsx'}
SHEET_T = "温度"
SHEET_C = "水分浓度"
SHEET_SINGLE = "Sheet1"
COL_HEADER_CORNER = "时间\\到药材中心的距离"
N_OUT_RADII = int(round(R0_CM / OUT_DR_CM)) + 1
OUT_RADII_CM = [round(k * OUT_DR_CM, 10) for k in range(N_OUT_RADII)]

# ---------------- 附件实测元数据与平台段延拓值（DATA_FACTS 实测统计） ----------------
ATT1_FILE = "附件1.xlsx"
ATT2_FILE = "附件2.xlsx"
ATT1_ROWS = 241
ATT2_ROWS = 145
ATT1_COLS = ['时间', '温度', '水分浓度']
ATT2_COLS = ['时间', '半径']
ATT1_DT_S = 60.0
ATT2_DT_S = 1800.0
ATT1_T_END_S = (ATT1_ROWS - 1) * ATT1_DT_S
ATT2_T_END_S = (ATT2_ROWS - 1) * ATT2_DT_S
ATT1_TA_RANGE = [28.0, 50.246]
ATT1_CA_RANGE = [0.01963, 0.05025]
ATT2_R_RANGE_CM = [1.198, 2.0]
PLATEAU_WINDOW_S = 10800.0
PLATEAU_N = 61
PLATEAU_TA = 49.9989
PLATEAU_CA = 0.049988
PLATEAU_TA_STD = 0.1531
AMBIENT_EXTRAP = "plateau_hold"   # 假设 A3：t>ATT1_T_END_S 用尾段平台均值延拓

# ---------------- 离散参数（第 10 节收敛性实测确定，非题面给定） ----------------
N_GRID = 40
PICARD_TOL = 1e-11
PICARD_MAX = 60
GAUSS_N = 8
DT_STAGE1 = 0.01
DT_STAGE2 = 1.0
DT_STAGE3 = 60.0
DT_COARSE_STAGE1 = 1.0
T_SWITCH_1 = P1_WINDOW_S
T_SWITCH_2 = ATT1_T_END_S
GRID_CONV = [20, 40, 80, 160]
GRID_CONV_TDRY = [20, 40, 80]
DT_CONV = [120.0, 60.0, 30.0, 10.0]
P1_DT_CONV = [2.0, 1.0, 0.5, 0.25, 0.1]
BESSEL_DT = [1.0, 0.5, 0.1, 0.05]
BESSEL_N = 80
BESSEL_TERMS = 10
CONV_TOL_PCT = 0.1
BI_LUMPED_MAX = 0.1
T_MAX_SEARCH_S = 200.0 * SEC_PER_HOUR
SENS_FRAC = 0.1
DS = 1.0 / N_GRID
OUT_STRIDE = int(round(OUT_DR_CM / R0_CM * N_GRID))

# ---------------- 建模开关（MODELING_REPORT 第二节假设的参数化） ----------------
TRANSPORT_MODEL = "fickian_concentration_driven"
LATENT_HEAT = 0            # 假设(3)：主模型不含蒸发潜热汇项，敏感性中开启对照
LATENT_HEAT_J_PER_KG = 2.40e6
DRIVE_REF = "A"            # 解释 A：表面干基含水率与空气含水率直接相减；B 为对照
SHRINK_MODE = "given_series"   # 假设(6)：R(t) 由附件 2 给定，禁止自建收缩律
FRAME = "material"         # 问题 4 主模型用物质坐标；lab 为对照分支

# ---------------- 派生无量纲数（用于排除集中参数与表面平衡模型） ----------------
BI_H_A2 = H_HEAT * R0 / LAM_A2
D_INIT_A2 = D0_A2 * math.exp(-DC_A2 / C_INIT)
BI_M_A2 = H_MASS * R0 / D_INIT_A2

# ---------------- 口径一致性断言（每次 import 复核，防长上下文漂移） ----------------
_F = json.loads(SOURCE_FACTS.read_text(encoding="utf-8"))
_S = _F["sub_problems"]


def _eq(a, b, tol=1e-12):
    return abs(float(a) - float(b)) <= tol * max(1.0, abs(float(b)))


assert _eq(R0_CM, _F["geometry"]["radius_cm"]), "R0_CM 偏离 PROBLEM_FACTS"
assert _eq(LENGTH_CM, _F["geometry"]["length_cm"]), "LENGTH_CM 偏离 PROBLEM_FACTS"
assert _eq(T_INIT, _F["initial_conditions"]["T_init_degC"]), "T_INIT 偏离 PROBLEM_FACTS"
assert _eq(C_INIT, _F["initial_conditions"]["C_init_kg_per_kg"]), "C_INIT 偏离 PROBLEM_FACTS"
assert _eq(C_THRESHOLD, _F["criterion"]["C_threshold_kg_per_kg"]), "C_THRESHOLD 偏离 PROBLEM_FACTS"
assert _eq(H_HEAT, _F["transfer_coefficients"]["h_convective_heat_W_m2K"]), "H_HEAT 偏离 PROBLEM_FACTS"
assert _eq(H_MASS, float(_F["transfer_coefficients"]["h_mass_expr"])), "H_MASS 偏离 PROBLEM_FACTS"
assert _eq(RHO_A2, _S["problem1"]["rho_const_kg_m3"]) and _eq(CP_A2, _S["problem1"]["cp_const_J_kgK"])
assert _eq(LAM_A2, _S["problem1"]["lambda_const_W_mK"]), "LAM_A2 偏离附录 2"
assert _eq(D0_A2, float(_S["problem1"]["D_prefactor_expr"])) and _eq(DC_A2, _S["problem1"]["D_conc_coeff"])
assert _eq(RHO3_A, _S["problem2"]["rho_intercept"]) and _eq(RHO3_B, _S["problem2"]["rho_slope"])
assert _eq(CP3_A, _S["problem2"]["cp_intercept"]) and _eq(CP3_B, _S["problem2"]["cp_slope"])
assert _eq(LAM3_A, _S["problem2"]["lambda_intercept"]) and _eq(LAM3_B, _S["problem2"]["lambda_slope"])
assert _eq(D0_A3, float(_S["problem2"]["D_prefactor_expr"])) and _eq(DC_A3, _S["problem2"]["D_conc_coeff"])
assert _eq(DT_ARR_A3, _S["problem2"]["D_arrhenius_coeff"]), "附录 3 Arrhenius 系数偏离"
assert _eq(RHO4_A, _S["problem4"]["rho_intercept"]) and _eq(RHO4_B, _S["problem4"]["rho_slope"])
assert _eq(CP4_A, _S["problem4"]["cp_intercept"]) and _eq(CP4_B, _S["problem4"]["cp_slope"])
assert _eq(LAM4_A, _S["problem4"]["lambda_intercept"]) and _eq(LAM4_B, _S["problem4"]["lambda_slope"])
assert _eq(D0_A4, float(_S["problem4"]["D_prefactor_expr"])) and _eq(DC_A4, _S["problem4"]["D_conc_coeff"])
assert _eq(DT_ARR_A4, _S["problem4"]["D_arrhenius_coeff"]), "附录 4 Arrhenius 系数偏离"
assert DECIMALS == int(_F["output_precision"]["decimals"]), "四位小数政策偏离题面"
assert _eq(S9_DAYS_MIN, _S["problem2"]["duration_days_min"]) and \
       _eq(S9_DAYS_MAX, _S["problem2"]["duration_days_max"]), "S9 时长区间偏离题面"
assert 0.0 < S9_DAYS_MIN < S9_DAYS_MAX, "S9 区间必须为正且有序"
assert HOURS_PER_DAY * S9_DAYS_MAX * SEC_PER_HOUR <= T_MAX_SEARCH_S, \
    "搜索窗须覆盖 S9 的 2—3 天量级"

# 量纲与口径
assert _eq(R0, 0.02) and _eq(LENGTH_M, 0.25), "cm→m 换算错误"
assert KELVIN == 273.15, "Arrhenius 项必须用绝对温度，偏移量固定 273.15"
assert T_INIT + KELVIN > 0 and PLATEAU_TA + KELVIN > 0, "开尔文温度必须为正"
assert APPENDIX_OF_PROBLEM[2] == APPENDIX_OF_PROBLEM[3] == 3, "问题 2、3 共用附录 3"
assert APPENDIX_OF_PROBLEM[1] == 2 and APPENDIX_OF_PROBLEM[4] == 4, "附录归属必须互斥"
assert set(APPENDIX_OF_PROBLEM.values()) == {2, 3, 4}, "附录编号只能取 2/3/4"
assert RHO_A2 != RHO3_A != RHO4_A, "三套附录密度常数不得混用"

# 物性单调性与取值域（规则 R2、R3）
assert D0_A2 > 0.0 and D0_A3 > 0.0 and D0_A4 > 0.0, "扩散前因子必须为正"
assert DC_A2 > 0.0 and DC_A3 > 0.0 and DC_A4 > 0.0, "浓度项系数为正，exp(-a/C) 随 C 减小而减小"
assert math.exp(-DC_A2 / C_INIT) > math.exp(-DC_A2 / C_THRESHOLD), "D 必须随含水率下降而减小"
assert math.exp(-DC_A3 / C_INIT) > math.exp(-DC_A3 / C_THRESHOLD), "附录 3 的 D 单调性错误"
assert math.exp(-DC_A4 / C_INIT) > math.exp(-DC_A4 / C_THRESHOLD), "附录 4 的 D 单调性错误"
assert 0.0 < C_INIT / (C_INIT + 1.0) < 1.0, "湿基含水率必须落在 [0,1)"
assert RHO3_A + RHO3_B * C_INIT > 0.0 and RHO4_A + RHO4_B * C_INIT > 0.0, "密度为正"
assert CP3_A + CP3_B * C_INIT / (C_INIT + 1.0) > 0.0, "比热容为正"
assert LAM3_A + LAM3_B * C_INIT / (C_INIT + 1.0) > 0.0, "热传导系数为正"
assert LAM4_A + LAM4_B * C_INIT / (C_INIT + 1.0) > 0.0, "热传导系数为正"
assert RHO4_A > RHO3_A > 0.0, "附录 4 绝干密度高于附录 3，与收缩致密化一致"

# 无量纲数：排除集中参数模型与表面平衡假设
assert 0.1 < BI_H_A2 < 100.0, "热 Biot 数不在分布参数区间，集中参数模型不可用"
assert BI_M_A2 > 1.0, "质量 Biot 数大于 1，表面不可视为瞬时平衡"

# 输出网格必须与节点对齐（禁止空间插值，MC-NODEOUT）
assert N_GRID % OUT_STRIDE == 0, "N_GRID 必须能被输出步长整除"
assert _eq(OUT_STRIDE * DS * R0_CM, OUT_DR_CM), "输出步长与网格不对齐"
assert N_OUT_RADII == 21 and _eq(OUT_RADII_CM[-1], R0_CM), "输出位置应为 0..2 cm 共 21 列"
assert len(P1_TABLE_TIMES_S) == 7 and len(P2_TABLE_TIMES_H) == 6, "论文表格时刻数偏离题面"

# 附件元数据与延拓值
assert _eq(ATT1_T_END_S, 14400.0) and _eq(ATT2_T_END_S, 259200.0), "附件时间跨度偏离实测"
assert ATT1_TA_RANGE[0] <= PLATEAU_TA <= ATT1_TA_RANGE[1], "平台温度超出实测范围"
assert ATT1_CA_RANGE[0] <= PLATEAU_CA <= ATT1_CA_RANGE[1], "平台空气含水率超出实测范围"
assert PLATEAU_CA < C_THRESHOLD < C_INIT, "空气含水率与药材含水率量级混淆（规则 R6）"
assert ATT2_R_RANGE_CM[1] <= R0_CM, "附件 2 半径不得超过初始半径"
assert T_MAX_SEARCH_S > ATT2_T_END_S * 0.5, "搜索窗须覆盖预期烘干时长"
assert 0.0 < DT_STAGE1 < DT_STAGE2 < DT_STAGE3, "分段时间步必须递增"
assert PICARD_TOL < 1e-9 and PICARD_MAX >= 20, "Picard 容差或迭代上限过松"

CALIBER_OK = True

PARAM_SUMMARY = {
    "R0_m": R0, "T_init_degC": T_INIT, "C_init": C_INIT, "C_threshold": C_THRESHOLD,
    "h": H_HEAT, "h_m": H_MASS, "N_grid": N_GRID, "Bi_h_appendix2": BI_H_A2,
    "Bi_m_appendix2": BI_M_A2, "plateau_Ta": PLATEAU_TA, "plateau_Ca": PLATEAU_CA,
}
