# 2026A 药材热风烘干建模

这是 2026 年全国大学生数学建模竞赛 A 题的参赛提交整理包，主题为药材热风烘干过程中的耦合传热传质、达标时间确定和收缩边界建模。

## 项目内容

- main_副本.pdf：项目根目录保留的原始提交 PDF。
- workspace/paper/：论文 LaTeX 工程、参考文献、论文 PDF 和排版改进记录。
- workspace/code/：四个问题的求解器、验证、灵敏度分析、模型对比、Excel 导出和约束审计代码。
- workspace/figures/：图表生成脚本、PDF 图表以及 JSON 计算结果。
- workspace/user_data/：附件 1、附件 2 和四个问题的结果表。
- workspace/RESULTS.md：详细结果、验证记录、灵敏度分析和已知偏差登记。
- workspace/_utils/：绘图辅助代码。

## 主要结论

统一模型使用圆柱坐标一维径向耦合传热传质方程、有限体积离散、全隐式后向 Euler 和 Picard 迭代。

- 问题 3（固定半径、变物性）：t_dry = 57.5316 h。
- 问题 4（变物性、动边界收缩）：t_dry = 51.1384 h，达标时半径比 R/R0 = 0.6000。
- 问题 3 和问题 4 的质量守恒残差分别约为 5.47e-13 和 3.32e-13。
- RESULTS.md 记录了网格收敛、时间步收敛、解析解交叉校验、灵敏度分析和 K01–K16 约束审计。

## 环境

代码依赖 Python 3.9 或更高版本，以及：

- NumPy
- pandas
- SciPy
- openpyxl
- Matplotlib

安装依赖：

~~~powershell
python -m pip install -r requirements.txt
~~~

## 复现说明

代码使用 workspace 作为工作目录，并大量依赖相对路径。输入数据和已生成的结果均已随项目保存。

当前归档缺少两个运行辅助文件：

- workspace/code/compute_checkpoint.py
- workspace/_utils/run_compute.py

因此，当前仓库适合阅读论文、检查代码、查看图表和复核已落盘结果；从干净环境直接运行完整流水线前，需要补回上述两个文件。已有结果不会被删除或覆盖。

论文源码入口为 workspace/paper/main.tex，完整结果说明见 workspace/RESULTS.md。

## 文件组织原则

项目原有目录和结果文件保持不变，只补充仓库级说明、依赖和忽略规则，便于后续协作和版本管理。原始说明保存在 README.original.txt。
