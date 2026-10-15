# 参赛提交材料说明

本压缩包为「参赛提交版」，只含提交必要文件：
- paper/    论文成品（main.pdf / Word）与 LaTeX 源码（main.tex、sections/、references.bib）
- code/     求解代码（*.py），依赖见 requirements.txt
- figures/  图表成品（PDF/PNG）与生成脚本
- user_data/ 代码所需的输入数据（赛题原文档不含在内）
- RESULTS.md 关键结果汇总
- AI工具使用详情.pdf 独立 AI 使用支撑材料（工作区已生成时随包保留）

## 如何运行代码
代码为多文件、互相 import，且大多用相对路径（如 code/xxx.py、figures/、_tmp/）。
请【在本 workspace 根目录】（即本 README 所在目录）运行，不要 cd 进 code/：
    pip install -r code/requirements.txt   # 依赖清单在 code/ 下（如无则忽略）
    python code/main.py                     # 入口脚本，按序跑完全部分析
提示：图表成品脚本（figures/gen_fig_*.py）可直接跑，用的是随包自带的结果数据；
      重跑 code/main.py 会重新生成中间结果并覆盖它们。

注意：
- 精简版不含运行时中间产物（_tmp/ 等），代码会在运行时自动重建；若某脚本报
  “找不到 _tmp / non-existent directory”，先在根目录手动建一个：mkdir _tmp
- 若要完整复现全部中间数据与产物，请用「导出全部」版本。

导出的是当前已存在的产物，不代表工作流已全部完成或已通过提交检查。
若包含 EXPORT_WARNINGS.txt，请先阅读其中的运行依赖/文件缺失提示；有提示时不保证代码可在软件外完整复现。
