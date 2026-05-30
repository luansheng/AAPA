# AAPA 开发约束（Agent Instructions）

本文件约束所有 AI Agent 参与 AAPA 包开发时的行为规范。

## 1. 项目概述

AAPA（Anchor-Assisted Pedigree Assignment）是一个 R 包，用于基于已知候选父母本和 anchor 个体的高通量全同胞家系归属。核心算法结合 Mendelian 冲突检测与 anchor 亲缘打分，支持 top-k 剪枝和多规则拒判机制。

算法细节参见 [AAPA-research-plan.md](AAPA-research-plan.md)。

## 2. 代码风格

- 遵循 tidyverse 代码风格（参见 `.lintr` 配置）。
- 行宽上限 120 字符。
- 函数命名：`snake_case`。
- 公共 API 函数使用 `aapa_` 前缀（如 `aapa_assign`），或语义明确的动词短语（如 `read_genotype`、`qc_filter`）。
- 内部辅助函数使用 `.` 前缀（如 `.compute_allowed_genotypes`、`.is_compatible`）。
- 所有导出函数必须有完整的 Roxygen2 文档：`@param`、`@return`、`@family`、`@examples`。
- 不要添加不必要的注释；代码应自解释。

## 3. 包结构

```
aapa/
  R/
    data_io.R     # 数据读取：基因型、父母本、anchor
    scoring.R     # 打分引擎：Mendelian 冲突、anchor 亲缘、综合分数
    assign.R      # 主流程：归属管线与拒判逻辑
    qc.R          # 质量控制：缺失率、MAF 过滤
    plot.R        # 可视化：分数分布、top-k、拒判诊断
  src/            # C++ 代码（Rcpp 加速，当前为占位）
  tests/testthat/ # 单元测试
  vignettes/      # 使用指南
```

### 模块职责

- **data_io.R**：负责所有外部数据的读取与格式转换。大文件读取**必须**使用 `data.table::fread()`，禁止使用 `read.csv()` / `read.table()` 等 base R 函数读取大型数据文件。
- **scoring.R**：实现冲突计算、亲缘打分和综合分数。核心打分公式 `S = -α·C + β·K` 的修改需先在 research plan 中讨论。
- **assign.R**：主入口函数 `aapa_assign()`，编排整个流程。拒判规则（4 条）的变更需附带模拟实验验证。
- **qc.R**：SNP 和样本的质控过滤。
- **plot.R**：结果可视化，依赖 `ggplot2`（Suggests），必须提供 base R 回退方案。

新增模块需讨论后再添加，不得随意拆分或合并现有模块。

## 4. 数据处理

- **大文件读写一律使用 `data.table`**：`data.table::fread()` 读取、`data.table::fwrite()` 写出。
- `data.table` 是 Imports 依赖，可直接使用。
- 数据读取后，根据下游需要转换为 matrix 或 data.frame。内部计算矩阵（如基因型矩阵）继续使用 base R matrix。
- 处理表格类数据时，优先使用 `data.table` 的语法进行分组、聚合、筛选等操作，以保证大数据场景下的性能。

## 5. 依赖管理

### Imports（核心依赖，严格控制）

- `data.table`：大文件高性能读写与数据处理
- `cli`：用户友好的消息输出
- `checkmate`：参数校验
- `stats`、`utils`：基础统计与工具函数

**不得随意新增 Imports 依赖**。如需引入新依赖，必须说明理由并评估对包体积和可移植性的影响。

### Suggests（开发/测试依赖）

`testthat`、`ggplot2`、`knitr` 等仅用于开发、测试和文档。

## 6. 测试要求

- 所有导出函数必须有对应的单元测试。
- 使用 `testthat` edition 3。
- 测试文件命名：`test-{模块名}.R`（如 `test-scoring.R`）。
- 测试覆盖率目标 ≥ 80%。
- 新增功能必须同时提交测试。
- 修改现有函数时，确保现有测试通过，必要时补充回归测试。

## 7. CI/CD

已配置以下 GitHub Actions workflow：

- `R-CMD-check.yaml`：R CMD check
- `lint.yaml`：代码风格检查
- `test-coverage.yaml`：测试覆盖率
- `pkgdown.yaml`：文档站点构建

所有 PR 必须通过 CI 检查。不得绕过或禁用 CI。

## 8. 算法修改约束

- 核心打分公式的修改需先更新 `AAPA-research-plan.md`。
- 拒判规则的增删或参数默认值变更需附带模拟实验验证结果。
- 新增打分方法必须添加到 `scoring.R`，保持与现有接口风格一致。
- 性能优化不得改变数值结果（需回归测试验证数值等价）。

## 9. 性能与工程约束

- **MVP 阶段不引入 Rcpp**。先用纯 R 实现正确逻辑，profiling 后再迁移热点到 C++。
- 优先向量化操作，尽量减少显式 for 循环。
- 大矩阵操作注意内存占用，必要时分块处理。
- 目标运行规模：`N=10³~10⁵` 个体、`F=10~10²` 家系、`M=10³~10⁵` SNP。

## 10. Git 与 PR 规范

- 提交信息简洁明确，说明改了什么。
- 每个 PR 关联 issue 或在描述中说明动机。
- 不允许直接 push 到 main 分支。
- PR 需通过所有 CI 检查后合并。
