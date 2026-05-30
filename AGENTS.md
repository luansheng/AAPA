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

新增模块需讨论后再添加。原则上不随意拆分或合并现有模块；若为降低复杂度或隔离热点而重构，需在 PR 中说明动机与边界。

## 4. 对象契约与不变量

### genotype matrix

- `genotype` 必须是 base R matrix；行为个体、列为 SNP/marker。
- 行名必须是唯一的个体 ID；列名必须是唯一的 marker ID。
- 基因型编码语义固定为 `0/1/2/NA`，且同一数据集内参考等位基因方向必须一致。
- 下游函数不得静默重排样本或 marker 顺序；如需对齐，必须在单一入口显式完成并校验。

### aapa_parents

- `aapa_parents` 的每个元素必须包含 `family_id`、`sire_id`、`dam_id`、`sire_geno`、`dam_geno`。
- `names(parents)` 必须与 `family_id` 一致，且在整个流程中保持稳定。
- 父母本 ID 必须能回溯到 `genotype` 的行名；父母本基因型长度必须与 marker 数一致。

### aapa_anchors

- `aapa_anchors` 必须至少包含 `individual_id`、`family_id`；`weight` 如存在，其语义必须固定且缺省值一致。
- `attr(anchors, "geno")` 必须与 `individual_id` 顺序严格对齐。
- anchor 的 `family_id` 必须存在于候选家系集合中；缺失或不一致应尽早报错，不允许静默跳过。

### aapa_result

- `aapa_result` 中的 `assignment`、`topk`、`conflict_matrix`、`kinship_matrix`、`score_matrix`、`params` 必须相互对应。
- `assignment$individual_id` 必须与结果矩阵行名可对齐；家系列名必须在各矩阵中保持一致。
- `assigned_family`、`status`、`reject_reason` 的语义必须稳定，不因内部重构而改变解释口径。

## 5. 分层规则

- **用户层**：导出函数、S3 方法、文档示例；负责稳定 API 和用户可理解输出。
- **workflow 层**：主流程编排、top-k 剪枝、拒判规则、结果装配。
- **computation 层**：Mendelian conflict、anchor kinship、综合打分等纯计算逻辑。
- **infrastructure 层**：文件 I/O、对象构造辅助、编译配置、C++ bridge。

约束：

- computation 层不得负责消息输出、S3 类恢复、报表格式化。
- workflow 层不得重复实现底层数值逻辑。
- infrastructure 层不得承载算法决策或拒判规则。
- C++ 代码只承担计算核，不承担输入校验、文档示例、S3 方法和用户消息输出。

## 6. 数据处理

- **大文件读写一律使用 `data.table`**：`data.table::fread()` 读取、`data.table::fwrite()` 写出。
- `data.table` 是 Imports 依赖，可直接使用。
- 数据读取后，根据下游需要转换为 matrix 或 data.frame。内部计算矩阵（如基因型矩阵）继续使用 base R matrix。
- 处理表格类数据时，优先使用 `data.table` 的语法进行分组、聚合、筛选等操作，以保证大数据场景下的性能。
- 对 `connection` 等特殊输入可保留受控例外，但不得把该例外扩展为大文件默认路径。

## 7. 依赖管理

### Imports（核心依赖，严格控制）

- `data.table`：大文件高性能读写与数据处理
- `cli`：用户友好的消息输出
- `checkmate`：参数校验
- `stats`、`utils`：基础统计与工具函数

**不得随意新增 Imports 依赖**。如需引入新依赖，必须说明理由并评估对包体积和可移植性的影响。

### Suggests（开发/测试依赖）

`testthat`、`ggplot2`、`knitr` 等仅用于开发、测试和文档。

## 8. API 与反膨胀规则

- 非必要不新增导出函数；新增公共 API 必须说明长期维护价值。
- 新 helper 引入前，先检查能否复用或扩展现有 helper。
- 不接受仅包一层旧逻辑的新 wrapper；若出现重复路径，优先合并并删除旧实现。
- 新增参数必须有稳定语义，不为一次性实验或临时兼容暴露接口。
- 修改内部结构时，应优先减少系统复杂度，而不是叠加兼容层。

## 9. 测试要求

- 所有导出函数必须有对应的单元测试。
- 使用 `testthat` edition 3。
- 测试文件命名：`test-{模块名}.R`（如 `test-scoring.R`）。
- 测试覆盖率目标 ≥ 80%。
- 新增功能必须同时提交测试。
- 修改现有函数时，确保现有测试通过，必要时补充回归测试。
- 测试至少覆盖四类场景：
  - **invariant tests**：对象结构、维度、命名、ID 对齐。
  - **regression tests**：历史 bug、边界条件、拒判规则。
  - **failure-path tests**：非法输入、缺失样本、未对齐对象、空结果。
  - **round-trip tests**：读取/构造对象后进入主流程，输出保持一致。
- 引入 C++ 或并发优化时，必须保留纯 R 对照测试：整数结果严格一致，浮点结果在约定容差内一致。

## 10. CI/CD

已配置以下 GitHub Actions workflow：

- `R-CMD-check.yaml`：R CMD check
- `lint.yaml`：代码风格检查
- `test-coverage.yaml`：测试覆盖率
- `pkgdown.yaml`：文档站点构建

所有 PR 必须通过 CI 检查。不得绕过或禁用 CI。

## 11. 算法修改约束

- 核心打分公式的修改需先更新 `AAPA-research-plan.md`。
- 拒判规则的增删或参数默认值变更需附带模拟实验验证结果。
- 新增打分方法必须添加到 `scoring.R`，保持与现有接口风格一致。
- 性能优化不得改变数值结果、`NA` 语义、排序语义与拒判语义（需回归测试验证数值等价）。
- 算法规则变更若影响热点实现，需同时说明对纯 R 参考实现与未来 C++ 路径的一致性要求。

## 12. 性能与工程约束

- **MVP 阶段以纯 R 为准**。先用纯 R 实现正确逻辑，profiling 后再迁移热点到 C++。
- 优先向量化操作，尽量减少显式 for 循环。
- 大矩阵操作注意内存占用，必要时分块处理。
- 目标运行规模：`N=10³~10⁵` 个体、`F=10~10²` 家系、`M=10³~10⁵` SNP。
- 性能优化必须遵循 **profile-first** 原则；未经 profiling 不得先行引入复杂优化。
- 纯 R 实现是 truth oracle；任何编译优化都必须以纯 R 路径作为数值对照。
- 仅将已确认热点下沉到 C++，例如 Mendelian conflict、anchor kinship、top-k 等计算核；I/O、参数校验、对象装配与用户消息留在 R 层。
- 性能相关 PR 必须附 benchmark，至少报告运行时间与内存开销；具体规模、阈值和记录格式写入 research plan 或开发文档，不写入本文件。

### Rcpp / C++ 边界

- 允许在后续阶段引入 `Rcpp`，但只用于纯计算核，不得把 C++ 作为唯一可运行路径。
- 编译路径必须保留可验证的串行/纯 R 回退。
- C++ 接口应尽量接收已对齐的 matrix / vector / index 数据，不在 C++ 内重复做高层对象解析。
- 编译优化不得改变结果对象结构和公开 API。

### 并发约束

- 包内并发必须是**可选加速**，无并发支持时必须可靠回退到串行。
- 包内核心并发优先考虑可选 OpenMP；`future`/`furrr` 一类框架仅用于外部实验脚本或 benchmark，不作为包内默认后端。
- 禁止嵌套并行；不得与未知外层并发或多线程 BLAS 叠加造成不可控资源竞争。
- 并行实现必须保证结果可重复，尤其是 top-k 排序和 tie 行为必须稳定。
- 工作线程中不得调用 R API，不得进行消息输出、对象构造或依赖外部可变状态。

## 13. 重构触发条件

出现以下情况时，应优先考虑结构性重构，而不是继续叠加补丁：

- 同类 bug 在多个文件重复出现。
- 新功能需要同时改动多个无关模块。
- 单个文件同时承担校验、编排、计算、输出等多层职责。
- 新贡献者难以判断应该复用哪个 helper。
- 性能优化需要在多个位置复制相同逻辑。

## 14. AI 贡献规则

- AI 修改应优先减少系统复杂度，而不是新增包装层或兼容层。
- 不得通过静默 fallback 掩盖数据问题、对象不一致或算法异常。
- 新增概念、对象或 API 时，必须在 PR 描述中说明动机、边界和影响面。
- 修改热点逻辑时，必须同步补充测试或回归用例。
- 若新 helper 能替代旧 helper，应一并清理旧路径，避免长期并存。

## 15. Git 与 PR 规范

- 提交信息简洁明确，说明改了什么。
- 每个 PR 关联 issue 或在描述中说明动机。
- 不允许直接 push 到 main 分支。
- PR 需通过所有 CI 检查后合并。

PR 提交前至少自查：

- 是否新增了不必要的导出 API。
- 是否引入了重复逻辑或多余 wrapper。
- 是否保持对象契约、维度和命名稳定。
- 是否补充了必要测试与回归用例。
- 是否为性能变更附带了 benchmark 或 research plan 说明。
