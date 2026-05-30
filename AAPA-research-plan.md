# AAPA：锚点辅助的水产全同胞家系系谱归属算法研究方案

## 1. 项目定位

本项目拟提出一种面向水产动物家系育种群体的高效系谱归属算法：

**AAPA：Anchor-Assisted Pedigree Assignment**

中文名称：

**锚点辅助的水产全同胞家系系谱归属算法**

该算法不以完全替代 COLONY 等通用系谱重构软件为目标，而是针对对虾、鱼类、贝类等水产育种中的特定高通量应用场景：

```text
已知父母本
+ 已知家系设计
+ 已知部分全同胞家系成员作为 anchor 个体
+ 大量测试个体需要快速归属到已有家系
```

在这种场景下，系谱重构问题可以从传统的全局无监督搜索转化为：

```text
在已知父母组合和 anchor 家系约束下，
对大量测试个体进行高效、可解释、可拒判的家系归属。
```

## 2. 科学问题

### 2.1 背景问题

水产动物家系育种中，常需要对大量候选个体进行系谱重构或家系归属。传统方法如 COLONY 能够在亲本未知或家系结构未知时进行通用亲缘推断，但在大规模育种应用中存在以下问题：

1. 全局搜索计算量大；
2. 随测试个体数、候选亲本数和 marker 数增加，运行时间快速上升；
3. 对已知父母本和已知家系 anchor 信息的利用不够充分；
4. 对生产育种场景中的快速批量分析不够友好；
5. 对低质量样本、未知来源样本、混样样本需要更明确的拒判机制。

### 2.2 核心科学问题

本研究关注的问题是：

> 在已知父母本和部分家系 anchor 个体的条件下，如何构建一种比通用全局系谱搜索更高效、更可解释、适合大规模水产家系育种的系谱归属算法？

进一步分解为：

1. 如何利用父母本 Mendelian 传递约束快速排除不可能家系？
2. 如何利用 anchor 个体刻画全同胞家系的遗传相似性？
3. 如何将父母本兼容性和 anchor 亲缘相似度整合为稳定的家系归属分数？
4. 如何在保持高准确率的同时显著压缩计算时间？
5. 如何识别未知家系、混样、低质量基因型和家系记录错误？

## 3. 与 COLONY 的关系

### 3.1 COLONY 的优势

COLONY 是成熟的通用系谱重构软件，优势包括：

- 可处理复杂亲缘关系；
- 支持 full-likelihood 方法；
- 可在亲本未知或部分未知时推断 sibship 和 parentage；
- 有较强的统计遗传学基础；
- 在许多物种和场景中已有应用。

### 3.2 COLONY 在本研究场景中的局限

在对虾等水产育种群体中，如果已知父母本和部分家系成员，直接使用通用全局重构存在一定冗余：

```text
COLONY 问题形式：
在未知或部分未知系谱结构下搜索最优全局亲缘配置。

AAPA 问题形式：
在已知父母组合和 anchor 约束下，将测试个体归属到已有家系。
```

AAPA 的目标不是解决更宽泛的问题，而是利用更强先验解决更具体的问题，从而提高计算效率和工程可用性。

## 4. AAPA 的核心思想

AAPA 由四个核心模块组成：

```text
父母本 Mendelian 兼容性
+ anchor 家系亲缘相似度
+ top-k 候选家系剪枝
+ 置信度校准与拒判机制
```

算法基本逻辑：

1. 对每个测试个体，计算其与所有已知父母组合的 Mendelian 冲突率；
2. 根据冲突率筛选 top-k 候选家系；
3. 在候选家系内，计算测试个体与该家系 anchor 个体的亲缘相似度；
4. 综合父母本兼容性、anchor 相似度和数据质量指标形成最终得分；
5. 根据 top1 分数、top1-top2 差值和冲突率阈值进行归属或拒判。

### 4.1 当前方案评估与优化方向

从研究设计角度看，AAPA 的优势在于问题定义足够清晰：它不试图重做通用 pedigree reconstruction，而是把强先验条件转化为一个可剪枝、可解释、可拒判的 supervised assignment 问题。

当前方案最需要补强的不是增加更多复杂模型，而是把以下四个环节前置并制度化：

1. 输入数据质控：父母本、anchor 和测试个体必须先通过一致性检查，否则后续分数不可解释；
2. 阈值校准：`tau_score`、`tau_gap`、`tau_conflict` 不宜只凭经验固定，应由模拟集或留出验证集确定；
3. unknown / ambiguous 评估：拒判机制必须作为主要结果之一评估，而不是作为附加异常处理；
4. 最小闭环实验：第一版应优先证明“父母本冲突率 + top-k + 基础拒判”能稳定工作，再加入 anchor score 和 Rcpp 加速。

因此优化后的开发顺序应为：

```text
数据校验
-> Mendelian conflict baseline
-> top-k 候选家系
-> 阈值校准与拒判
-> anchor score 增益验证
-> Rcpp / OpenMP 加速
```

## 5. 算法输入与输出

### 5.1 输入文件

推荐输入包括：

```text
1. genotype matrix
   行为个体，列为 SNP/marker，取值为 0/1/2/NA

2. parents table
   每个家系对应的父本和母本 ID

3. anchors table
   anchor 个体及其已知家系标签

4. test sample list
   待归属测试个体列表

5. marker metadata
   可选，包括染色体、位置、MAF、缺失率、质量分数等
```

### 5.2 parents table 示例

```text
family_id   father_id   mother_id
F001        P001        P101
F002        P002        P102
F003        P003        P103
```

### 5.3 anchors table 示例

```text
sample_id   family_id
A0001       F001
A0002       F001
A0003       F002
A0004       F002
```

### 5.4 输出结果

每个测试个体输出：

```text
sample_id
assigned_family
status
top1_family
top1_score
top2_family
top2_score
score_gap
mendelian_conflict_top1
anchor_score_top1
missing_rate
effective_marker_count
```

其中 `status` 可包括：

```text
Assigned
Ambiguous
UnknownFamily
LowQualityGenotype
PossibleContamination
PossibleParentRecordError
```

### 5.5 输入数据校验与最小质控

AAPA 的分数解释依赖输入表之间的一致性。正式计算前应先执行最小质控，并输出独立 QC 报告。

必须检查：

1. `parents` 中的 `father_id`、`mother_id` 均存在于 genotype matrix；
2. `anchors` 中的 `sample_id` 均存在于 genotype matrix，且 `family_id` 存在于 parents table；
3. `test_samples` 不与父母本和 anchor 重复，除非显式允许重复验证；
4. genotype matrix 的样本 ID 和 marker ID 唯一；
5. 基因型编码统一为 `0/1/2/NA`，并明确 `0/1/2` 对应同一参考等位基因方向；
6. 每个样本和 marker 的缺失率低于预设阈值；
7. 每个家系至少有可用父母本，anchor 数低于下限的家系应被标记，而不是静默参与完整评分；
8. 每个 anchor 应先计算其与记录父母本的 Mendelian conflict rate，冲突率异常的 anchor 不应直接用于家系原型。

建议输出：

```text
sample_qc table
marker_qc table
family_qc table
anchor_qc table
excluded_samples
excluded_markers
warnings
```

## 6. 方法模型

### 6.1 Mendelian conflict rate

对于测试个体 `i` 和家系 `f`，设家系 `f` 的父母本为：

```text
father_f, mother_f
```

在 marker `m` 上，根据父母本基因型判断子代可能基因型集合：

```text
PossibleGenotypes(f, m)
```

如果测试个体基因型不在可能集合中，则记为一次 Mendelian conflict。

定义：

```text
ConflictRate(i, f) =
  ConflictCount(i, f) / EffectiveMarkerCount(i, f)
```

其中 `EffectiveMarkerCount` 为测试个体和父母本均有可用基因型的 marker 数。

父母本兼容性可以定义为：

```text
ParentCompatibility(i, f) = 1 - normalized ConflictRate(i, f)
```

也可以使用带分型错误率的 likelihood：

```text
MendelianLogL(i, f) =
  sum_m log P(g_i,m | g_father,m, g_mother,m, error_rate_m)
```

第一版建议先使用 conflict rate，便于解释和调试；后续可加入 likelihood 模型。

### 6.2 Anchor kinship score

对于测试个体 `i` 和家系 `f` 的 anchor 集合：

```text
AnchorSet(f) = {a1, a2, ..., ar}
```

计算测试个体与每个 anchor 的相似度：

```text
Similarity(i, a)
```

可选相似度包括：

```text
IBS similarity
KING kinship
GRM relationship
rare allele sharing
```

家系层面的 anchor 分数定义为：

```text
AnchorScore(i, f) =
  mean(top-r Similarity(i, a), a in AnchorSet(f))
```

使用 top-r 而非全部均值，目的是降低个别低质量 anchor 或错误标签 anchor 的影响。

### 6.3 综合归属分数

基础综合分数：

```text
FinalScore(i, f) =
  w1 * ParentCompatibility(i, f)
+ w2 * AnchorScore(i, f)
+ w3 * FamilyPrior(f)
- w4 * MissingPenalty(i)
- w5 * ConflictPenalty(i, f)
```

其中：

- `ParentCompatibility` 表示父母本遗传传递兼容性；
- `AnchorScore` 表示与该家系 anchor 个体的亲缘相似度；
- `FamilyPrior` 可选，来自预期家系规模或抽样比例；
- `MissingPenalty` 对高缺失个体进行惩罚；
- `ConflictPenalty` 对超出合理错误率的 Mendelian conflict 进行惩罚。

如果没有可靠家系规模先验，第一版可令：

```text
w3 = 0
```

第一版不建议把线性权重解释为已经训练好的统计模型。更稳妥的做法是将其视为 ranking score，并通过验证集或模拟数据校准拒判阈值。推荐先使用简单、可复现的默认形式：

```text
ParentScore(i, f) = 1 - min(ConflictRate(i, f) / tau_conflict_high, 1)

FinalScore(i, f) =
  w1 * ParentScore(i, f)
+ w2 * AnchorScoreScaled(i, f)
- w4 * MissingPenalty(i)
```

其中 `AnchorScoreScaled` 应在同一测试个体的候选家系内缩放到可比较范围，避免 IBS、KING 或 GRM 分数尺度不同导致权重不可解释。

### 6.4 top-k 候选家系剪枝

对每个测试个体，先用 Mendelian conflict rate 对所有家系排序，只保留 top-k：

```text
CandidateFamilies(i) =
  top-k families with lowest ConflictRate(i, f)
```

后续 anchor score 和综合打分只在候选家系内计算。

这一步是 AAPA 提高效率的关键：

```text
父母本全量筛选：N * F * M
朴素 anchor 精细评分：N * k * A * M
预计算或 bitset 加速后 anchor 评分：约 N * k * A
```

其中 `k << F`。因此 top-k 剪枝减少的是精细评分所需的候选家系数；如果 anchor 相似度仍逐 marker 现算，则复杂度中仍然包含 `M`。正式 benchmark 应分别报告朴素实现和 bitset / 预计算实现的速度。

### 6.5 拒判规则

对测试个体 `i`，设最终得分最高的家系为 `f1`，第二高为 `f2`。

归属条件：

```text
FinalScore(i, f1) >= tau_score
FinalScore(i, f1) - FinalScore(i, f2) >= tau_gap
ConflictRate(i, f1) <= tau_conflict
AnchorScore(i, f1) >= tau_anchor
MissingRate(i) <= tau_missing
EffectiveMarkerCount(i, f1) >= tau_marker
```

若不满足，则根据失败原因标记为：

```text
Ambiguous
UnknownFamily
LowQualityGenotype
PossibleContamination
PossibleParentRecordError
```

### 6.6 阈值校准策略

拒判阈值应通过数据驱动方式确定。建议至少准备一个 calibration set，可以来自模拟数据、半模拟数据或真实数据中的留出家系。

推荐校准目标：

```text
tau_conflict: 控制真实同家系个体的误拒率
tau_gap: 控制近缘家系之间的误分配率
tau_anchor: 控制 unknown family 被强制分配的比例
tau_missing: 控制低质量样本进入归属流程的比例
tau_marker: 保证每个判断有足够有效 marker 支撑
```

阈值选择不应只追求最高 assignment accuracy，而应同时优化：

```text
assigned accuracy
unknown detection AUROC
ambiguous precision
false assignment rate
rejection rate
```

论文中应报告不同阈值下的 trade-off 曲线，尤其是 `false assignment rate` 与 `rejection rate` 的关系。这比单一准确率更符合育种生产场景，因为错误归属通常比拒判更严重。

## 7. 可选增强模块

### 7.1 Marker informativeness weighting

不同 SNP 对区分家系的贡献不同。可以根据以下指标给 marker 加权：

- MAF；
- 父母本之间的区分度；
- 家系间 allele frequency 差异；
- 缺失率；
- 分型错误率；
- 与其他 SNP 的 LD 程度。

加权后：

```text
WeightedConflictRate(i, f) =
  sum_m weight_m * conflict_m / sum_m weight_m
```

### 7.2 Family-specific error model

不同家系、批次或测序芯片可能有不同错误率。可建立：

```text
error_rate_m
error_rate_batch
error_rate_family
```

用于更准确地计算 Mendelian likelihood。

### 7.3 Anchor robustness weighting

对 anchor 个体赋予质量权重：

```text
AnchorWeight(a) =
  function(missing_rate, genotype_quality, consistency_with_parents)
```

家系 anchor 分数改为加权 top-r 均值。

### 7.4 图传播模块

如果基础 AAPA 在近亲家系或高缺失数据中区分力不足，可引入轻量图传播：

```text
节点：测试个体、anchor 个体、家系
边：IBS/kinship 相似度、候选家系边、已知 anchor 标签边
方法：label propagation 或 personalized PageRank
```

该模块作为增强，不作为第一版必需组件。

### 7.5 深度学习模块

若数据规模足够大且 anchor 标签充分，可尝试：

- genotype embedding；
- contrastive learning；
- GraphSAGE；
- MLP classifier；
- family prototype network。

但深度学习应作为后续研究方向，第一阶段不建议依赖它作为核心算法。

## 8. 复杂度分析

设：

```text
N = 测试个体数
F = 已知家系数
M = marker 数
A = 每个家系平均 anchor 数
k = 每个测试个体保留的候选家系数
```

### 8.1 父母本冲突率计算

```text
O(N * F * M)
```

该步骤可高度并行：

- 按测试个体并行；
- 按家系并行；
- 按 marker 分块；
- 使用 bitset/SIMD；
- 使用 Rcpp/OpenMP。

### 8.2 anchor 精细评分

如果只对 top-k 家系计算：

```text
O(N * k * A * M)
```

实际可通过预计算 anchor 表征或 kinship 降低成本：

```text
O(N * k * A)
```

### 8.3 最终排序与拒判

```text
O(N * k)
```

### 8.4 与 COLONY 的区别

AAPA 的优势不是单个 likelihood 公式一定更简单，而是避免全局组合搜索：

```text
COLONY：搜索可能的全局 sibship / parentage 配置
AAPA：在已知家系集合内进行个体级候选排序
```

因此在已知父母本和 anchor 的大规模育种场景中，AAPA 有望显著压缩运行时间。

## 9. 实验设计

### 9.1 数据类型

建议使用三类数据：

1. 模拟数据；
2. 半模拟数据；
3. 真实对虾或其他水产育种数据。

### 9.2 模拟数据

模拟因素包括：

```text
家系数：50, 100, 200, 500
每家系测试个体数：20, 50, 100, 500
SNP 数：500, 1k, 5k, 10k, 50k
anchor 数：1, 2, 5, 10, 20
分型错误率：0.1%, 0.5%, 1%, 2%, 5%
缺失率：1%, 5%, 10%, 20%
近亲父母本比例：低、中、高
unknown family 混入比例：0%, 5%, 10%, 20%
```

### 9.3 半模拟数据

从真实父母本基因型出发，根据 Mendelian 传递模拟子代，保留真实群体的 allele frequency、LD 结构和缺失模式。

半模拟数据的优势是更接近真实应用场景。

### 9.4 真实数据

真实数据用于最终验证：

- 对虾家系育种数据；
- 已知父母本；
- 部分已知家系成员；
- 大量待测个体；
- 若有人工记录或历史系谱，可作为参考 truth。

### 9.5 对照方法

至少包括：

1. COLONY；
2. KING / PLINK 最近亲缘归属；
3. 仅父母本 Mendelian conflict；
4. 仅 anchor kinship；
5. AAPA 完整模型。

可选：

6. AAPA + marker weighting；
7. AAPA + graph propagation；
8. AAPA + likelihood error model。

### 9.6 评价指标

准确性：

```text
assignment accuracy
top-k accuracy
precision
recall
F1-score
unknown detection AUROC
ambiguous detection rate
```

效率：

```text
runtime
memory usage
scalability with N
scalability with F
scalability with M
speedup relative to COLONY
```

可靠性：

```text
calibration curve
confidence-score relationship
error rate under missing data
error rate under genotype error
performance under related parents
```

### 9.7 消融实验

必须设计消融实验，证明每个模块的贡献：

```text
AAPA without anchor score
AAPA without parent compatibility
AAPA without top-k pruning
AAPA without rejection
AAPA with / without marker weighting
AAPA with different anchor numbers
```

### 9.8 最小闭环实验

在进入大规模 benchmark 前，建议先完成一个最小闭环实验，用于验证核心假设是否成立。

建议设置：

```text
家系数：20-50
每家系 anchor 数：3-5
每家系测试个体数：20-50
SNP 数：1k-5k
分型错误率：0.5%-2%
缺失率：1%-10%
unknown family 比例：5%-10%
```

必须输出：

```text
conflict-rate-only baseline
anchor-only baseline
AAPA combined score
top-k recall
assigned accuracy
unknown detection AUROC
false assignment rate
runtime
```

进入下一阶段的最低判据可以设为：

```text
top-k recall 明显高于 top1 recall
AAPA combined score 优于任一单模块 baseline
unknown 样本不会被大量强制归属
运行时间随测试个体数近似线性增长
```

## 10. 软件实现路线

### 10.1 技术选择

建议采用：

```text
R package + Rcpp
```

理由：

- Rcpp 生态成熟；
- 便于快速形成 R 包；
- 适合论文实验、统计分析和可视化；
- 与 testthat、roxygen2、devtools 等 R 包开发工具兼容；
- 后续可将核心模块进一步迁移到 Rust 或独立 CLI。

### 10.2 推荐包名

候选包名：

```text
aapa
aquaped
anchorPed
pedAAPA
```

建议首选：

```text
aapa
```

### 10.3 R 包结构

```text
aapa/
  DESCRIPTION
  NAMESPACE
  R/
    assign.R
    qc.R
    scoring.R
    simulate.R
    benchmark.R
    plot.R
    io.R
  src/
    mendelian.cpp
    anchor_kinship.cpp
    topk.cpp
    utils.cpp
  inst/
    extdata/
      small_example/
  tests/
    testthat/
      test-mendelian.R
      test-anchor-score.R
      test-assignment.R
      test-rejection.R
  vignettes/
    aapa-workflow.Rmd
    benchmark-colony.Rmd
```

### 10.4 R 用户接口

核心函数：

```r
result <- aapa_assign(
  geno = geno_matrix,
  parents = parents_table,
  anchors = anchors_table,
  test_samples = test_ids,
  top_k = 5,
  threads = 8,
  conflict_threshold = 0.02,
  score_gap_threshold = 0.10
)
```

其中 `threads` 必须是**可选加速参数**：

- `threads = 1` 或无并发支持时，必须可靠回退到串行；
- 并发仅影响性能，不得改变结果对象结构、排序语义和拒判语义；
- 若后续增加 compiled backend，用户层接口保持稳定，不暴露底层实现细节。

输出为 data frame：

```r
sample_id
assigned_family
status
top1_family
top1_score
top2_family
top2_score
score_gap
mendelian_conflict
anchor_score
missing_rate
effective_marker_count
```

### 10.5 MVP 版本

第一版最小可用功能：

1. 输入 dosage matrix；
2. 输入 parents table；
3. 输入 anchors table；
4. 计算 Mendelian conflict rate；
5. 输出每个测试个体 top-k 候选家系；
6. 设置 conflict threshold 拒判。

该版本可以先用 R 实现，确保公式和结果正确。

### 10.6 Rcpp 加速版本

第二版将瓶颈迁移到 Rcpp：

- Mendelian conflict 计算；
- top-k 排序；
- IBS / anchor similarity；
- 缺失值统计；
- 多线程并行。

建议优先加速：

```text
N * F * M 的 conflict rate 计算
```

这是全流程最可能的计算瓶颈。

#### 10.6.1 实现边界

Rcpp 版本的边界应尽早固定，避免后续接口反复重构：

- R 层负责输入校验、ID 对齐、对象构造、拒判规则和用户消息输出；
- C++ 层只负责纯计算核，如 Mendelian conflict、anchor similarity、top-k 选择和缺失计数；
- 传入 C++ 的对象应尽量是已对齐的 matrix / vector / index，不在 C++ 内重复解析高层 data.frame / list 语义；
- 编译路径不得成为唯一可运行路径，必须始终保留纯 R 参考实现。

#### 10.6.2 benchmark 设计原则

性能 benchmark 不只是报告“更快”，还必须回答“是否仍然正确、可复现、可扩展”。

固定原则：

1. **profile-first**：先用纯 R 版本做 profiling，再决定优化哪一个核；
2. **正确性先于速度**：每个 compiled kernel 必须与纯 R truth oracle 做数值等价比较；
3. **分层 benchmark**：分别测试单核函数、完整流程、不同 `top_k`、不同 anchor 数和不同 missing rate；
4. **同时报告时间与内存**：避免只看运行时间而忽略峰值内存；
5. **记录环境信息**：R 版本、编译器、CPU、线程数、BLAS 环境都应写入 benchmark 结果。

#### 10.6.3 benchmark 场景矩阵

建议固定三档 benchmark 场景，用于开发期和论文期复现：

```text
Small:
  families = 20
  offspring_per_family = 50
  anchors_per_family = 3
  SNPs = 1,000

Medium:
  families = 50
  offspring_per_family = 200
  anchors_per_family = 5
  SNPs = 10,000

Large:
  families = 100
  offspring_per_family = 1,000
  anchors_per_family = 5-10
  SNPs = 50,000-100,000
```

每档场景至少再做以下扰动：

- missing rate：1%、5%、10%；
- error rate：0.1%、0.5%、1%；
- unknown family 比例：0%、5%、10%；
- `top_k`：1、3、5、10；
- anchor 数：1、3、5、10。

#### 10.6.4 benchmark 指标模板

每次 benchmark 至少输出以下字段：

```text
scenario_id
n_individuals
n_families
n_markers
n_anchors
missing_rate
error_rate
top_k
backend            # pure-r / rcpp-serial / rcpp-openmp
threads
elapsed_sec
peak_mem_mb
assigned_accuracy
rejection_rate
false_assignment_rate
unknown_detection_auroc
top1_recall
topk_recall
result_identical   # TRUE/FALSE 或 tolerance pass/fail
```

开发期建议把结果整理成两张固定表：

1. **correctness table**：比较 pure R 与 compiled backend 的结果一致性；
2. **performance table**：比较不同 backend 与线程数下的 runtime / memory。

#### 10.6.5 数值等价与容差

推荐把纯 R 版本作为 truth oracle，并采用以下比较口径：

- 整数结果（例如计数、rank、top-k family ID）要求严格一致；
- 浮点结果（例如 conflict rate、anchor score、composite score）允许在预设容差内一致；
- `NA` 传播语义必须完全一致；
- top-k 的 tie 行为必须固定，不能因并发或 partial sort 改变输出顺序；
- `status` 与 `reject_reason` 必须与纯 R 路径解释一致。

具体容差阈值、比较函数和失败时的诊断脚本建议维护在 benchmark 文档或测试辅助文件中，而不是分散到多个实验脚本。

#### 10.6.6 并发策略

包内并发建议遵循“**默认串行，可选多线程**”的路线：

- 首选可选 OpenMP，对外仅表现为 `threads` 参数；
- 无 OpenMP 支持时自动退回串行，不报伪错误；
- `future` / `furrr` 更适合外部实验脚本，不建议作为包内默认 backend；
- 不建议在包内核心路径优先引入 `RcppParallel`，除非 OpenMP 兼容性或维护成本明显不可接受；
- 禁止嵌套并行，避免与用户外层并发、BLAS 多线程叠加。

建议优先并行的维度：

```text
按 test individual 并行
或按 family block 并行
```

这两种方式最容易保持任务独立、写回位置固定，便于保证结果可重复。

#### 10.6.7 线程安全约束

并发实现必须满足：

- 工作线程中不调用 R API；
- 工作线程中不进行消息输出、S3 对象构造或依赖外部可变状态；
- 所有输出缓冲区预先分配，并按固定索引写回；
- 若使用 partial sort / top-k 选择，必须显式定义 stable tie-break 规则；
- benchmark 必须覆盖 `threads = 1, 2, 4, 8` 等典型设置，并验证结果完全一致。

#### 10.6.8 阶段性交付物

Rcpp / 并发阶段的最小交付物不只是代码，还应包括：

```text
1. 纯 R 与 compiled backend 的对照测试
2. Small / Medium / Large 三档 benchmark 结果
3. 串行 vs 多线程 的速度与内存报告
4. 线程数与加速比曲线
5. 已知限制与回退策略说明
```

### 10.7 后续 Rust / CLI 可能性

如果后续需要产业级部署或超大规模分析，可新增：

```text
Rust backend
standalone CLI
```

但第一阶段不建议同时维护 Rcpp 和 Rust，避免分散研发精力。

## 11. 预期创新点

### 11.1 场景创新

将水产全同胞家系育种中的系谱分析问题，从通用无监督重构重新定义为：

```text
已知父母本和 anchor 约束下的大规模家系归属问题。
```

### 11.2 方法创新

提出父母本 Mendelian 兼容性与 anchor 亲缘相似度的联合评分框架。

### 11.3 计算创新

提出 top-k 候选家系剪枝策略，大幅减少精细评分和后续优化的计算量。

### 11.4 应用创新

面向对虾等水产育种群体，提供可解释、可拒判、可批量运行的系谱归属方法。

### 11.5 软件创新

开发 R 包实现，提供从数据输入、质控、归属、benchmark 到可视化的完整研究工具链。

## 12. 论文框架建议

### 12.1 题目候选

英文：

```text
AAPA: Anchor-Assisted Pedigree Assignment for Large-Scale Full-Sib Family Reconstruction in Aquaculture Breeding
```

中文：

```text
面向水产全同胞家系育种的锚点辅助高效系谱归属算法
```

### 12.2 摘要逻辑

1. 水产育种需要大规模系谱归属；
2. 通用系谱重构方法在大规模应用中计算成本较高；
3. 已知父母本和 anchor 个体提供了强先验；
4. 本研究提出 AAPA；
5. AAPA 结合 Mendelian 兼容性、anchor 亲缘相似度和 top-k 剪枝；
6. 在模拟和真实数据上比较准确率与速度；
7. AAPA 在保持准确率的同时显著降低运行时间。

### 12.3 主要章节

```text
Introduction
Materials and Methods
  Study scenario
  AAPA algorithm
  Simulation design
  Real dataset
  Comparative methods
  Evaluation metrics
Results
  Accuracy comparison
  Runtime comparison
  Effect of SNP number
  Effect of anchor number
  Robustness to missingness and genotyping error
  Unknown family detection
Discussion
  Practical implications
  Limitations
  Future extensions
Conclusion
```

## 13. 阶段计划

### 阶段 1：算法原型

目标：

- 完成 R 原型；
- 验证 Mendelian conflict 计算；
- 输出 top-k 家系；
- 建立小规模模拟数据。

成果：

```text
R scripts
toy dataset
初步准确率和运行时间
```

### 阶段 2：AAPA 完整模型

目标：

- 加入 anchor kinship；
- 加入综合评分；
- 加入拒判；
- 完成中等规模模拟 benchmark。

成果：

```text
AAPA v0.1
模拟实验结果
消融实验结果
```

### 阶段 3：Rcpp 加速

目标：

- 将 conflict rate 和 anchor score 核心计算迁移到 Rcpp；
- 支持可选多线程；
- 优化内存使用。

成果：

```text
AAPA R package prototype
速度 benchmark
纯 R vs compiled 一致性报告
串行回退说明
```

### 阶段 4：真实数据验证

目标：

- 在真实对虾或水产育种数据上验证；
- 与 COLONY / KING / PLINK baseline 比较；
- 评估 unknown 和 ambiguous 样本。

成果：

```text
真实数据结果
论文图表
软件说明文档
```

### 阶段 5：论文和软件发布

目标：

- 完成论文；
- 整理 R 包；
- 编写 vignette；
- 发布 GitHub 版本。

成果：

```text
manuscript
R package
example dataset
user guide
```

## 14. 风险与应对

### 14.1 父母本记录错误

风险：

```text
错误父母本会导致真实子代出现高 Mendelian conflict。
```

应对：

- 设置 PossibleParentRecordError 状态；
- 检查 anchor 与父母本的一致性；
- 支持父母本替代候选检测。

### 14.2 anchor 数量不足

风险：

```text
家系原型不稳定，anchor score 噪声较大。
```

应对：

- 分析不同 anchor 数的性能；
- 使用 top-r 聚合；
- 对 anchor 做质量加权；
- 当 anchor 不足时提高父母本兼容性权重。

### 14.3 近亲父母本导致家系难区分

风险：

```text
top1 和 top2 家系分数接近。
```

应对：

- 使用 score gap 拒判；
- 输出 Ambiguous；
- 加入更多高信息 SNP；
- 使用 marker weighting。

### 14.4 分型错误和缺失率高

风险：

```text
冲突率升高，误判或拒判增加。
```

应对：

- 建立错误率模型；
- 设置有效 marker 数阈值；
- 使用 likelihood 代替硬冲突；
- 对低质量样本标记 LowQualityGenotype。

### 14.5 unknown family 混入

风险：

```text
未知来源个体被强行分配到已知家系。
```

应对：

- 设置 UnknownFamily 类别；
- 使用 top1 score 和 score gap；
- 使用 anchor score 下限；
- 在实验中专门评估 unknown detection。

## 15. 推荐的近期工作

近期优先完成以下任务：

1. 固定 AAPA 的输入输出格式和 QC 报告格式；
2. 准备一个带 truth label、unknown family、缺失和分型错误的小型模拟数据集；
3. 用 R 实现 Mendelian conflict rate，并输出 conflict matrix；
4. 输出每个测试个体的 top-k 家系候选，并计算 top-k recall；
5. 实现基础拒判和阈值校准脚本；
6. 加入 anchor-only 和 conflict-only baseline；
7. 设计 benchmark 表格模板，固定 accuracy、rejection、runtime 三类指标；
8. 在确认瓶颈后再决定 Rcpp 加速接口。

建议把第 7-8 项具体化为：

- 固定 `Small / Medium / Large` 三档 benchmark 数据规模；
- 固定 `pure-r / rcpp-serial / rcpp-openmp` 三类 backend 标签；
- 固定 `threads = 1, 2, 4, 8` 的测试矩阵；
- 固定 correctness table 与 performance table 两类输出；
- 先完成串行 compiled kernel，再增加并发，避免把“编译优化”和“并发优化”混为一次变更。

第一阶段不要同时引入深度学习、图神经网络和复杂文件格式解析，应优先把核心科学问题和基准实验跑通。

## 16. 总结

AAPA 的核心价值在于：针对水产家系育种中“已知父母本 + 已知 anchor + 大量测试个体”的现实场景，提出一种比通用全局系谱重构更直接、更高效、更可解释的系谱归属算法。

该研究路线兼顾科学创新和工程落地：

```text
科学上：
提出 anchor-assisted pedigree assignment 框架。

方法上：
整合 Mendelian 兼容性、anchor 亲缘相似度、top-k 剪枝和拒判机制。

工程上：
开发 R package + Rcpp 原型，支持大规模数据分析和论文实验复现。
```

推荐将 AAPA 作为主算法，将 COLONY、KING、PLINK 等作为对照和 baseline，从而形成具有明确应用场景、算法贡献和软件支撑的科研成果。

