---
layout: post
title: "ChIP-seq FRiP：计算方法、归一化伪影与结果验证"
date: 2026-10-20
description: "从 paired-end fragment 与 FRiP 计算，到 NSD1、DOT1L 抑制、全局组蛋白修饰变化及 ChIP-seq 归一化伪影的辨析"
tags: [技术, 生物信息学, ChIP-seq, FRiP, 归一化, 表观遗传]
---

FRiP（**Fraction of Reads in Peaks**）是 ChIP-seq、CUT&Tag 等富集类测序常用的质量控制指标，用于描述有多少测序 reads（或 DNA fragments）落在被识别为显著富集的 **peaks** 内。

## 1. FRiP 到底表示什么？

最直观的理解是：

> **所有有效测序片段中，有多大比例位于 peak calling 识别出的富集区域？**

计算公式：

```text
FRiP = Fragments overlapping peaks / Total usable fragments

FRiP (%) = FRiP × 100%
```

假设一个**完全虚构**的示例样本中，有 1,000 个通过质控的 DNA fragments，其中 300 个与 peak 区域重叠，则：

```text
FRiP = 300 / 1,000 = 0.30 = 30%
```

这表示 30% 的有效片段落在 peak 区域内。**它不是“30% 的基因被修饰”，也不是“组蛋白修饰丰度为 30%”。**

## 2. Paired-end 中 fragment 是什么？

对于 paired-end 测序，一个 DNA fragment 对应同一个插入片段产生的 **R1/R2 两端 reads**。通过它们在参考基因组上的比对坐标，可以推断完整插入片段的区间。

```text
Genome:    --------------------------------------------------->

R1:            [>>>>>>]
R2:                                  [<<<<<<]

Fragment:      [============================]
                |                          |
               start                      end
```

上图中，R1 与 R2 中间的区域即使没有被直接测序，也属于根据 paired-end 比对推断的 fragment 区间。

这里的“合并”是**根据坐标重建 DNA 插入片段的范围**，而不是把两条 reads 的碱基序列直接拼接起来。

**计数单位非常重要：**

- 如果按 reads 计数，R1 和 R2 通常是两条 reads。
- 如果按 paired-end fragments 计数，一个正确配对的 R1/R2 只算 **1 个 fragment**。
- 分子和分母必须采用一致的单位与过滤规则，不能用“位于 peaks 内的 reads / 全部 fragments”。

尽管 FRiP 名字中包含 *Reads*，在 paired-end 分析中也经常以 fragments 为单位计算；报告结果时应明确采用的是哪一种口径。

## 3. 怎样判断 fragment 是否位于 peak 中？

只要一个 fragment 与任意 peak 有重叠（例如至少 **1 bp**），就可以视为一个 **fragment in peaks**。

```text
Peak A:       [-------------]
Peak B:                              [-----------]

Fragment 1:       [=======]                         yes
Fragment 2:                     [=======]           no
Fragment 3:                           [========]    yes
```

特别注意：**一个 fragment 即使同时覆盖两个相邻或重叠的 peaks，也只能计数一次。** 实际计算时通常先合并互相重叠的 peak 区间，或者使用“有任意重叠即记录一次”的去重逻辑。

因此，FRiP 的分子不是 peak 数量，也不是覆盖碱基的总数，而是**与一个或多个 peaks 重叠的有效片段个数**。

## 4. 计算 FRiP 的基本流程

1. **Peak calling**：从 ChIP-seq 数据得到 peak 区间，例如 BED 或 narrowPeak。不同实验目标可能使用 broad-peak 或 narrow-peak 方法。
2. **BAM 过滤**：预先定义有效片段，例如保留正确配对、唯一或高质量比对的片段，排除 PCR duplicates、secondary/supplementary alignments 和 QC-fail reads。
3. **统计分母**：计算所有符合过滤规则的有效 fragments 总数。
4. **统计分子**：计算这些 fragments 中有多少与至少一个 peak 重叠；每个 fragment 最多计数一次。
5. **计算比例**：用分子除以分母，必要时乘以 100% 转为百分比。

示例伪代码：

```text
usable_fragments = filter_and_reconstruct_fragments(BAM)
peak_regions = merge_overlapping_peaks(PeakBED)

total = count(usable_fragments)
in_peaks = count(fragments_overlapping_any_peak(usable_fragments, peak_regions))

FRiP = in_peaks / total
```

**可重复性的关键**是明确记录：peak calling 工具及阈值、BAM 过滤参数、paired-end 计数口径、peak 坐标格式，以及采用的是每个重复样本自己的 peaks 还是一套共同的 peaks。

## 5. FRiP 高或低意味着什么？

FRiP 反映测序片段向**被识别的富集区域**集中的程度，可辅助判断 ChIP 信号与背景的分离情况。但它**不能单独决定实验是否合格**。

影响 FRiP 的因素包括：

- **ChIP 富集效率与背景噪声**：背景较高时，reads 可能更分散。
- **Peak calling 参数**：峰越多、覆盖的基因组区域越大，FRiP 往往更容易升高。
- **目标的基因组分布**：狭窄的转录因子结合位点与覆盖大片染色质的 histone marks 具有不同的 peak 宽度。
- **数据过滤方式**：MAPQ、duplicates、配对过滤等都会改变分母和分子。
- **真实生物学变化**：处理后信号重新分布，甚至广泛丢失富集，都可能改变 FRiP。

因此，**broad marks 与 narrow marks 的 FRiP 不适合直接作数值高低排名**。不同数据集、不同软件及阈值算出的 FRiP 也未必可比。

## 6. FRiP 能否证明组蛋白修饰整体升高或下降？

**不能。**

例如在一个纯教学场景中，某处理可能使修饰在整个基因组中广泛下降，导致可识别 peaks 减少、FRiP 降低。但反过来，FRiP 降低也可能由 peak calling 阈值、背景噪声、有效片段过滤和区域分布变化造成。

更重要的是，FRiP 是**样本内部的相对富集比例**，并不是每个细胞中组蛋白修饰总量的绝对测量。

如果要支持“某种 histone modification 全局减少”的结论，可结合 Western blot、定量质谱或使用外源 spike-in 校准的 ChIP-seq；同时检查信号轨迹、重复一致性和适当归一化后的区域信号。

## 7. 拓展案例：NSD1 缺失、DOT1L 抑制与 ChIP-seq 的异常 peaks

下面讨论的是**公开文献支持的生物学机制及通用分析思路**，不是对任何特定样本或未公开实验结果的描述。

### 7.1 NSD1 KO 为什么可能出现 H3K36me2 loss 与 H3K27me3 gain？

NSD1 是催化 H3K36me2 的甲基转移酶之一。H3K36me2 与 PRC2 催化的 H3K27me3 存在拮抗关系，因此在部分细胞模型中，NSD1 缺失会导致：

- **H3K36me2 loss**：原本受到 NSD1 维持的部分染色质区域失去 H3K36me2。
- **H3K27me3 gain**：H3K36me2 形成的拮抗作用减弱后，H3K27me3 可能出现在新区域，或者向既有 Polycomb 结构域周围扩展（*spreading*）。
- **修饰交界处的重排**：在相同基因组区域看到 H3K36me2 下降、H3K27me3 上升，是值得进一步检验的空间关联，但不自动证明两者存在直接因果关系。

这里需要区分两个层面：**局部 peak gain** 与 **全局修饰总量增加** 不是同一个概念。H3K27me3 的 spreading 可能使信号从较集中的 peaks 转向较广泛的区域；是否表现为 FRiP 下降，取决于 peak caller 是否覆盖新增信号、峰宽、过滤阈值等因素。因此不能仅凭 FRiP 下降就认定发生了 spreading。

这种 NSD1–H3K36me2–H3K27me3 拮抗关系已在不同细胞模型中被报道，但**不同细胞类型和处理时长可能产生不同结果**，仍需在具体模型中验证。

### 7.2 EPZ（DOT1L 抑制剂）为什么可能同时出现 H3K79me2 loss 和 apparent gain？

DOT1L 负责催化 H3K79 甲基化；EPZ-5676 等 DOT1L 抑制剂可以降低细胞内 H3K79me2 水平。公开研究已用 Western blot 和外源参考校准的 ChIP-seq（**ChIP-Rx**）证明，传统相对归一化有时会掩盖这种全局性下降。

一个典型风险是：当大量强富集区域的 H3K79me2 真实减少时，其余相对稳定、微弱富集或接近背景的区域，可能在测序文库中占据更高的**相对比例**。于是用总 mapped reads 进行 CPM/RPM 归一化后，个别区域可能表现为正向 log2FC，即 **apparent gain（表观上调）**。

```text
真实每细胞信号（未知绝对尺度）：

强富集区域：  DOT1L inhibition → H3K79me2 loss
弱富集区域：  可能不变、轻微减少或真实重分布
                          ↓
         ChIP 后测序与文库组成变化
                          ↓
         Total-library normalization
                          ↓
     部分弱信号区域可能出现 positive log2FC
                          ↓
              需要独立验证
```

**注意：** “可能出现”不等于所有 gain peaks 都是伪影。真实的区域重分布、抗体非特异性结合、测序背景、peak 检测阈值和统计模型偏差，均是需要分别评估的解释。

CPM/RPM 直接按文库深度缩放；DESeq2 的 median-of-ratios 与 edgeR 的 TMM 使用不同的组成校正策略，**不能简单视为同一种算法**。当绝大多数分析特征受到同向变化时，常见归一化方法的参考不变假设可能不适用，但 TMM 和 median-of-ratios 并非在任何全局变化下都会失效。

### 7.3 FRiP 能提供什么线索？哪些结论不能推导？

如果 DOT1L 抑制后 H3K79me2 的 FRiP 下降，可以将它作为 **富集峰内片段占比下降** 的线索，与全局修饰减少的机制假设一起考虑。

但不能进行以下推断：

- **不能把 `1 − FRiP` 直接叫作“非特异性背景比例”**：peak 之外也可能有真实但未达到阈值的修饰信号；peak 集合本身也可能随着处理改变。
- **不能用两组 `1 − FRiP` 的比值直接预测特定 locus 的 log2FC**：该比值只概括各自 peak 外的片段比例，不包含该区域的读段分布、峰宽、背景、归一化系数等信息；数值偶然接近也不构成证明。
- **不能由 FRiP 下降直接推出全局绝对丰度下降**：FRiP 缺少每细胞修饰丰度的外部尺度。
- **不能仅凭 Menin 等共因子 FRiP 下降认定蛋白从染色质解离**：还可能涉及 ChIP 效率、背景噪声和峰集合差异，需要结合其他证据。

### 7.4 如何用生信分析验证 normalization artifact？

建议按照由易到难的顺序：

1. **先确认是否存在真正的 spike-in**：核查实验设计与原始数据中是否加入经定量控制的外源细胞、染色质或 DNA。若满足参考标定的实验前提，可用 spike-in 计算 scaling factor 并重新比较区域信号。只有怀疑数据存在 spike-in 才应尝试该方法；不能事后把污染 reads 当作有效 spike-in。
2. **使用固定 genomic regions 重新计数**：分别建立 H3K36me2、H3K27me3、H3K79me2 的统一 peak union 或全基因组 bins；所有条件在**完全相同区域**计数，不仅比较每个样本自己调用的 peaks。这样能减少 peak 集合变化带来的偏差。
3. **比较不同归一化策略**：保留原始 fragment counts，分别检查总文库深度 CPM、DESeq2 size factors、TMM 等方法下的 log2FC 分布。重点观察所谓 gain loci 是否只在某一种归一化方法下出现。若无可靠外源尺度，这属于**敏感性分析，不等于证明绝对变化**。
4. **比较 MA plot 和区域分布**：在同一固定区域集合上画 MA plot、log2FC boxplot 和信号散点图，检查弱信号区域是否更容易出现正向波动；同时统计 peaks 数量、union peak 覆盖长度及 FRiP。
5. **检查 IGV 或 genome browser 轨迹**：挑选具有代表性的 H3K79me2 loss、apparent gain、相对稳定区域，采用一致的归一化口径、明确的 Y 轴尺度，并结合 Input 轨迹判断是否接近背景。未归一化的 raw counts 也应考虑不同样本的有效测序深度。
6. **做独立验证**：如果科学问题是“每个细胞中的组蛋白修饰是否整体减少”，优先使用 Western blot、定量质谱或经过合理实验设计的 spike-in ChIP-seq。若未加 spike-in，单靠常规 ChIP-seq 不能可靠恢复全局绝对尺度。

### 7.5 结果怎样才具有说服力？

| 观察 | 可以支持的解释 | 仍然不能证明 |
| --- | --- | --- |
| NSD1 KO 后，部分区域 K36me2 降低、K27me3 增加 | 局部组蛋白修饰拮抗和空间重排 | 所有区域都出现相同变化 |
| K27me3 peak 边缘或邻近 bins 信号扩展 | 与 H3K27me3 spreading 相符 | FRiP 单独可以确认 spreading |
| DOT1L 抑制后 H3K79me2 FRiP 下降 | 富集在已识别 peaks 内的片段比例减少 | 全局每细胞 H3K79me2 绝对下降 |
| 弱信号区域的 gain 对归一化方法敏感 | 这些 gain 可能具有组成偏差或背景贡献 | 所有 gain 都是假的 |
| 外源校准后原有 gain 消失或反向 | 支持原结果受到归一化选择影响 | 仅凭一次重算就能完全排除实验偏差 |

**核心原则：先用固定区域确认变化是否稳定，再分析归一化是否改变方向；最终用适当的外部参照或独立实验判断全局绝对变化。**

### 相关公开文献

- [Drosos et al., Molecular Cell, 2022](https://doi.org/10.1016/j.molcel.2022.04.015)：NSD1 缺失与 H3K36me2–H3K27me3 拮抗及 Polycomb 域扩展。
- [Orlando et al., Cell Reports, 2014](https://doi.org/10.1016/j.celrep.2014.10.018)：利用 EPZ-5676 与 ChIP-Rx 演示传统 ChIP-seq 归一化对全局 H3K79me2 loss 的遮蔽。
- [DOT1L inhibition reveals a distinct subset of enhancers dependent on H3K79 methylation, Nature Communications, 2019](https://doi.org/10.1038/s41467-019-10844-3)：DOT1L 抑制后 H3K79me2/3 的全局变化及外源参考验证。

---

## 8. FRiP 与其他 ChIP-seq 质控指标有什么不同？

| 指标 | 主要回答的问题 |
| --- | --- |
| Mapping rate | 有多少 reads 可以比对到参考基因组？ |
| Duplicate rate | 有多少 reads 可能属于重复测序或 PCR 扩增片段？ |
| Peak number | 在设定阈值下鉴定到多少富集区域？ |
| FRiP | 有效 reads/fragments 中有多大比例落在 peaks 内？ |
| Replicate concordance | 生物学重复的富集模式是否一致？ |

**总结：FRiP 是“有效测序片段落入 peaks 的比例”，适合用作富集质量的辅助证据；它不是绝对修饰丰度，也不能替代重复一致性或全局定量实验。**

*注：本文数值、坐标示意与伪代码均为教学示例，不对应任何真实研究样本或未公开课题数据。*
