---
layout: post
title: "RcisTarget：从一组基因反推潜在转录调控因子"
date: 2026-09-10
description: "从 motif ranking、AUC 和 NES 理解 RcisTarget 的转录因子调控富集分析流程"
tags: [技术, 生物信息学, 转录因子, RcisTarget, 多组学]
---

RcisTarget 通过分析候选基因调控区域中的 transcription factor binding motif 富集情况，寻找可能共同调控一组基因的转录因子。

## 1. RcisTarget 是做什么的？

在转录组分析中，我们经常会得到一组具有共同生物学特征的基因，例如：

* 差异表达基因 DEG
* 共表达模块 genes
* 某个细胞类型的 marker genes
* 某个疾病、衰老或处理相关的候选基因

接下来的一个自然问题是：

> **这些基因为什么会一起发生变化？它们是否受到某些共同转录因子调控？**

RcisTarget 就是用来回答这个问题的。

它通过分析候选基因调控区域中的 **transcription factor binding motif** 富集情况，寻找可能共同调控这些基因的转录因子 TF。

其核心过程可以概括为：

```text
Gene set
   ↓
Motif enrichment
   ↓
Enriched motifs
   ↓
Motif → TF annotation
   ↓
Candidate transcription factors
   ↓
Potential target genes
```

RcisTarget 官方定义的核心任务也是识别一个 gene set 或 genomic region set 中显著富集的 DNA motifs。

---

# 2. 从第一性原理理解 RcisTarget

## 2.1 转录因子如何调控基因？

转录因子 TF 是能够识别特定 DNA 序列的蛋白。

例如某个 TF 偏好结合：

```text
ACGTGCAA
```

这样的 DNA 序列模式称为：

**motif**

如果很多基因都受到同一个 TF 调控，那么这些基因附近的 promoter 或 enhancer 中，就可能反复出现这个 TF 对应的 motif。

因此，如果我们有一组基因：

```text
Gene A
Gene B
Gene C
Gene D
...
```

发现它们的调控区域中都富集某个 motif，就可以推测：

```text
这些 genes
      ↑
   某个 TF
```

可能存在共同调控关系。

这就是 RcisTarget 的基本出发点。

---

# 3. RcisTarget 最核心的设计：Ranking

RcisTarget 并不是拿到一个 gene set 后，临时扫描每个基因的 DNA 序列。

它提前建立好了一个巨大的 **motif ranking database**。

假设研究 motif M1。

数据库已经计算出：

```text
             motif M1

Gene A       rank 1
Gene F       rank 2
Gene K       rank 3
Gene B       rank 4
Gene Z       rank 5
...
Gene X       rank 20000
```

排名越靠前，说明这个基因附近的调控区域越支持 motif M1。

然后假设我们输入的 gene set 是：

```text
Gene A
Gene B
Gene F
Gene K
Gene M
```

其中：

```text
Gene A → rank 1
Gene F → rank 2
Gene K → rank 3
Gene B → rank 4
```

大量候选基因集中在 motif M1 ranking 的顶部。

RcisTarget 就会判断：

> motif M1 在这个 gene set 中明显富集。

反过来，如果候选基因随机散落在 ranking 中：

```text
rank 300
rank 2500
rank 8000
rank 15000
```

就没有明显 motif enrichment。

因此，RcisTarget 本质上是在问：

> **我的 gene set 是否异常集中在某个 motif ranking 的顶部？**

---

# 4. 为什么使用 Ranking，而不是直接看 motif 数量？

因为不同基因的调控区域长度、motif 数量和序列组成都不同。

直接比较：

```text
Gene A 有几个 motif
Gene B 有几个 motif
```

容易受到很多技术因素影响。

RcisTarget 将 motif information 转化为：

```text
Gene ranking
```

然后分析 gene set 在 ranking 中的位置。

因此整个问题从：

```text
这个 gene 有没有 motif？
```

变成：

```text
我的这组 genes 是否整体集中在这个 motif 的高排名区域？
```

这是一个典型的 **rank-based enrichment analysis**。

它和 GSEA 的思想非常接近。

GSEA 是：

```text
ranked genes
↓
pathway enrichment
```

RcisTarget 是：

```text
motif-specific ranked genes
↓
gene set enrichment
```

---

# 5. RcisTarget 的主要功能

RcisTarget 实际上完成三件事情。官方 `cisTarget()` workflow 也是由这三个步骤组成。

```text
Step 1
Motif enrichment
calcAUC()

        ↓

Step 2
Motif → TF annotation
addMotifAnnotation()

        ↓

Step 3
Identify enriched target genes
addSignificantGenes()
```

最终得到：

```text
Gene set
↓
Motif
↓
TF
↓
Target genes
```

---

# 6. 输入是什么？

RcisTarget 最重要的输入有三个。

## 6.1 Gene set

例如：

```r
geneLists <- list(
    aging_genes = c(
        "MT1",
        "BIRC3",
        "GPX4",
        "FOSL2"
    )
)
```

输入可以来自：

```text
DEGs
Module genes
Marker genes
Pathway genes
Disease-associated genes
Aging-associated genes
```

RcisTarget 支持一个或者多个 gene sets。

---

## 6.2 Motif ranking database

第二个输入是预先建立好的 motif ranking database：

```r
motifRankings <- importRankings(
    "xxx.genes_vs_motifs.rankings.feather"
)
```

这个数据库描述的是：

```text
每一个 motif
        ↓
所有 genes 的 ranking
```

因此可以想象成：

| Gene  | motif1 | motif2 | motif3 |
| ----- | -----: | -----: | -----: |
| GeneA |      1 |    532 |   3200 |
| GeneB |      4 |     13 |    640 |
| GeneC |     85 |      2 |    120 |

数字代表该 gene 在相应 motif 中的排名。

数据库必须与研究物种和所选择的调控区域匹配，例如 promoter 附近或者更大的 TSS-centered regulatory region。官方文档也明确要求根据 organism 和 search space 选择相应 ranking database。

因此：

**数据库选错，比参数调得不好更加严重。**

---

## 6.3 Motif annotation database

第三个输入用于完成：

```text
motif
↓
TF
```

因为一个 motif 本身只是一个 DNA sequence pattern，例如：

```text
motif_1234
```

需要进一步知道：

```text
motif_1234
↓
RELA
NFKB1
```

RcisTarget 的 motif annotation 包含直接注释、orthology 推断和 motif similarity 推断等不同证据等级。

---

# 7. 核心计算：AUC

对于每一个 motif，RcisTarget 都会检查：

> 输入 gene set 在该 motif ranking 顶部富集得有多明显？

然后计算 **AUC：Area Under the Curve**。

简单理解：

```text
ranking top
│
│  输入 gene
│  ↑ ↑ ↑ ↑
│
└──────────────── ranking
```

如果输入 genes 大量集中在前面：

```text
AUC 高
```

如果输入 genes 随机散布：

```text
AUC 低
```

因此：

**AUC 衡量某个 motif 对当前 gene set 的 enrichment strength。**

官方 `calcAUC()` 就是针对每个 gene set × motif 组合计算 AUC。

---

# 8. NES：RcisTarget 最重要的统计量

仅看 AUC 还不够。

RcisTarget 会把某一个 motif 的 AUC 与数据库中所有 motifs 的 AUC 分布进行比较：

$$
NES =
\frac{AUC_{motif}-mean(AUC)}
{SD(AUC)}
$$

这就是：

**Normalized Enrichment Score，NES**

所以 NES 可以简单理解为：

> 当前 motif 的 enrichment 比普通 motif 高出多少。

例如：

```text
motif A    NES = 1.2
motif B    NES = 2.1
motif C    NES = 5.6
```

motif C 的 enrichment 最强。

RcisTarget 默认通常使用：

```text
NES ≥ 3
```

作为显著 motif 的筛选标准。

因此实际看结果时：

**NES 通常是首先关注的指标。**

---

# 9. 如何找到真正贡献 enrichment 的 genes？

一个 gene set 可能有几百个 genes。

并不是所有 genes 都真正支持某个 motif。

例如：

```text
输入 200 genes

其中真正集中在 NF-κB motif ranking 顶部的
只有 35 genes
```

RcisTarget 可以进一步通过：

```r
addSignificantGenes()
```

寻找这些 highly ranked genes。

这些 genes 可以理解为：

> **真正推动该 motif enrichment 的核心 genes。**

它与 GSEA 中的：

```text
leading-edge genes
```

非常类似。

官方输出中的 `enrichedGenes` 就表示在对应 motif ranking 中高度富集的 genes。

---

# 10. RcisTarget 的最终输出

最终结果通常类似：

| motif   | NES |  AUC | TF_highConf | TF_lowConf | enrichedGenes        |
| ------- | --: | ---: | ----------- | ---------- | -------------------- |
| motif_1 | 6.2 | 0.19 | RELA        | NFKB1      | BIRC3;NFKBIA;TNFAIP3 |
| motif_2 | 5.1 | 0.17 | FOS         | JUNB       | FOSL2;JUN            |
| motif_3 | 4.4 | 0.14 | STAT1       | STAT2      | GBP2;IRF1            |

官方结果中最主要的信息包括 NES、AUC、高低置信度 TF annotation，以及显著富集的 genes。

因此最终可以得到一个调控关系：

```text
RELA
 ↓
BIRC3
NFKBIA
TNFAIP3
```

或者：

```text
STAT1
 ↓
GBP2
IRF1
...
```

---

# 11. 最重要的参数

真正需要关注的参数并不多。

| 参数               | 作用                         | 常用设置                   |
| ---------------- | -------------------------- | ---------------------- |
| `nesThreshold`   | motif enrichment 筛选阈值      | 默认 3                   |
| `aucMaxRank`     | 用 ranking 前多少 genes 计算 AUC | 常见约 1–10%              |
| `motifRankings`  | motif ranking database     | 根据物种和调控区域选择            |
| `motifAnnot`     | motif → TF annotation      | 与 ranking database 匹配  |
| `geneErnMethod`  | 识别 enriched genes 的方法      | `aprox` / `iCisTarget` |
| `geneErnMaxRank` | recovery curve 最大 ranking  | 常见默认 5000              |
| `nCores`         | 并行计算线程数                    | 根据计算资源设置               |

其中最值得关注的是：

```text
motifRankings
nesThreshold
aucMaxRank
```

`cisTarget()` wrapper 中 `nesThreshold=3`，`aucMaxRank` 默认使用 ranking 中约前 5% 的 genes；官方文档建议 AUC 搜索范围通常可在约 1–10% 内调整。不同版本或直接调用 `calcAUC()` 时，应以当前安装版本的帮助文档为准。

---

# 12. geneErnMethod 怎么选？

识别 enriched genes 时主要有：

```text
aprox
iCisTarget
```

`aprox`：

```text
速度快
适合多个 gene sets
```

`iCisTarget`：

```text
计算更严格
计算量更大
```

官方文档将 `iCisTarget` 描述为更接近原始 iRegulon/i-cisTarget 方法，而 `aprox` 是速度更快的近似实现。

一般批量分析：

```r
geneErnMethod = "aprox"
```

通常已经足够。

---

# 13. 一个最简单的 RcisTarget workflow

核心流程实际上只有几步：

```r
library(RcisTarget)

# 1. 输入 gene set
geneLists <- list(
    target_genes = c("BIRC3", "GPX4", "FOSL2", "GBP2")
)

# 2. motif ranking database
motifRankings <- importRankings(
    "genes_vs_motifs.rankings.feather"
)

# 3. motif enrichment
motifs_AUC <- calcAUC(
    geneLists,
    motifRankings
)

# 4. motif → TF
motifEnrichmentTable <- addMotifAnnotation(
    motifs_AUC,
    nesThreshold = 3,
    motifAnnot = motifAnnotation
)

# 5. 找到真正贡献 enrichment 的 genes
results <- addSignificantGenes(
    motifEnrichmentTable,
    geneSets = geneLists,
    rankings = motifRankings,
    method = "aprox"
)
```

也可以直接使用 `cisTarget()` 完成完整 workflow。官方 wrapper 本身就是依次执行 motif enrichment、TF annotation 和 significant gene identification。

---

# 14. 结果应该怎么解析？

假设得到：

```text
motif         NES     TF        enrichedGenes

motif_A       6.1     RELA      BIRC3,NFKBIA,TNFAIP3
motif_B       5.4     NFKB1     BIRC3,NFKBIA
motif_C       4.8     STAT1     GBP2,IRF1
```

可以分三层理解。

第一层：

```text
NES = 6.1
```

说明 motif_A 在你的 gene set 中明显富集。

第二层：

```text
motif_A → RELA
```

说明这个 motif 被注释为 RELA 可能识别的 motif。

第三层：

```text
BIRC3
NFKBIA
TNFAIP3
```

是最主要贡献该 motif enrichment 的候选 target genes。

因此可以提出：

```text
RELA
 ↓
BIRC3 / NFKBIA / TNFAIP3
```

可能是这个 gene set 背后的一个共同转录调控程序。

---

# 15. 一个非常重要的注意点：RcisTarget 不能证明 TF 活性

RcisTarget 找到：

```text
RELA motif enrichment
```

代表的是：

> 这些 genes 的调控区域具有 RELA 相关 motif 的富集证据。

它并不能直接证明：

```text
RELA 当前一定被激活
```

也不能直接证明：

```text
RELA 一定真实结合了这些 genes
```

因为 motif enrichment 是一种 **regulatory potential inference**。

更加完整的证据链通常需要结合：

```text
RcisTarget
motif evidence
        +
TF expression/activity
        +
ChIP-seq / CUT&Tag
        +
gene expression
```

例如：

```text
RcisTarget
RELA motif enriched
        ↓
DoRothEA / decoupleR
RELA activity increased
        ↓
BIRC3 expression increased
```

这样的证据明显强于单独使用 RcisTarget。

---

# 16. RcisTarget 和常见 TF 分析方法的区别

可以把几个常见工具放在一起理解：

```text
RcisTarget
gene set → motif → TF

DoRothEA
target gene expression → TF activity

SCENIC
co-expression
     +
RcisTarget motif enrichment
     ↓
TF regulon

ChIP-seq / CUT&Tag
实验测量 TF 实际 DNA binding
```

所以 RcisTarget 最适合回答的问题是：

> **这一组具有共同变化的 genes，背后可能受到哪些共同 TF 调控？**

---

# 17. 一句话理解 RcisTarget

RcisTarget 的本质可以概括为：

```text
给我一组 genes
        ↓
看它们是否集中在某些 motif ranking 顶部
        ↓
找到 enriched motifs
        ↓
motif 映射到 TF
        ↓
推断潜在的上游 transcriptional regulators
```

因此它本质上是一个：

**基于 motif ranking 的转录因子调控富集分析工具。**

它把普通的：

```text
“这组 genes 有什么功能？”
```

进一步推进到：

```text
“这组 genes 为什么会一起变化，
可能是谁在调控它们？”
```

这就是 RcisTarget 最主要的价值。
