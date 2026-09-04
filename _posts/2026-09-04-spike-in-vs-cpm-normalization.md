---
layout: post
title: "Spike-in normalization 与 CPM normalization"
date: 2026-09-04
description: "比较 ChIP/CUT&Tag 分析中 Spike-in normalization 与 CPM normalization 的原理和适用场景"
tags: [技术, ChIP, CUT&Tag, Spike-in, CPM, normalization]
---

## Spike-in normalization

不是把每个样本的**总 reads**变成最小 spike-in 数，而是按照 spike-in reads 的比例缩放，使每个样本在校正后都相当于具有相同数量的 spike-in：

`SF_i = min(spike-in reads) / sample spike-in reads`

校正后：

`spike-in reads_i × SF_i = min(spike-in reads)`

实验物种 reads 也按照相同比例缩放，但校正后实验物种的总 reads 不一定相同。这个差异可以保留真实的全局 ChIP/CUT&Tag 信号变化。

## CPM normalization

CPM 是把每个样本的实验物种总 reads 统一换算为 100 万：

`CPM = region coverage / experimental mapped reads × 10^6`

校正后所有样本的实验物种总信号理论上都是 100 万，因此会消除样本之间的全局信号总量差异。

## 直观区别

- **CPM**：以样本自身的实验物种 reads 为尺子，假设所有样本总体信号相同。
- **Spike-in**：以加入的外源 reads 为尺子，允许实验物种总体信号发生真实变化。

所以更准确地说：

> Spike-in normalization 把所有样本缩放到相同的“外源参考量”；CPM 把所有样本缩放到相同的“实验物种总 reads = 100 万”。
