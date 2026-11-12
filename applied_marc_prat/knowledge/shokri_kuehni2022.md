# Shokri-Kuehni, Sahimi, Shokri (2022) —— 盐水蒸发预测的 perspective

*A personal perspective on prediction of saline water evaporation from porous media.*
Drying Technology 40(4), 691–696. 文件：`literature/A personal perspective on prediction of saline water evaporation from porous media.pdf`

**性质：短评/展望，无新数据、无方程。**价值在于问题清单与文献指针。

## 论文提出的核心缺口

> 目前没有任何定量模型能够在不引入调节参数的情况下预测"演化中的多孔结晶盐"对蒸发通量的影响。

这正是我们模型的目标定位，可作为 Introduction 的引用锚点。

## 可引用的实验证据

- **结晶盐多孔且湿润**：SEM（Fig 1，砂柱、3.5 M、5 天）；热成像冷点（Fig 2，20 % NaCl、水位在表面下 4 cm 的**wicking 型**实验，
  Shokri-Kuehni et al. 2017 GRL）说明蒸发优先经盐层发生。
- **总蒸发对浓度不敏感**：0.5 / 1.5 / 2.5 M 砂柱累计蒸发相近（Fig 3），归因于表面盐层参与蒸发，抵消 a_w 降低。
- 盐层可形成枝状/树突结构，显著增大蒸发面积。

## 作者建议的方向

- 分子动力学研究成核（受限孔、双水模型）；
- **立方晶体堆积的计算模拟**给出盐层的 ε、k（Fig 4，ref 31/32）—— 这与 Prat 2026 用立方体堆积拟合 a_v(ε_s) 的思路一致，
  也是我们 (k_s, P_cth) 闭合缺口的可能来源；
- 结果向 PNM 与连续尺度上传。

## 有用的文献指针（PDF 未收录）

| 引用 | 内容 | 对我们的用处 |
|---|---|---|
| Piotrowski et al. 2020, Geosciences 10:423 | 蒸发形成盐壳的**气体渗透率**测量 | crust 蒸汽阻力闭合的实验依据 |
| Nachshon et al. 2018, GRL 45:6100 | NaCl 壳的三维结构与对蒸发的影响 | crust 微结构 |
| Dashtian, Shokri, Sahimi 2018, AWR 112:59 | 蒸发致盐析的 PNM（含相关性与非均质） | PNM 耦合参照 |
| Jambhekar et al. 2016, TiPM 114:341 | 自由流–多孔介质耦合的动力学沉淀模型（DuMuX） | 与本组 DuMuX 工作衔接 |
| Rufai & Crawshaw 2018 | 润湿性变化与盐沉积导致的渗透率损伤 | 盐层 k 演化 |
| Shokri-Kuehni et al. 2017 AWR 105:154 | 盐种类与环境条件（RH、T）的影响 | RH 张力的补充数据 |
