# 模型开发知识库（索引）

本目录把 `literature/` 中每篇论文能给模型提供的**物理图景、方程、参数、实验数据、局限**整理成可查的条目。
每个数值都标明来源与性质：**测量 / 文献引用 / 模型假设 / 反标定 / 本项目数值验证**。
写模型、选参数、做验证前先查这里，再回原文核对页码或式号。

## 文件

| 文件 | 内容 |
|---|---|
| `porous_salt_parameters.md` | **跨论文的结晶盐多孔介质参数汇编**（ε_s、孔径、k_s、P_cth、a_v、d_c/d_sub、D_v,eff、a_w），含冲突与取值建议 |
| `eloukabi2013.md` | Eloukabi, Sghaier, Ben Nasrallah, Prat 2013, IJHMT 56:80 —— drying；crusty–patchy 转变；毛细柱模型式(12)(13) |
| `veran_tissoires_prat2014.md` | Veran-Tissoires & Prat 2014, JFM 749:701 —— **wicking**；成核位置；§6.2 格子生长模型；§6.3 Kelvin 与稀疏开孔 |
| `shokri_kuehni2022.md` | Shokri-Kuehni, Sahimi, Shokri 2022, Drying Technol. 40:691 —— perspective；开放问题；盐堆积模拟建议 |
| `ben_amara2025.md` | Ben Amara, Sghaier, Kharaghani 2025, Drying Technol. 43:1992 —— drying；粗层覆盖细基底；结构光测高；**主要验证数据** |
| `wijnhorst_prat_shahidzadeh2025.md` | Wijnhorst, Prat, Shahidzadeh 2025, arXiv 2508.18779 —— 盐结构分级多孔结构；分形维数；**毛细-黏性不限高**的分级模型 |
| `prat2026.md` | Prat 2026, Langmuir —— 局部生长律 V_n = G·F(Da*)·j；代码实现见 `../prat2026/` |
| `wijnhorst2024.md` | Wijnhorst, Van der Sloot, Pel, Shahidzadeh 2024, PRApplied 21:064055 —— drying 截锥砂岩；Pe 1.5 给非阻塞薄膜、Pe 2.7 给 cauliflower 合并成阻塞壳；**顶部失连后结晶转到壳底、盐腿抬壳**；表面晶体湿/干循环 |
| `rodriguez_rodriguez2026.md` | Rodríguez-Rodríguez et al. 2026, arXiv 2607.11587 —— 毛细管单向蒸发的盐爬行；**φ_eff/φ_sat ≈ 1.20 的异质成核阈值与潮解滞后**；吸液回湿间歇；爬行型沉积致密细晶、自由生长疏松粗晶 |
| `nucleation_supersaturation.md` | 专题页（非单篇）：首次结晶的过饱和（PRL 2012 的 1 到 Naillon 2017 的 1.72）、成核之后浓度如何回到饱和（Naillon 2015/2017/2018、Hidri 2013、Derluyn 2014、DuMuX 各模型）、k_R 的单位陷阱、400 µm 孔的弛豫时间估算；**由检索子代理整理，未逐条复核** |
| `piotrowski2020.md` | Piotrowski, Huisman, Nachshon, Pohlmeier, Vereecken 2020, Geosciences 10:423 —— drying；**分离盐壳的气体渗透率直接测量**，NaCl 4.77e−12 m²（240 µm 砂） |

## 每篇论文属于哪个 regime、支撑哪个模块

| 论文 | regime | 基底 | RH / T | 支撑的模块 |
|---|---|---|---|---|
| Eloukabi 2013 | drying | 玻璃珠 5–300 µm，6 档 | 3 % / 30 °C | 供液上限判据；crusty–patchy 随粒径与 j 的转变；patch 增强蒸发 |
| VT & Prat 2014 | **wicking** | 玻璃珠（1 mm 为主）+ 5–50 µm 结壳实验 | 7 % / 室温 | 成核位置（局部 Pe）；屏蔽；格子生长算法；干涸门控；结壳的蒸汽阻力解释 |
| Shokri-Kuehni 2022 | perspective | 砂 | — | 问题清单；盐层湿润证据；盐堆积模拟思路 |
| Ben Amara 2025 | drying | 玻璃珠 1–50 与 400–600 µm，分层 | 7 % / 21 °C | **验证数据**：质量损失、高度演化、覆盖形貌、粗层厚度幂律 |
| Wijnhorst 2025 | drying（砂岩）+ 烧杯 | 砂岩孔 30 µm / 无基底 | 45 % / 21 °C | 盐结构内部微结构与 (k_s, P_c) 的解耦；分级模型；生长停止=盐量耗尽 |
| Prat 2026 | 局部闭合，regime 无关 | 砂 0.58 mm；珠 90 µm | — | j 到 V_n 的定量律；Da* 与基底粒径的关系 |
| Piotrowski 2020 | drying | F32 石英砂 240 µm，k=2.1e−11 | 20–40 % / 25 °C | patchy 型 k_s 的唯一直接测量；壳/砂渗透率比 0.23；"crust"一词跨 regime 混用的警示 |
| Wijnhorst 2024 | drying | Prague 砂岩，孔 30 µm，截锥 | 45 % / 21 °C | Pe 决定薄膜 vs 合并壳；失连后壳底结晶与盐腿抬壳（脱离状态的机制）；湿/干循环 |
| Rodríguez-Rodríguez 2026 | 毛细管，无基底 | 玻璃管 | 10–70 % | 成核阈值 S ≈ 1.2；回湿滞后；生长模式→微结构 |

## 跨论文一致的结论（可直接作为模型假设）

1. **外部盐是多孔介质且在生长期保持湿润**：2013、2014、2022（热成像冷点）、2025 Ben Amara（扫描时"大部分盐霜因毛细泵吸保持湿润"）、2025 Wijnhorst（micro-CT）、2026 一致。
2. **盐结构孔径小于基底孔径**：2014、2025 Wijnhorst（d_c/d_sub≈0.14，孔≈晶粒/4）、2026 一致；结壳情况下小到纳米（2014 §6.3 推断，Licsandru 2022 实测）。
3. **粗基底给 patchy、细基底给 crust**，转变在平均粒径 100–200 µm 附近（2013 Fig 8：≥200 µm 回收率≈100 %；Ben Amara：1–50 结壳、400–600 patchy）。
4. **patchy 不会铺满表面**的机制是屏蔽 + 离子改道（2014 §5、Ben Amara §3.4 复述）。
5. **crust 使蒸发降约 14 倍但不完全阻断**（2014 Fig 23，5–50 µm 珠；Ben Amara 引用同一数字）。
6. **patchy 可使蒸发率高于纯水**（2013 Fig 2f，300 µm 珠、C0=25 %），因盐结构外表面扩大交换面积，抵消 a_w≈0.75 的降幅。

## 论文之间的张力（模型必须显式处理）

| 张力 | 一方 | 另一方 | 对模型的含义 |
|---|---|---|---|
| 毛细-黏性竞争是否限制高度 | 2013 式(12)、2014 §6.2：黏性阻力限高并导致结壳 | Wijnhorst 2025 分级模型：P_c(h_e)/P_c,max < 1，**不限高**，高度由盐量决定；Prat 2026 讨论：j 变小时压降也变小，质疑该机制解释生长减缓 | 供液门控必须用**解耦的 (k_s, P_cth)**：分级结构里粗内孔给 k_s、细外孔给 P_cth，两者不能由同一孔径经 Carman–Kozeny 推出。单尺度细孔（crust）才会触发黏性限高 |
| 盐结构孔隙率 | 测量：0.43（Dai）、0.40、0.25（Licsandru 结壳）、≈0.20（Wijnhorst 分级） | 模型假设：0.04（2013）、0.10（2014 P_cth 用） | 2013/2014 的假设值只服务于定性演示，不能当参数用 |
| 生长减缓/停止的原因 | 基底干涸耗尽盐量（Wijnhorst 2025、Prat 2026 首选） | 盐霜自身干涸（旧文献）；黏性限高（2013/2014）；盐霜与基底间形成缝隙（Licsandru 2022：慢蒸发期缝隙增大但盐霜仍饱和） | wicking 模型里"停止"不能靠盐量耗尽，需要显式的门控机制；缝隙形成是尚无模型的候选 |
| RH 对形貌的作用 | Wijnhorst：RH 20 % 出薄膜结壳、RH 45 % 出 cauliflower | Ben Amara / 2014：RH 7 % 下粗珠仍 patchy | 形貌由 (孔径, j) 共同决定，2013 §4 已预言"同一粒径改 j 也能转变"；模型的判据应是 Π 这类含 j 的无量纲数，不是孔径单独 |
| C_sat 的写法 | 质量分数 0.264（Prat 2026、Ben Amara）、0.266 @30 °C（2013） | 0.36（Wijnhorst，是 g 盐/g 水） | 换算时勿混用 |

## 尚无任何论文闭合的环节（模型开发的真正空白）

- d_sub 到 (k_s, P_cth) 的定量映射（2014 明说"所知甚少"，2026 仍说"stereological parameters 仍是研究课题"）。Piotrowski 2020 给了 240 µm 砂一个点，纳米孔 crust 仍无直接测量。
- crust 的蒸汽有效扩散系数与 a_w 随厚度/干燥程度的演化（Ben Amara 只说"变干变密后更不透"；Shokri 2022 指向 Piotrowski 2020 的盐壳气体渗透率测量与 Nachshon 2018 的三维壳结构）。
- 成核判据 C_nuc（2014 认为需过饱和；2013 取 C_sat 即结晶）。
- 分级结构里 d_c(z) 的演化律（Prat 2026 "Remaining topic"）。
- 盐霜与基底之间缝隙的形成机制。

## 本项目已有的数值事实（来自 `../salt_grow`，见各文件"本项目验证"节）

- Prat 2026 生长律实现验证到 0.08 %；ρ 在形貌计算中约掉。
- H* = ρ k_s P_cth/(μ j̄) 是格子模型的控制群；Carman–Kozeny 闭合被 Licsandru 实测证伪（高估 k 三个数量级）。
- 纯"生长律 + Laplace 场"不产生结壳（覆盖率恒定），干涸门控是唯一正则化。
- 均匀 F 不改形貌；纯 Laplacian 生长无内禀尺度，形貌不随网格收敛。
- 三维自屏蔽远弱于二维（E 1.10 对 2.04）。
