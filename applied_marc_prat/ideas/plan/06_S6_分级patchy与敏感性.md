# S6 分级 patchy 分支与敏感性矩阵

**Regime：evaporation-wicking。** 分级规则来自 Wijnhorst 2025 的 drying 砂岩与烧杯实验，假设微结构不依赖 regime（[A]）；Piotrowski 的 k_s 实测也是 drying 样品。

## 1. 目的

给 `PATCHY_POROUS_SALT` 加上按生成高度分级的 birth 规则，使"粗内孔给 k_s、细外孔给 P_cth"自动成立；
然后系统扫描两个分支的参数，把 I§7 的"哪些参数真正影响结果"从预期变成数值事实。

## 2. 依据

- 物理：Wijnhorst 2025，相邻层晶粒与孔径减半、孔 ≈ 晶粒/4、ε_s ≈ 0.20；分级模型式(1)（P_c = 2γcosθ/a_n）与式(3)（逐层压降求和）。
  Prat 2026 自复制 regime（r = 1，ε_s = 0.43）作为另一极端。d_c,0 = 0.14 d_sub（M，Wijnhorst 2024）。
- 参数（I§3.1）：d_c(n) = d_c,0 r^n；a_1(n) = d_c(n)/4；k_s(n) = c_k a_1(n)³/(12 d_c,0)（A，立方律）或 CK；P_cth(n) = 2γcosθ/a_1(n)；r_p(n) = a_1(n)/2；
  层序 n(z) 由 h_1 = 0.1 mm、α = 2（A，Wijnhorst 算例）；减半下限 d_c ≥ 1 µm（A）。
  前因子 c_k 是 B：锚 Piotrowski 2020 的 240 µm 砂 4.77e−12 m²（M，气测，patchy 型）。
- 门控分母在 patchy 里由重力主导（S4），预期 k_s、P_cth 对 patchy 形貌不敏感；这是本步要证实或推翻的。

## 3. 前置

S5 全部通过。

## 4. 实现

- `constitutive.py`：`PatchyPorousSalt(d_sub, r, eps_s, c_k, h_1, alpha, d_c_min)`，`birth(z)` 返回 (eps_s, d_c, k_s, P_cth, r_p, gen)；每个 `SaltConstitutive` 同时给出默认初始态（A 修订）：patchy 稀疏种子、crust 整层膜。
- `advance.py`：翻转时调用 `birth(z_voxel)` 写入体素属性数组；`darcy.py`/`gating.py` 改读体素属性场（S4 里是均匀标量）。
- `ledger.py`：结构级有效量 k_eff（对整个 `SALT_WET` 域施加单位压差解一次 Darcy）、有效毛细阈值（界面体素 P_cth 最小值）。
- `runs/S6_sensitivity/`：参数矩阵与结果表。

## 5. 验收（`tests/test_S6_hierarchy.py`）

| 测试 | 判据 | 依据 |
|---|---|---|
| 一维分级柱 | n 层串联压降 Σ μ U h_i/k_i 与 Wijnhorst 式(3) 解析一致 <1e−10 | 解析 |
| Π 的层分布 | 界面层 Π 用最细层 P_cth、压降由最粗层 k_s 主导；n = 6 时 Π ≫ 1 | W25 Fig 5c |
| c_k 标定 | d_sub = 240 µm 的 k_s(0) 落在 4.77e−12 的 3× 内，c_k 记录 | Piotrowski 2020 |
| k_eff | 均匀 k_s 时 k_eff = k_s；分级时 k_eff 介于最粗层与最细层之间且由粗层主导 | 定义 |
| 敏感性矩阵 | 见下表；每项给出 h̄、η∞、E、f_dry 对参数的对数灵敏度 | I§7 |

敏感性矩阵（每格一次全耦合运行，δ = 15 mm，RH 7 %）：

| 分支 | 强影响（预期） | 弱影响（预期） | 判据 |
|---|---|---|---|
| patchy | ε_s（0.20–0.43）、d_sub 经 Da*（25/90/500 µm）、δ（2.5/7.5/15 mm）、种子密度（S3 已扫，此处只取两档） | k_s（×1e−3–×1e3）、P_cth（×0.1–×10）、r、h_1、α | 弱影响项对 h̄、η∞、E 的变化 <5 % 则记"不敏感"；否则升级为标定量 |
| crust | k_s·P_cth 乘积（×0.5–×2）、φ_o（×0.5–×2）、ε_s（0.15–0.30）、初始态（整层膜 / 密集种子） | k_s 与 P_cth 各自拆分、r_p（>10 nm）、ε_s,min（0.02–0.10） | 同上 |

## 6. 完成后

可以声称：patchy 形貌与速率对 (k_s, P_cth) 不敏感、对 ε_s 与 Da* 敏感（若证实）；crust 只受乘积与开度控制；
因此模型的自由度只有 c_k、k_s·P_cth、φ_o 三个，其余参数固定不影响结论。
不能声称：分级律 h_1、α、r 的取值正确（只能说结果对其不敏感）；d_sub → 微结构切换。

## 7. 陷阱

- 层序按体素生成高度 z 而不是按生成时间；侧向铺展的低处体素永远是第 0 层。
- 若 patchy 对 P_cth 竟敏感，先查重力项是否被 L_s 放大（基底 Jurin 高度可能先于盐结构触限，2014 估 1 mm 珠 4.5 cm）。

## 8. 工作量

2 人日。
