# 成核过饱和与成核之后的浓度：别人怎么做（2026-09-19）

**来源与可信度。** 本页由一个检索子代理读开放获取全文（CORE / arXiv / OATAO）和本地文件整理，主会话**没有逐条对照 PDF 复核**；
要在代码里用到其中任何方程或数值之前，先按 CLAUDE.md 的规则回到原文核对。标"未核"的条目连子代理也没有读到原文。
本地可核对的：`../salt_grow/paper.txt`（JFM 2014）、`literature/PhysRevLett.108.054502 (1).pdf`、`../salt_grow/reference/prat2026.txt`。

体系：以下实验多为毛细管 / 微流控通道里的**干燥**；用到蒸发-芯吸模型时只借"成核与晶体生长动力学"这一局部物理。

## 1. 首次结晶需要多大的过饱和

| 来源 | S = C/C_sat（或 m/m_sat） | 备注 |
|---|---|---|
| Veran-Tissoires, Marcoux, Prat 2012, PRL 108, 054502 | **1**："Crystallization occurs when the saturation concentration C_sat is reached (supersaturation effects are negligible for NaCl)" | 不推断过饱和 |
| Veran-Tissoires & Prat 2014, JFM 749, §4.4 Table 3 | 1.05–1.26 | 由首次析出时间反推 |
| Hidri, Sghaier, Eloukabi, Prat, Ben Nasrallah 2013, Phys. Fluids 25, 127101 | ≈ 2（拟合），自称"more qualitative than quantitative" | 热点处；表面平均仍低于 C_sat |
| Naillon, Duru, Marcoux, Prat 2015, J. Cryst. Growth 422, 52 | 全局 1.39 ± 0.1，局部 1.57 ± 0.1 | 1 mm 玻璃管，22 °C，RH 7 % |
| Naillon, Joseph, Prat 2017, J. Cryst. Growth 463, 201 | 1.72（99 次，1.48–2.12） | 5×5、20×20 µm² PDMS 通道 |
| Desarnaud, Derluyn, Carmeliet, Bonn, Shahidzadeh 2014, JPCL 5, 890 | 1.6 ± 0.2（> 100 次），与毛细管尺寸（50–1000 µm）、形状、润湿性无关 | 低于 1.6 可"保持数天" |
| Rodríguez-Rodríguez 2026（本库另页） | 1.20 | 玻璃壁上的异质成核 |

## 2. 成核之后浓度怎么回到饱和

一致的图像（Hidri 2013 §IX；Naillon 2015 / 2017 / 2018；Derluyn & Prat 编 2024 书第 3 章 §3.4）：过饱和以**有限速率**被消耗，分两段——
先是表面反应控制（10–100 ms），之后晶体表面处于 c_eq，溶液靠离子扩散（加生长诱导的对流）向这个 Dirichlet 汇弛豫。
书中原话："once the crystallization has started locally, the excess ions are very quickly consumed. This question is still not resolved."

- **速率律**（Naillon 2017 式 (1)(2)(8)）：J_R = k_R (c_i − c_eq)，一阶；**c 是质量浓度 [kg/m³]**，k_R [m/s]；晶面速度 w_cr = (k_R/ρ_cr)(c_i − c_eq)。
  与扩散阻力串联：k_G = k_D k_R/(k_D + k_R)。测得 k_G = 2.33e-3 m/s（10 ms 平均），作为 **k_R 的下界**报告；
  此前文献的 1e-5–1.3e-4 m/s 被证明都是扩散控制的 k_G（k_G = 2√(D/π t_m)）。
- **陷阱**：Prat 2026 把它写成 w_cr = k_r (C_i − C_eq)，C 是**质量分数**，与 2017 的定义差一个 ρ_l/ρ_cr ≈ 0.55。`prat2026/params.py` 的 `KR_T21` 沿用 Prat 2026 的写法。
- **成核之后 Prat 组自己用的闭合是平衡闭合**：Naillon 2015 式 (17)（低 Pe 极限，"precipitates immediately … x_s ≈ x_sat"）；
  JFM 2014 式 (6.5) j_cr = j·C_sat。PRL 2012 与 JFM 2014 的孔网络都只算到第一颗晶体为止；Prat 组**从未**把成核后的动力学律放进 PNM。
- 带动力学律的其他模型：Jambhekar 2015（DuMuX，q = k_p φ ρ S_l (x − x_max)，k_p = 1 s⁻¹，"equilibrium"，无阈值）；
  Jambhekar 2016（r = k_n A_n |1 − Ω^θ|^η，无阈值）；Schollenberger 2025（两相 PNM，s = k A_react (Ω − 1)，k = 1e-3 mol/(m² s)，无阈值，与本机 DuMuX `pnm2pncmin2` 相同）；
  **Derluyn, Moonen, Carmeliet 2014, JMPS 63, 242**：唯一同时带"阈值 > 1"和动力学律的——U_thr = 1 + (U_start − 1) exp(−υ S̄_cr)，U_start = 1.5 或 2，
  **阈值随附近晶体积累弛豫回 1**（成核要过饱和、已有晶体上生长不要）。Espinosa 2008 的形式与常数仅经 Derluyn 2014 转述，未核。

## 3. 弛豫时间的量级（400 µm 立方孔）

两个估算，几何假设不同，都是估算：

- 子代理：**晶体在孔内**，50 µm 种子，准稳态扩散到球形汇 τ_D = V/(4π D R_eff)，晶体长大时积分：半衰 40 s，e 折 ≈ 54 s，90 % ≈ 105 s；
  Da = k_R R_eff/D = 50–140，**扩散控制**，k_R 取无穷大结果变化 < 1 %。
- 主会话（README E1 §15.9）：**晶体在孔口之外**（外部盐霜），孔里的过饱和靠向上平流冲刷（ρV/ṁ_b ≈ 25 min）和经被盖住的开口面积扩散（V r/(D A_盖)：种子 55 min → 全盖住 1 min）交给结构：几十分钟。

两者都说明"一步重置"在形式上不对；对小时量级的输运过程，它在效果上是否可接受取决于晶体在哪里。

## 4. 本项目的取法

2026-09-19 起 `S_NUC = 1.0`（用户决定）：不设过饱和，到 C_sat 即出现第一颗晶体。先例：PRL 2012、DuMuX 三个模型。
若以后要恢复 S_nuc > 1：发表过的做法是 Derluyn 式的阈值滞后 + 一阶弛豫 dC/dt = −(C − C_sat)/τ（τ 取扩散时间，不要用小 k 的表面反应律——那会把输运阻力算两遍）。
