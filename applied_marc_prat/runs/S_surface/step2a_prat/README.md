# Step 2A：Prat (2026) 单尺度生长律的独立验证（一维反应–输运参考解）

分支 `verify/step2a-prat-single-scale`，父提交 ea48aac（`verify/step1-grid-refinement-10um`）。
`run.py`（约 3 s）→ `results/step2a_results.json`、`step2a_table.csv`、`step2a_summary.md`；`plot.py` → `figures/`；测试 `tests/test_step2a_prat.py`。
未耦合气相求解器或 PNM；未改任何生产内核。

## 1. 物理模型与方程

已存在的、均质、全湿的多孔 NaCl 结构在给定蒸发通量与饱和盐水供给下一维竖直生长。Prat (2026) 式号：

    (1)  ∂(ρεC)/∂t + ∂(ρCU)/∂z = ∂(ρ ε D_eff ∂C/∂z)/∂z − a_v ρ k_r (C − C_eq)
    (2)  −ρ_cr ∂ε/∂t = a_v ρ k_r (C − C_eq)
    (3)  ∂(ρε)/∂t + ∂(ρU)/∂z = −a_v ρ k_r (C − C_eq),  ρ 常数 (4)
    (6)  ρU = [ρε + ρ_cr(1−ε)] dh/dt + j        z = h
    (7)  ρCU − ρ ε D_eff ∂C/∂z = [ρCε + ρ_cr(1−ε)] dh/dt     z = h
    (20) dh/dt = β (ρ/ρ_cr) k_r (C_top − C_eq);   C(0) = C_eq

解析生长律（式 29/32）：dh/dt = [C_eq/(1−C_eq)]·[j_w/(ρ_cr(1−ε_new))]·F(Da*)，F = √Da*/(1+√Da*)，Da* = β²(1−ε)²(1−C_eq)k_r/(ε D_eff a_v)。
该形式所需的假设：(i) Da = h²(1−C_eq)a_v k_r/(ε D_eff) ≫ 1，使 tanh√Da → 1（式 28 → 29）；(ii) 供给饱和，C(0) = C_eq，盐输送率恰为 C_eq j/(1−C_eq)；(iii) Pe ≪ Da，带内对流项可略；(iv) 带内 ε、a_v、D_eff 取出生值（线性化），但 (2)–(3) 的孔隙率下降保留在液相平衡中——这正是 (1−C_eq) 因子的来源；(v) 稳态行波，C_top − C_eq ≪ C_eq；(vi) 无毛细供给限制、等温。

## 2. 参数（Prat 2026 复现约定，≈21 °C）

| 参数 | 值 | 标签 |
|---|---|---|
| C_eq | 0.264 | [L] 论文 Materials and Methods |
| ρ_cr | 2170 kg/m³ | [L] 式 2 下 |
| ρ（饱和盐水） | 1200 kg/m³ | [B] 反标定（`prat2026/params.py`） |
| ε_new | 0.43 | [A] 出生孔隙率假设（论文引 Dai 的实测效应孔隙率；不声称终态等于它） |
| k_r | 2.3e−3 m/s | [L] Naillon |
| D_s | 1.3e−9 m²/s | [L] 论文正文（不是 PNM 的 1.54e−9） |
| β | 2 | [L] |
| j_w | 1.6e−5 kg m⁻² s⁻¹ → U_evap = 1.3333e−8 m/s | [A] 给定 |
| d_sub | 112.5 µm | [A] 名义值，非由表面孔半径导出 |
| d_c = 0.14 d_sub | 15.75 µm | [L, 经验比] 非验证的局部孔–晶映射 |
| h0 | 0.3 mm | [A] |

脚本断言这些值与 `prat2026/params.py` 及参考解 `PratParams` 的默认值一致。

## 3. 方法

- **Task A**：`prat2026/core.py`（经 `saltv2.growth.prat`）。`core.specific_surface` 用论文取整的系数 36：a_v = 36(1−ε)/d_sub；未取整为 5(1−ε)/d_c = (5/0.14)(1−ε)/d_sub。两者都报告；参考解用未取整式（式 16）。
- **Task B**：`layered_reference.py` 是 `review/prat-layered-validation` @ d80e791 的 `reviews/prat_layered_validation/layered_reference.py` 的逐字复制（sha256 `5f001da3…0ee75`，按路径加载，不是 `saltv2` 模块）。变体 V0：`evolve_eps=True, liquid_consumption=True, local_props=False, transport_eps_fixed=True`，即 Prat 自己的线性化。界面随动坐标，有限体积，几何拉伸网格（n_band 格/δ_C），后向 Euler，Picard 到 1e−10。dh/dt 来自式 (20) 与数值 C_top，**不**用式 29。稳态：`homogeneous_steady`，n_band 20/40/80；瞬态：`ReferenceSolver(h0 = 0.3 mm)`，24 h，dt = 0.5 δ_C/v。
- **Task C**：同一稳态解的各项收支。

## 4. 结果

**A. 解析（未取整 a_v；括号内为系数 36）**：a_v = 1.8095e5 (1.8240e5) 1/m；D_eff = 8.5247e−10 m²/s；δ_C = 1.0939 (1.0896) µm；Da* = 33.167 (32.904)；F = 0.85205 (0.85155)；C*_top = 1.5542e−6；v = 3.95346e−9 m/s = 14.2324 µm/h，v/U_evap = 0.29651；盐输送 5.7391e−6 kg m⁻² s⁻¹ → 前沿 4.8900e−6 (F)、内部 8.4910e−7 (1−F)；ε_final 闭合 0.33103。36 对 5/0.14：a_v +0.80 %，Da* −0.79 %，F 与 v **−0.059 %**。大 Da 检验：h0 处 Da = 7.58e4，Pe = 1.09e−2，tanh√Da = 1 到机器精度，式 28 与 29 差 2e−16。参考轨迹 h(t) = 0.3 mm + 14.2324 µm/h · t。

**B. 独立参考解（V0）**

| n_band | N | v/U_evap | 对式 29 | F（v/v_max） | F_num = P_i/(P_i+S) | C*_top 对式 31 | ε_final | 1−(1−ε)/F_num | 盐残差 | 固体残差 |
|---|---|---|---|---|---|---|---|---|---|---|
| 20 | 123 | 0.2965340 | +8.4e−5 | 0.852122 | 0.852122 | | 0.331082 | 0.331082 | 1.3e−10 | 4.7e−11 |
| 40 | 146 | 0.2965201 | +3.7e−5 | 0.852082 | 0.852082 | +3.7e−5 | 0.331050 | 0.331050 | 4.4e−10 | 1.9e−10 |
| 80 | 169 | 0.2965156 | +2.1e−5 | 0.852069 | 0.852069 | | 0.331040 | 0.331040 | 5.1e−10 | 3.2e−11 |

剖面对式 19 指数层的最大偏差 /C*_top：2.2e−4 / 3.9e−5 / 1.4e−5；Richardson 阶 1.6；窗口内孔隙率 [0.33105, 0.42756]，单调、可容许。
瞬态（h0 = 0.3 mm，24 h，625 步，0.4 s）：v 在 16 δ_C/v = 74 min 内进入稳态的 1e−6；启动损失 −16.4 nm（图 1c）；24 h 生长量对线性轨迹误差 −1.1e−5；累计 24 h 盐残差 3.9e−10、固体残差 1.2e−7。

**C. 收支（稳态，n_band 40，单位面积单位时间）**
- 盐：供给 ρW_L C_eq = 5.73913e−6 = 式 29 的 C_eq j/(1−C_eq)（−2e−10）；前沿 P_i 4.89021e−6，内部 S 8.48922e−7，溶解盐积累 0（稳态），残差 4.7e−10。M_front/M_precip = 0.8520819 对 F = 0.8520507（+3.7e−5）；M_int/M_precip = 0.1479181 对 1−F。
- 固体：前沿 ρ_cr(1−ε_new)v = P_i（2e−16）；内部固体体积率 S/ρ_cr = 3.9121e−10 m/s = 孔体积损失率 (ε_new − ε_final)v（−9e−10）。
- 液体：W_L − W_i = 7.0743e−10 m/s = 被 (3)–(4) 消耗的盐水体积 S/ρ（−1e−9）：**内部沉淀消耗的液体由底部多流入的 W 补足**，这就是方程里"丢失的孔体积"对应的液体去向。W_i = U_evap + (ρ_cr/ρ)(1−ε_new)v = 1.3333e−8 + 4.0752e−9；到达界面的水 ρW_i − P_i = j（精确，式 6+7；反向扩散的水与向界面扩散的盐相抵）。新孔空间 ε_new v = 1.70e−9 m/s 由随固体平流的盐水充填（2.04e−6 kg m⁻² s⁻¹），不在 W 内。
- 孔隙率：ε_final 0.331050 = 1−(1−ε_new)/F_num（3e−10）；对解析闭合 0.331026 差 7e−5（即 F 的离散误差）。

这些关系成立的条件：稳态行波（溶解盐无积累）、C(0) = C_eq 的饱和供给、ρ 常数、带内 a_v/D_eff 取出生值。ε_final = 1−(1−ε_new)/F 是其中前沿/内部分配的直接推论，不是独立输入。

**常密度闭合的体积不一致（识别，未修改）**：式 (3)–(4) 对每沉淀 S 的盐移走盐水体积 S/ρ，而物理上只有盐离开液相，其中水的份额 (1−C_eq)S = 6.25e−7 kg m⁻² s⁻¹ = j 的 3.9 % 实际留在（变稀的）液体里。这是 Prat 模型式 (4) 的理想化，不是求解器的守恒误差；它改变的是底部补给需求（≈ U_evap 的 4 %），不改变盐的分配与 F。真实的盐水体积变化应由 NaCl 的偏比容（量级 0.3 L/kg，文献值，本仓库 `materials/` 未收录，此处未用于任何计算）给出。

## 5. 图

- `figures/fig1_growth_analytic_vs_numerical.{png,svg}`：(a) h(t) 参考解对 h0 + v_Prat t；(b) 启动瞬态 v/v_Prat 对 t/(δ_C/v)（对数轴）；(c) h 与线性轨迹之差（nm），启动偏移 −16 nm 后以 +3.7e−5 的速度差线性漂移。
- `figures/fig2_reaction_zone_partition.{png,svg}`：(a) C*/C*_top 对 ξ/δ_C，解对式 19；(b) 前沿/内部分配对 F、1−F，内嵌 ε(ξ) 从 0.43 降到 0.3310 对闭合值。

## 6. 判定

| 项 | 目标 | 结果 | 判定 |
|---|---|---|---|
| 生长速度 解析 vs 数值 | 0.5 % | +3.7e−5（n40），+2.1e−5（n80），一阶收敛向 ~1e−5 平台 | PASS |
| 前沿/内部分配 | 0.5 % | +3.7e−5 | PASS |
| 界面过饱和 C*_top、剖面 | – | +3.7e−5；剖面 3.9e−5 C*_top | PASS |
| 盐残差 | < 1e−8 | 4.4e−10（稳态）、3.9e−10（24 h 累计） | PASS |
| 固体残差 | < 1e−6 | 1.9e−10（稳态）、1.2e−7（24 h 累计） | PASS |
| 孔隙率可容许 | – | [0.331, 0.4276] ⊂ (0, 1)，单调 | PASS |
| 反应带网格敏感性 | – | n_band 20→80：v 变 6e−5，ε_final 变 4e−5 | PASS |
| 大 Da 假设（h0） | – | Da 7.6e4，式 28 = 式 29 到 2e−16 | PASS |

**Step 2A 总判定：PASS。** 在 Prat 自己的线性化（V0）下，式 29/31/19 是方程组 (1)–(3)、(6)、(7)、(20) 稳态解的精确结果，误差 ~1e−5 级且随网格收敛。

**VERIFIED**：上述八项；a_v 系数 36 与 5/0.14 的差对 v 为 −0.059 %；ε_final = 1−(1−ε_new)/F 是分配的推论；液体与孔体积收支在方程内闭合，内部沉淀消耗的液体由底部流入补足。
**ASSUMED**：ρ = 1200 反标定；ε_new = 0.43 为出生孔隙率假设；d_c/d_sub = 0.14 经验比；a_v = 5(1−ε)/d_c 立方晶粒堆积；D_eff = √ε D_s；k_r、D_s 文献值（21 °C）；常密度盐水（式 4）；线性化（带内 ε、a_v、D_eff 取出生值——审查分支已证明按局部 ε 评价的自洽闭合在 112.5 µm 上使 v 高 0.5 %、细晶粒上差异更大，那是另一种闭合，不在本步范围）。
**NOT VERIFIED**：三维局部 Prat 生长（式 36 沿法向读取 j）；实际成核；实验中的晶粒尺寸选择；毛细供给极限与结壳；ε_new 的数值；与任何实验数据的一致性。

## 7. 复现

    PY=~/Documents/Humboldt/salt/salt_grow/.venv/bin/python
    OMP_NUM_THREADS=1 $PY runs/S_surface/step2a_prat/run.py
    $PY runs/S_surface/step2a_prat/plot.py
    OMP_NUM_THREADS=1 $PY -m pytest tests/test_step2a_prat.py -q
