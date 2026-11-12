# saltv2 对照 Prat (2026) 的审计：ε_birth、(1−F) 份额、salt_unsourced、分级结构、供液失效

日期 2026-10-25。审计对象：主检出 `codex/localized-nucleation-seed` @7179312 的 `saltv2/`（行号均指该版本）。
`../salt_grow_v2-baseline` @72ee682 与 `../salt_grow_v2-hierarchical` @6faf96f 在本文涉及的逻辑上与主检出相同
（差别只在：baseline 把 D_s、k_r、T、`salt_structure` 作为 case 输入传给 Prat 律；hierarchical 多一个未提交的
`layer_av` 开关，见 §4.3）。未改任何生产代码；反例脚本在 `reviews/prat_pnm_audit_checks/`。

Prat 原文：`../salt_grow/reference/prat2026.txt`（pdftotext，公式区乱码；式号按正文核对）与 `prat2026/core.py`
（逐式实现，37 个测试通过）。式号下文记作 P(n)。

## 0. 结论

| # | 问题 | 判定 | 性质 |
|---|---|---|---|
| 1 | ε_birth = 1 − (1−ε_s0)/F 是物理还是记账 | **记账 + 一个自洽的解读**。推导在 Prat 模型内成立（§1.1），但它对运输**零影响**：k_s、P_cth、r_p 与 F 无关，供液门从不读体素场（§1.2，反例 1）。只影响质量守恒、反应带容量 (ε_birth−ε_min)ρ_cr δ_C 和 D_eff ∝ ε^1.5。 | 未验证的物理假设 [A] + 已验证的实现事实 |
| 2 | (1−F) 盐是否重复计数或错放 | **无重复计数**。总盐一次进 `m_cell`，体积按 F 份额推进，翻转时 ε 由盐/体积反算；台账三条平衡独立闭合（残差 ≤ 1e-15，34 + 10 个测试通过）。两处小不一致：δ_C、Da* 用 ε_s0 = 0.43 而材料实际 ε_birth = 0.23–0.40（§2.3）；侧向生长的体素取"所占列"的 F 而非供液孔的 F（§2.4）。 | 实现无缺陷；参数解读有分岔 |
| 3 | `salt_unsourced` 是否允许非物理沉积 | **是，按构造允许**：无孔供液的簇仍按气相场通量蒸发、按 β(C_sat) 沉积，水与盐都无来源；上限 = 该簇自填到 ε_min 的盐量（§3）。生产记录中占 `salt_grown` 的 0.08–0.30 %；构造反例可达 100 %（反例 3）。 | 已记录的建模决定，非隐藏泄漏 |
| 4 | 分级结构是否违反 Prat 假设 | **违反两条**：(a) Prat 是单尺度 (d_c = 0.14 d_sub, a_v 常数)，主检出对每一层都用第 1 层的 F、δ_C，而自身层叠结构的晶粒每层减半 → 按 Prat 自己的 a_v = 5(1−ε)/d_c，F 在第 3–5 层跌破结壳极限 1−ε_s（反例 4）；(b) ε_s = 0.43 (Dai，单尺度) 与 r = 0.5 (Wijnhorst，分级) 配对，违反 `BirthClosure` 自己的 docstring（0.20 对应 r = 0.5）。hierarchical 的 B′ 修 (a) 但引入新假设（钳位 + 致密出生）。 | 参数与假设缺陷 |
| 5 | 供液失效、持续蒸发、沉淀是否自洽 | **不自洽**：供液门判定"h 处弯月面保不住"后，该面仍以气相场全通量蒸发，水照样从网络抽（`Q_evap` 无条件计入 `U_of` 和 `Q_draw`），盐全进反应带。c_k = 1e-9 时门给出可达高度 0.04 体素，簇仍蒸发到封闭为止（反例 5）。且门用的是出生时的解析层叠 (k_n(c_k), P_c(a_n))，与体素实际 k_s/P_cth 场脱节（§1.2）。 | 已记录的建模决定 (plan/16 M2) + 一个实现缺陷（死字段） |

优先修复见 §7。

## 1. ε_birth = 1 − (1 − ε_s0)/F

### 1.1 推导与 Prat 的对应（物理模型）

Prat 的盐结构内盐输运 P(1)：对流-扩散 + 体沉淀汇 a_v k_r (C − C_eq)；顶部边界 P(20)–(22)：界面沉淀 β(1−ε_s)k_r (C_top − C_eq) 造新层。
Da ≫ 1 时解为厚 δ_C = √(ε_s D_s,eff/((1−C_eq) a_v k_r)) 的指数边界层（P(19) 后文）。顶部过饱和 P(31)：

    C*_top = U_evap C_eq / [ β(1−ε_s)(1−C_eq)k_r + √((1−C_eq) a_v k_r ε_s D_s,eff) ]

分母两项就是两个盐汇：界面汇 β(1−ε_s)k_r 与带内汇 a_v k_r δ_C（a_v k_r ∫(C−C_eq)dz = a_v k_r C*_top δ_C）。
界面份额

    F = β(1−ε_s) / [β(1−ε_s) + a_v δ_C] = √Da* / (1 + √Da*)，   Da* = β²(1−ε_s)²(1−C_eq)k_r/(ε_s D_s,eff a_v)   P(27),(29)

（`prat2026/core.py:118-152, 249-259`）。所以 P(29)/(32) 的 (1−F) 份额在 Prat 模型里**沉在反应带孔壁上**，模型本身却保持 ε_s 不变
（原文 "our model considers that the precipitation only leads to new layers of similar porosity … self-replication regime"，Remaining Topic 段；Prat 承认 "one could expect that the salt precipitation on top leads to pore blocking or at least porosity reduction"）。

项目的闭合（`saltv2/growth/crust.py:117-139, 183-206`；`saltv2/coupling/driver3d.py:265-270, 572, 589`）：面交出 S = β j A；界面份额 F S 以 ε_s0 推进 V_n = F S/(ρ_cr(1−ε_s0))；带内份额 (1−F)S 在材料穿过带（时间 δ_C/V_n）时沉入 → 固体密度 = F S/(V_n A) + (1−F)S/(V_n A) = S/(V_n A) → 1 − ε_birth = (1−ε_s0)/F。
**在 Prat 的方程内这是精确的一维稳态结果**，前提是 ε_s0 解读为"界面新生层孔隙率"（ε_new）。对应的供液需求 P(33) 也相应改写（`driver3d.py:594-596`：ε_b + (ρ_cr/ρ_L)(1−ε_b)），与 Prat 一致。

### 1.2 实现中它改变了什么（已验证）

`reviews/prat_pnm_audit_checks/check1_eps_birth_is_bookkeeping.py`，d_sub = 34 µm，F = 0.759，ε_birth = 0.249：

| 字段 | F = None | F = 0.759 | 比 |
|---|---|---|---|
| eps_s | 0.430 | 0.249 | 0.58 |
| k_s [m²] | 1.770e-16 | 1.770e-16 | **1.00** |
| P_cth [Pa] | 1.26e5 | 1.26e5 | **1.00** |
| r_p [m] | 5.95e-7 | 5.95e-7 | **1.00** |
| D_eff [m²/s] | 5.65e-6 | 2.49e-6 | 0.44 |
| skin_capacity [kg/m²] | 5.33e-4 | 2.98e-4 | 0.56 |

原因：
- `BirthClosure.permeability` (`crust.py:170-176`) 与 `capillary_threshold` (`:178-181`) 只含 d_sub、n、c_k，**不含 ε**（ε_s0 = 0.20 与 0.43 给同一 k_s）。
- `CrustField.update_properties` (`crust.py:385-398`) 以 ratio = ε/ε_birth 演化 k、r、P_c；驱动在 `driver3d.py:977-979` 把 eps_s 和 eps_birth 同时覆盖为 `eps_from_salt` 再调用它 → ratio = 1，无变化。
- 供液门 `supply.first_failure_height` (`saltv2/growth/supply.py:33-61`) 用解析层叠 k_n = c_k 2^n a_n³/(12 d_1)、P_c = 2γ/a_n，签名里没有 F、没有 ε、没有体素场；`driver3d.py:651-652` 只传 (U, d_sub, s_foot, c_k, h_top)。
- `grep` 全 `saltv2/`：`.k_s`、`.P_cth`、`k_birth`、`P_birth` 除 `crust.py` 外**无任何读取者**。`CrustField.k_s / P_cth` 是只写不读的死字段。

所以 ε_birth 对"毛细供液 → 冻结"这条反馈链没有任何作用；它只通过反应带容量缩短封闭时间（0.56 倍）和降低 D_eff（0.44 倍）影响结果。
这与 2026-09-06 的任务 B 文档一致（"no grid length"），但**不能**把它说成"把 (1−F) 份额的物理后果带入模型"。

### 1.3 解读分岔（未验证的假设）

宏观 dh/dt = S/(ρ_cr(1−ε_final))，只约束 (1−ε)/F（`ideas/Prat生长律适用性审查_2026-10-25.md` §2.2 已指出）。
- ε_s0 = ε_new = 0.43：ε_final = 0.23–0.40 随 d_sub 变化；结壳极限 F ≤ 0.57 存在（`crust.py:200-205` 抛错），d_sub,crit = 6.0 µm（Prat 参数 D_s = 1.3e-9）。
- ε_s0 = ε_final = 0.43（Dai 测的是终态）：ε_new = 1 − F(1−0.43) = 0.57–0.68，极限消失；推进率同上式，与 F 无关。
Dai 的 0.43 是终态测量，所以**现行代码的解读与其数据来源相反**。两种解读的 dh/dt 相同，但 δ_C、Da*、带容量、D_eff 都不同。这是整条结壳路线的分岔点，现有数据不能判定。

## 2. (1−F) 份额的簿记：无重复计数

### 2.1 路径（实现）

1. 面盐 `salt_face = β_c · max(j,0)` (`driver3d.py:588`)，推进 `Vn_face = salt_f F_grid/(ρ_cr(1−ε_s0))` (`:589`)。
2. 候选格同时接收**全部盐** `salt_rate`（`:818-820`）与**体积率** `rate = saltF_rate/(ρ_cr(1−ε_s0))`（`:821-824`）；`m_cell += salt_rate dt`（`:965`），几何引擎按 `rate` 推进并用 `carry` 让溢出盐随体积同比例迁移（`geometry.py:245-302`）。
3. 翻转时 `eps_from_salt = 1 − m_cell/(ρ_cr V)`（`:973`），覆盖闭合的 ε；(0,1) 外抛错（`:974-976`）。
4. 前进面的 (1−F) 份额**不**经过 4b（`:762-775` 只处理 `frozen_f`）、不进 `dm_ext_rate`、不进 `salt_front_external`。

恒等式（解析核对）：Σ rate·ε_b,cell = Σ V_n A − Σ salt_f A/ρ_cr = Σ ε_b,grid V_n A，所以网络水账（`Q_fill_face`，`:596`）与实际新孔体积（`S_brine`，`:838`）一致，即使 ℓ > 0 平滑后 F_cell 变化。

### 2.2 台账证据（已运行）

- `tests/test_kinetic_share.py`（6 个）、`test_P4_front_and_ledger.py`、`test_P3_state_kernel.py`、`test_E1_supply.py`：**34 passed** (330 s)。
- `tests/test_P6_closed_loop.py`：**10 passed**；`test_regression.py`：3 failed（`ref_3d_small`、`ref_2d_small`、`ref_2d_cv`，差异 1e-15 相对，CLAUDE.md 记录的 frodo 平台钉扎），其余通过。
- 本文所有反例运行中 `salt_network`、`salt_interface`、`salt_precipitation` 残差 ≤ 4e-16。
- `prat2026/tests`：37 passed。

三条平衡的"独立性"值得保留：`salt_released_by_evaporation` 从气相场测、`salt_grown` 从生长律测，份额不为 1 会在这里现形（2026-09-06 前正是这样漏掉 5–24 %）。

### 2.3 不一致 A：δ_C 与 Da* 在 ε_s0 评价，材料在 ε_birth

`driver3d.py:255, 272, 288, 292` 全部传 `eps_s=EPS_S` = 0.43；出生材料 ε_birth = 0.23–0.40。带容量 (ε_birth − ε_min)ρ_cr δ_C (`crust.py:448`) 用致密 ε 但 δ_C 用疏松 ε：混合。
在 ε_new 解读下带内材料从 ε_s0（顶）过渡到 ε_birth（底），容量是下界；在 ε_final 解读下两处都应为 0.43。无论哪种，现在的组合都不是其中之一。

### 2.4 不一致 B：侧向生长的 F 取"所占列"

`exposed_faces` 的 mask 在气相候选格上（`geometry.py:152-164`），`F_grid` 按候选格所在列的孔径取值（`driver3d.py:540, 589`），注释说 "F of the COLUMN the new voxel will occupy"。
Prat 的 d_c = 0.14 d_sub 是晶粒对**供液孔**的继承；一个簇被孔 A 供液却长到孔 B 的列上，体素取 F_B、d_sub_B 的层叠。cv = 0 时差 0.3 %，cv > 0 时是假设 [A]。
`F_cluster`（`:468-471`，按供液孔平均）已算出但**未使用**（死代码，仅注释 `:845` 提到）。

### 2.5 D_s 不一致（参数）

主检出 Prat 律用 `prat2026/params.py:54` D_s = 1.3e-9，网络盐输运用 `materials/brine.py:34` D_NACL = 1.54e-9。baseline 已把 D_NACL 与 k_r(T) 传入（`driver3d.py:281-319` 那里）。Da* 差 18 %，F 差约 2 %（34 µm：0.760 → 0.744）。

## 3. `salt_unsourced`

### 3.1 机制（实现）

`driver3d.py:445-449`：簇的供液孔 = 湿底孔 (`cl_pores`)，否则足迹下任何孔 (`cl_pores_any`)，否则 `[]`。
`[]` 的簇：`beta_c = BETA_SAT`（`:533`，"keeps BETA_SAT, as it always has"）；`cl_pores` 空 → `h_new = 0`（`:682`）→ 全部面冻结 → 全部 β j 进反应带（4b）→ `S_k[c]` 无 `feed` → `S_unsourced`（`:856-862`）→ 台账 `salt_unsourced`（`:956`，`salt_network` 加项 `:309`）。
其水：`W_c` 从气相场累加（`:487-501`）但只在 `feed` 中的簇才进 `W_pore`（`:505-506`）；`Q_draw` 也只经 `feed` 进 `q`（`:712-713`）。**水和盐都没有源**，气相场照样给它 a_w(C_sat) 的 Dirichlet 面。

封闭后：湿体素向封闭体素蒸发 (`into_crust`，`:609-615`)，其盐经 `S_ext_k` 仍计到该簇 → 继续 unsourced，直到体素本体填到 ε_min（4c）。
所以上限 = (ε_birth − ε_min) ρ_cr V_cluster（簇把自己填实），幻影水 = 该盐/β。README 1503 行记录了这是有意选择（"试过让它拿不到盐，结果它永不封闭、永远蒸发"）。

### 3.2 数值

生产记录 `runs/E1_pre/crust_out/vapor_crust_*.npz`（λ 32/16/8，c_k 1e-5，2014 条件）：`salt_unsourced/salt_grown` = 3.0e-3、2.2e-3、1.3e-3、8.4e-4、1.2e-3；λ = 1：0。

构造反例 `check3_unsourced_cluster.py`：把成核格放到无孔列（monkeypatch `PoreOpenings.nucleation_seed`），pitch_factor 10 使簇 4 h 内到不了孔列：
`salt_supplied = 0`，`salt_unsourced = 5.46e-12 kg = salt_released_by_evaporation`，`water_evaporated_on_salt = 1.52e-11 kg`，4 个体素封闭，`M_liquid_in` 只等于裸孔蒸发。
（该反例须屏蔽 `salt_interface` 检查：种子盐 `per_pore(seed_mass)` 按足迹求和，无孔种子不从任何缓冲扣除，平衡立即以 −16 % 报警——这本身说明生产路径不会无声地产生无孔种子。）

### 3.3 判定

不是泄漏（具名、守恒检查通过），但**是非物理沉积**：违反模型自己的前提"结构的盐只来自被覆盖孔的结晶"（`driver3d.py:201-213`）。物理上该簇应失水（弯月面退缩、面停止蒸发）。替代"永不封闭、永远蒸发"的正确分支是让它退出湿集合（干透），不是给它幻影盐。

## 4. 分级结构 vs Prat 假设

### 4.1 Prat 的假设清单

单尺度立方晶粒 d_c = 0.14 d_sub（P(16)–(17)）；ε_s、a_v 随高度不变（自复制）；一维平面前沿、液饱和、准稳态；Da ≫ 1（tanh → 1，P(28)→(29)）；Pe ≪ Da；结构内 C ≈ C_eq。Remaining Topic 明确把 Wijnhorst 的多尺度结构列为未覆盖："need for a model predicting the cubic crystal size variation … in order to make use of eq 36"。

### 4.2 主检出：每层都用第 1 层的 F（反例 4）

`check4_layer_Da_star.py`：按 Prat a_v = 5(1−ε_s)/d_c 与项目层叠 d_c(n) = 0.14 d_sub 0.5^(n−1)（`crust.py:164-168`），ε_s = 0.43，Prat 参数：

| d_sub [µm] | F₁ | F 跌破 0.57 的层 | 该层 h_n [µm] | δ_C,1 / a₁(1) |
|---|---|---|---|---|
| 27.5 | 0.740 | 4 (F=0.50) | 220 | 0.56 |
| 34 | 0.760 | 4 (0.53) | 272 | 0.51 |
| 56 | 0.802 | 5 (0.50) | 896 | 0.39 |
| 112.5 | 0.852 | 6 (0.50) | 3600 | 0.28 |
| 1080 | 0.947 | > 7 (0.69 @ 7) | — | 0.09 |

主检出在所有层用 F₁（`driver3d.py:288, 540`）：对第 n ≥ 2 层，**它既不是 Prat（Prat 没有层），也不是自身层叠的推论**。`knowledge/prat2026.md` 最后一行已写明"a_v 与 Da* 应按最外层晶粒算，这是假设"，代码没有执行它。
另两点：δ_C,n ≈ a₁(n)（反应带与孔半宽同量级，第 2 层起带比孔厚）→ Prat 的连续介质带在分级结构里是亚 REV 量；Da(h_n) = h_n²/δ_C² ≥ 2.6e3，tanh → 1 没问题。

### 4.3 ε_s 与 r 的配对（参数缺陷）

`crust.py:103` docstring："eps_s = 0.43 (r = 1) or 0.20 (r = 0.5)"；`crust.py:60` R_HALVING = 0.5；`materials/salt.py:4` EPS_S = 0.43 → 主检出跑的是 **r = 0.5 与 ε = 0.43 的组合**，docstring 没有的组合。
Da*、δ_C、ε_birth、ρ_cr(1−ε_s)V_n 全部受影响（Da* ∝ (1−ε)/ε^1.5：ε 0.43 → 0.20 使 Da* 乘 (0.80/0.57)·(0.43/0.20)^1.5 ≈ 4.4）。baseline 引入 `EPS_S_HIERARCHICAL = 0.20` 与 `salt_structure` 记录（SINGLE_SCALE / HIERARCHICAL），并在 72ee682 的未提交修改里把 0.20 的出处降级为 [D]（"a/(d+a) in slit geometry, NOT a measured porosity"）——正确的方向，主检出未合并。

### 4.4 hierarchical 的 B′（`layer_av`，未提交）

实现 Da*_n = Da*_1 r^(n−1)（`crust.py` diff `layer_da_factor`），F 钳位在 F_min = (1−ε_s0)/(1−ε_min)，钳位处出生即 ε_min 并立刻 `seal`（`driver3d.py` diff），`shrink_at_birth` 把 r_p 按 √(ε_birth/ε_s) 缩。`tests/test_layer_av.py`：2 passed (3.6 s)。
这修正了 §4.2 的 (a)，但引入三个 Prat 之外的假设：(i) 带内沉淀的 (1−F) 份额在 F < F_min 时"把孔填实"而不是让 δ_C、a_v 自洽变化；(ii) 钳位 + 立即封闭是一个无中间态的开关；(iii) 致密出生体素的 r_p 用固定孔数律缩到 76 nm，与 Licsandru 的 3 nm 相差 25 倍。B′ 的可证伪预言（芯吸下转变粒径与通量无关）见 `ideas/plan/18_…md` §8。

## 5. 供液失效、蒸发、沉淀的一致性

### 5.1 门的构造（实现）

固定点 `driver3d.py:617-737`：`U_of(c,h) = (Q_evap[c] + cum(h))/A_open` (`:645-649`)，`Q_evap` = 簇全部湿面的水（`:591, 611`，"unconditional"），只有生长面的填孔项随 h 松弛；`reach = first_failure_height(U, d_sub, s_foot, c_k, h_top)`。
失效判据（`supply.py:51-56`）：s_foot + ρ g h + Σ μ U h_n/k_n > 2γ/a_n，即"h 处最细孔的弯月面保不住"。

### 5.2 失效后发生什么（反例 5）

`check5_supply_failure_vs_evaporation.py`，小箱 (λ 32, 2×2 孔, RH 7 %, 12 h)：

| c_k | 首次冻结 | 蒸发-only 需求下的可达高度 | 湿盐面蒸水 [kg] | E_total 冻结→末 | 封闭 |
|---|---|---|---|---|---|
| 1e-4 | 从不 | ∞ | — | — | 无 |
| 1e-7 | 0.63 h | ∞（需求 3.2e4 Pa < P_c 1.3e5 Pa） | 1.87e-10 | 1.000 → 0.000 | 全封，f_blocked 1 |
| 1e-9 | 0.05 h（第 1 行） | **0.29 µm = 0.04 体素**（需求 3.2e6 Pa ≫ P_c） | 1.70e-12 | 1.000 → 1.000 | 4 种子 0.06 h 封闭 |

c_k = 1e-9：门说连第一个体素顶都供不上，簇仍按 a_w(C_sat) 面蒸发、水记入 `Q_draw` 从孔抽、盐 5.4e-14 kg 全进带、封闭前 E 不变。
这是 plan/16 M2 的决定（"frozen face keeps evaporating"），注释 `driver3d.py:620`。**但它与门的前提矛盾**：判据的物理含义是液体退缩，退缩了就没有 a_w(C_sat) 的湿面。模型没有"退缩/干透"状态，于是把失效面当作湿面继续蒸发，直到沉淀把带填满——封闭时间 τ = (ε_birth − ε_min)ρ_cr δ_C/(β j) 用的是一个门已宣布供不上的 j。

### 5.3 第二处脱节

门从不读 `CrustField.k_s / P_cth`（§1.2）。封闭体素 ε = ε_min 的 k 按 `update_properties` 应降 ~(ε_min/ε_b)³ ≈ 1e-3–1e-4，P_cth 升 √(ε_b/ε_min) ≈ 4；站在封闭基底上的簇 (`cl_pores_any`，`:427-449`) 的水"穿过封闭体素"从孔抽（`Q_draw` 全额），门用的仍是出生层叠。沉淀 → 渗透率 → 供液的反馈在代码里**不存在**，尽管字段都在。

### 5.4 沉淀侧

冻结面全盐进带（4b）、穿壳面先填前沿体素本体再传给后面的带（4c）、基底盐进 z=0（4d）——三者互斥、守恒、dt 由三者之和限制（`:882-885`）。无问题。
`beta_pore` 的释放预算 `(1+γ)`（`:523-527`）正确：γ = ε_b F ρ_L C_sat/(ρ_cr(1−ε_s0)) 是每 kg 晶体带入新孔的溶解盐，量纲一致。

## 6. 可证实的缺陷 vs 未验证的假设

**实现缺陷（可证实）**
- D1 `CrustField.k_s / P_cth / k_birth / P_birth` 只写不读；供液门与体素场脱节（§1.2, §5.3）。
- D2 `F_cluster` 死代码（`driver3d.py:468-471`）。
- D3 主检出 ε_s = 0.43 与 r = 0.5 配对违反闭合自己的 docstring（§4.3）。
- D4 δ_C、Da* 在 ε_s0 评价而带容量在 ε_birth（§2.3）——两种解读下都不自洽。
- D5 D_s 两个值（§2.5）。
- D6 每层用 F₁（§4.2）与 `knowledge/prat2026.md` 自己声明的接口不符。

**建模决定（已记录、可争论）**
- M1 冻结面持续蒸发（plan/16 M2）——与门的前提矛盾（§5.2）。
- M2 无源簇按 β(C_sat) 沉积（§3）。
- M3 (1−F) 份额进出生孔隙率（任务 B）——Prat 方程内自洽，但依赖 ε_new 解读（§1.3）。
- M4 侧向体素取所占列的 F（§2.4）。

**未验证的物理假设**
- A1 ε_s0 是 ε_new 还是 ε_final（§1.3，决定结壳极限是否存在）。
- A2 分级结构上 Da* 按外层晶粒（B′）；δ_C 亚 REV（§4.2）。
- A3 B′ 的钳位/致密出生/r_p 缩放（§4.4）。

## 7. 优先修复

1. **D1**（小改，无物理改变）：要么让 `reach` 读簇内体素的 k_s/P_cth（把 `first_failure_height` 改成沿体素列积分，用 `crust.k_s` 与 `crust.P_cth` 的最小值；这会改所有记录），要么删掉死字段并在 docstring 写明"供液门只用出生层叠"。两者选一，不能保持现状——现状让读者以为沉淀→渗透率→供液闭环存在。配一个测试：`update_properties` 后 `reach` 必须变化（或字段不存在）。
2. **D3**：主检出合并 baseline 的 `salt_structure` 记录（SINGLE_SCALE: ε 0.43, r 1；HIERARCHICAL: ε 0.20 [D], r 0.5），禁止 0.43 × 0.5。改记录，需批量重跑。
3. **A1 分岔实验**（在修任何结壳物理之前）：同一小算例跑两种设置（ε_s0 = 0.43 作 ε_new；ε_s0 = 1 − F(1−0.43) 作 ε_new 使 ε_final = 0.43），对比 dh/dt（应相同）、封闭时间、E 平台。结果决定 §1.3、§2.3（D4）往哪边修。
4. **M1**：给失效面一个"退缩"分支：reach < z_face 的面退出湿集合（不蒸发、不沉积），而不是继续全通量蒸发。这改变 E(t) 与封闭时间，是物理决定，需用户拍板；至少把 `check5` 的表写入 README 作为现状声明。
5. **M2**：无源簇改为干透（同上一个分支），`salt_unsourced` 应归零；在此之前报告中标明它的上限 (ε_b − ε_min)ρ_cr V_cluster。
6. **D6/A2**：主检出要么接受 B′ (`layer_av`)，要么在 `BirthClosure` 里把 "F 对所有层相同" 写成显式假设并在 `describe()` 记录；`check4` 的表进 knowledge。
7. **D5**：主检出把 D_NACL、k_r(T) 传入 Prat 律（baseline 已做）。
8. **D2、M4**：删 `F_cluster`；`F_grid` 的列选择在 docstring 标 [A]。

## 8. 复现

```
PY=~/Documents/Humboldt/salt/salt_grow/.venv/bin/python
OMP_NUM_THREADS=1 $PY -m pytest tests/test_kinetic_share.py tests/test_P4_front_and_ledger.py tests/test_P3_state_kernel.py tests/test_E1_supply.py tests/test_P6_closed_loop.py -q
OMP_NUM_THREADS=1 $PY reviews/prat_pnm_audit_checks/check1_eps_birth_is_bookkeeping.py
PYTHONPATH=prat2026 $PY reviews/prat_pnm_audit_checks/check4_layer_Da_star.py
OMP_NUM_THREADS=1 $PY reviews/prat_pnm_audit_checks/check5_supply_failure_vs_evaporation.py     # ~3 min
OMP_NUM_THREADS=1 $PY reviews/prat_pnm_audit_checks/check3_unsourced_cluster.py                  # ~2 min
$PY reviews/prat_pnm_audit_checks/read_records_unsourced.py
```
