# V2 耦合模型的数学表述（V1：evaporation–wicking）

目标：把 Veran-Tissoires & Prat (2014) 的几何演化框架与 Prat (2026) 的局部生长律接成一个闭环模型，
使基底孔径 $d_{\rm sub}$ 通过 **(k_s, P_{cth})** 与 **Da\*** 两条通道决定盐结构形貌，形貌再反过来改变蒸发。
本文只写 wicking 情形（基底饱和、底部持续供液）；drying 扩展在 §10。

每条方程标注来源：**[2014]** Veran-Tissoires & Prat、**[2026]** Prat、**[2013]** Eloukabi、**[W25]** Wijnhorst-Prat-Shahidzadeh、
**[A]** 本项目假设、**[V]** 本项目在 `../salt_grow` 中已数值验证。参数出处见 `../knowledge/porous_salt_parameters.md`。

---

## 1. 域与状态

二维竖直截面 $(x,z)$，$z=0$ 为基底表面，$z$ 向上。域 $[0,L_x]\times[0,\delta_{\rm eff}]$。实验中样品管顶部开口于封闭箱，$\delta$ 是多孔面到管口的凹陷深度（2014：2.5/7.5/15 mm；Ben Amara：15 mm），**无盖板**；
$\delta_{\rm eff}=\delta+$ 管口外等效扩散长度（由 2014 Table 1 标定），顶边强加 $X_{v\infty}$。结构到达 $z=\delta$ 时模型几何失效（盐可伸出管口），计算停止。
每个体素 $i$ 带状态：

| 状态 | 含义 | 持有的场 |
|---|---|---|
| `GAS` | 空气 | $X_v$ |
| `SALT_WET` | 饱和盐水的多孔盐 | $p_\ell$；作为 $X_v$ 的 Dirichlet 边界 |
| `SALT_DRY` | 已干涸的多孔盐（正式语义 `SALT_DRY_OPEN`）：可透汽、不供液、不生长 | $X_v$（以有效扩散系数 $D_{\rm eff}$ 参与气相解） |
| `SALT_SEALED` | 孔隙被析盐填满（$\varepsilon_s\le\varepsilon_{s,\min}$）的盐：对蒸汽与液体都是墙 | 不进入活动扩散域（零通量） |
| `SALT_PARTIAL` | 正在被填充的界面格，填充度 $\phi_i\in[0,1)$ | 同 `GAS`，附加 $\phi_i$ |

基底 $z<0$ 在 V1 中不离散，只通过 $z=0$ 的边界条件出现（§6）。这正是日后接 PNM 的数据交换面。

每个盐体素在生成时记录 `eps_s, d_c, k_s, P_cth, r_p, t_birth, generation`，由 `structure_type`
（`PATCHY_POROUS_SALT` / `NANOPOROUS_CRUST`）的生成规则赋值；其中 `eps_s` 在内部湿前沿体素上随沉积密化减小（§5.1，2026-09-02 加入），其余不变。
`structure_type` 同时给定默认初始态：patchy 稀疏种子、crust 整层一格厚盐膜（成核密度是 d_sub 的第二条通道，V1 为输入，§7）。盐–基底脱离/缝隙状态暂不考虑（2026-09-02 决定）。
参数如何赋值、进入哪条方程、谁标定、影响多大，见 `v2_盐域参数在体素模型中的实现.md`。

## 2. 气相：蒸汽输运 [2014 式 6.1–6.4]

准稳态、稀释极限、纯扩散：

$$\nabla\cdot\big(\rho_g D_{\rm eff}\,\nabla X_v\big)=0,\qquad
D_{\rm eff}=\begin{cases} D_v & \text{GAS}\\ D_v\,\varepsilon_s/\tau_s & \text{SALT\_DRY}\ [A]\end{cases}$$

边界条件：

- 顶 $z=\delta_{\rm eff}$：$X_v=X_{v\infty}$（由 RH 定）。
- 侧边：周期。
- 湿盐面与无盐湿基底：$X_v = X_{ve}=a_w\,X_{v,\rm sat}(T)$，$a_w=a_w(C_{eq})\approx0.75$。
- **干盐内部的湿前沿**（`SALT_DRY` 与 `SALT_WET` 的交界）：$X_v = a_w\,a_{w,K}\,X_{v,\rm sat}$，Kelvin 因子
  $$a_{w,K}=\exp\!\left(-\frac{2\gamma\cos\theta\,M_v}{\rho_\ell R T\,r_p}\right)\quad[\text{2014 式 6.9 的量纲修正版，见注}]$$

局部蒸发通量（kg m$^{-2}$ s$^{-1}$）：$j=-\rho_g D_{\rm eff}\nabla X_v\cdot\mathbf n$，在所有湿面上取值。

> 注：2014 式(6.9) 与 Ben Amara 式(6) 印出的 Kelvin 式缺 $1/\rho_\ell$，量纲不闭合；上式为正确形式。
> 数值：$2\gamma\cos\theta M_v/(\rho_\ell RT)\approx1.1$ nm，故 $r_p$=1 nm 给 $a_{w,K}\approx0.33$，10 nm 给 0.90，100 nm 给 0.99。

把干盐当作"参与扩散的低扩散率区"而不是单独的 Robin 边界 [A]：壳的厚度、连通性、Kelvin 降活度都由同一套方程自动处理，
无需额外的壳阻力参数。$\tau_s$ 取 $\varepsilon_s^{-1/2}$（与 [2026] 液相曲折度同一假设）。

## 3. 局部生长律 [2026 式 32、36]

在每个湿界面点：

$$V_n = G\,F(Da^*)\,j,\qquad G=\frac{C_{eq}}{\rho_{cr}(1-\varepsilon_s)(1-C_{eq})},\qquad
F=\frac{\sqrt{Da^*}}{1+\sqrt{Da^*}}$$

$Da^*$ 由该点下方基底的 $d_{\rm sub}$ 决定（`prat2026.core.Da_star`）。$V_n$ 是**多孔盐结构/空气界面**的法向速度，
不是单晶生长速度；一格 = 一块多孔盐，固相体积分数 $1-\varepsilon_s$。

[V] 均匀 $F$ 只改时钟不改形貌，故 $d_{\rm sub}$ 对形貌的作用不在这一步。

## 4. 盐结构内供液 [2014 式 6.6–6.7；需求项来自 2026 式 33]

在 `SALT_WET` 内，不可压 Darcy：

$$\nabla\cdot\mathbf q=0,\qquad \mathbf q=-\frac{k_s}{\mu}\big(\nabla p_\ell-\rho_\ell\mathbf g\big)$$

边界条件：

- 底部（与基底接触的盐格，$z=0$）：$p_\ell=p_{\ell,0}(x)$（§6）。
- 湿界面：$\mathbf q\cdot\mathbf n=U_{\rm dem}$，
  $$U_{\rm dem}=\frac{j}{\rho_\ell}+\Big[\varepsilon_s+\frac{\rho_{cr}}{\rho_\ell}(1-\varepsilon_s)\Big]V_n$$
  第一项补蒸发，第二项填新生盐层的孔隙并携带将要结晶的盐。[2026 式 33] 忽略第二项在 $d_{\rm sub}$=25 µm 时低估需求约 25 %。
- 干盐的外表面（无蒸发）：$\mathbf q\cdot\mathbf n=0$。
- **盐/干盐交界（内部湿前沿）**：$\mathbf q\cdot\mathbf n=U_{\rm front}$，见 §5.1（2026-09-02 修正：此前写零通量，与 §2 的前沿蒸发矛盾，质量不守恒）。

$k_s$ 与 $P_{cth}$ 是**两个独立输入** [V，W25]：$k_s$ 由结构内最粗的孔道决定，$P_{cth}$ 由最细的外层孔决定。
不得用同一孔径经 Carman–Kozeny 同时给出两者。

## 5. 毛细兼容与干涸门控 [2014 式 6.8]

界面处所需毛细压力 $p_c(x)=p_{\rm atm}-p_\ell(x)$。定义

$$\Pi(x)=\frac{P_{cth}}{p_{\rm atm}-p_\ell(x)}$$

- $\Pi\ge1$：可供液，界面按 §3 生长。
- $\Pi<1$：该界面格转为 `SALT_DRY`：$V_n=0$，$\mathbf q\cdot\mathbf n=0$，气相中按 §2 的低扩散率区处理。

求解顺序（每个事件步）：解 §4 → 找 $\Pi<1$ 的界面格 → 标干 → 重解 §4，直到无新增干格。
这是 [2014] 的逐步标记法；严格形式是互补条件 $(P_{cth}-p_c)\ge0,\ (U_{\rm dem}-\mathbf q\cdot\mathbf n)\ge0$，两者乘积为零，
V1 不做。

**V1 假设 [A]**：状态机 `SALT_WET → SALT_DRY_OPEN → SALT_SEALED` **不可逆**（2026-09-04，plan/17 §3.1 冻结）。回湿、溶解、壳上二次成核都要求连续饱和度变量 $S$ 与吸液支 $p_{c,s}(S)$，两者都无数据，留到 V2。
（历史：2026-09-02 曾改为允许回湿——Wijnhorst 2024 的湿/干循环、Rodríguez-Rodríguez 2026 的 S ≈ 1.20 析出 / S = 1 溶解滞后支持"有滞后"。该分支在 `pre_experiment_vapor.run(allow_rewet=True)` 中保留，仅供复现旧结果；它需要 χ 与 τ_dry 两个纯数值参数才能压住湿/干抖动，见 plan/15 §7.2。）[2026] 讨论指出 wicking 下 $j$ 变小后压降也小，该机制能否解释后期减缓存疑；
Licsandru 2022 观察慢蒸发期盐霜仍饱和但与基底出现缝隙。回湿与缝隙都留作 V2 议题。

$\Pi$ 的一维闭式就是 [V] 已验证的控制群 $H^*=\rho_\ell k_s P_{cth}/(\mu \bar j)$；
[W25] 的结论换成这个语言：cauliflower 的 $\Pi\gg1$，crust 的 $\Pi\sim1$。

### 5.1 内部湿前沿：供液与沉积密化 [2026 式 33 同一推导 + A]

干盐体素与湿盐体素的交界面上有 Kelvin 降活度后的蒸发 $j_{\rm front}$（§2）。该面不移动（无新生孔隙），
局部析盐比例沿用 [2026] 的 $F(Da^*)$（其余 $1-F$ 经液相扩散回结构内部），于是

$$U_{\rm front}=\frac{j_{\rm front}}{\rho_\ell}\Big[1+F\frac{C_{eq}}{1-C_{eq}}\Big],\qquad
\dot S=\frac{F\,j_{\rm front}\,C_{eq}}{(1-C_{eq})\,\rho_{cr}}$$

$U_{\rm front}$ 是该面的液相 Neumann 边界（式 33 去掉移动界面的 $\varepsilon_s V_n$ 项），$\dot S$ 是单位面积的固相沉积速率（m/s）。
沉积进入前沿湿体素的孔隙：$\mathrm d\varepsilon_s/\mathrm dt=-\dot S/a$。$\varepsilon_s$ 降到 $\varepsilon_{s,\min}$（首轮扫 0.02 / 0.05 / 0.10 [V]）时该体素转 **`SALT_SEALED`**（2026-09-04，plan/17 §3.1、§4.4 更正：此前写作转 `SALT_DRY`，与"干盐可透汽"的语义冲突——被填满的体素对蒸汽也是墙）
（带着降低后的 $\varepsilon_s$，故其 $D_{\rm eff}$ 更低），前沿后退一格。

后果：壳厚平台之后 $E(t)$ 继续缓慢下降、壳缓慢增厚，即 Ben Amara 2025 的"壳变干变密"，不需要时变的开孔分数。
V1 假设 [A]：密化只改 $\varepsilon_s$，不改 $k_s$ 与 $r_p$；$\varepsilon_{s,\min}$ 要报敏感性。
Regime：前沿液体在 wicking 下由储液池经湿盐补给，永不中断；drying 下同一机制存在但供液终会中断（§10）。

## 6. 基底边界（= PNM 耦合面）

V1 把基底缩成 $z=0$ 上的两个给定场：

| 量 | V1 取法 | PNM 阶段由谁提供 |
|---|---|---|
| $p_{\ell,0}(x)$ | 常数 $p_{\rm atm}-\rho_\ell g L_s-\mu \bar j L_s/(\rho_\ell k_{pm})$（储液池在下方 $L_s$，[2014 §2] 几何） | PNM 表面孔的液压 |
| $C_0(x)$ | $C_{eq}$（wicking，$C_0\approx C_{\rm sat}$） | PNM 表面孔的浓度 |
| 无盐处蒸发 | Dirichlet $X_{ve}$（§2） | 同左，PNM 接收 $j_i$ 作为其蒸发边界 |
| 盐覆盖处 | 盐格底部 $p_\ell=p_{\ell,0}$ | PNM 接收 $U_{\rm dem}$ 作为该孔的抽液边界 |

反馈的完整含义：盐结构改变 (a) 裸露基底的 $j_i$（屏蔽）、(b) 被覆盖孔从"蒸发边界"变为"抽液边界"、(c) 抽液量含填孔项。
基底内部的浓度与流动在 PNM 里算，本模型不重复。

## 7. 成核

- **V1**：给定种子（位置、数量），**种子密度是形貌的第二条输入通道**：[2014 §5]、[2013 §4] 把 crusty/patchy 首先归于成核位置数量（粗基底 Pe 高、浓度峰局部化；细基底均匀、处处成核）。
  默认与 `structure_type` 绑定：patchy 间距 $\approx2d_{\rm sub}$ 量级 [2013 §4]；crust 整层盐膜（Ben Amara "crystalline veil"）。连片湿层不降蒸发，$E$ 的下降仍靠 §5。
- **V2**：由 PNM 表面浓度触发，$C_i>C_{\rm nuc}$，$C_{\rm nuc}>C_{\rm sat}$ [2014 §5 强调首次成核需过饱和]；
  已有盐结构的持续生长用 $C\approx C_{eq}$ [2026]。两个判据分开。

## 8. 界面推进与时间步

采用 `../salt_grow/salt_deterministic.py` 的确定性填充规则 [V]（式 36 直读版）：

- 界面格 $i$ 的容量 $V_{\max}=(1-\varepsilon_s)a^2$（单位厚度），已填 $\phi_iV_{\max}$。
- 固相沉积速率 $\dot V_i=(1-\varepsilon_s)\,\ell_i\,V_{n,i}$，$\ell_i$ 为该格的湿界面长度，$V_{n,i}$ 由该格上各湿面通量平均。
- 事件步 $\Delta t=\min_i (1-\phi_i)V_{\max}/\dot V_i$；所有在此 $\Delta t$ 内达容量的格同时翻转为 `SALT_WET`，其余按比例增加 $\phi_i$。
- 翻转后重解 §2、§4、§5。

[V] 一维平层下该规则逐位重现式(32)；对称种子严格保持镜像对称；被包住的气袋无通量永不填充（已知格子行为）。

**体素尺寸** [A]：$a=d_c=0.14\,d_{\rm sub}$。理由：纯 Laplacian 生长无内禀尺度、形貌不随网格收敛 [V]，
必须由物理长度正则化；多孔盐的"颗粒"就是晶粒。报告形貌时只能声称"晶粒尺度"分辨。

## 9. 无量纲量与输出

| 量 | 定义 | 用途 |
|---|---|---|
| $\Pi$ 或 $H/H^*$ | §5 | 形貌选择 |
| $Da^*(d_{\rm sub})$ | §3 | 速率修正 |
| $E=J_{\rm tot}/J_{\rm ref}$ | 总蒸发/无盐时 | 反馈强度；patchy 期望 $\approx1$ 至 $>1$ [2014 Fig 27、2013 Fig 2f]，crust 期望 $\ll1$ |
| 覆盖率 $\eta$ | 被盐覆盖的基底长度分数 | patchy 恒定、crust→1 |
| $\bar h,\ h_{\max}$ | 平均/最大高度 | 对 Ben Amara Fig 7/12 |
| 箱计数维数 | 形貌 | 对 [W25] 1.71（2D） |
| $r_p/\ell_K$，$\ell_K=2\gamma M_v/(\rho_\ell RT)$ | Kelvin 数 | 壳阻力敏感性 |

## 10. drying 扩展（V2）

加基底域 $z<0$：饱和度 $S(\mathbf x,t)$、$p_\ell$、$C$ 的 Richards + 对流扩散 + 沉淀汇；保持曲线用 Hidri 2013 的 van Genuchten
（$n=10$，$S_c=0.1$，$P_{c,\rm ref}=6(1-\varepsilon)\gamma/(\varepsilon d)$）[W25 式 2]。§6 的 $p_{\ell,0}$、$C_0$ 改由该域给出，
其余方程不变。生长停止在 drying 下由盐量与液压耗尽自然出现 [W25 式 6，2026 讨论]。

## 11. 参数清单

| 参数 | V1 默认 | 性质 | 标定/敏感性 |
|---|---|---|---|
| $\varepsilon_s$ | 0.43（patchy）、0.25（crust） | M | 扫 0.20–0.43 |
| $d_c/d_{\rm sub}$ | 0.14 | M | 固定 |
| $k_r, D_s, E_a, \beta$ | `prat2026/params.py` | L | 固定 |
| $k_s$ | 待定：patchy 用最粗层立方律 [W25]，crust 用 1e−22~1e−21 m² [V] | **B** | 对 [2014 Fig 26] 序列与 Ben Amara #1/#4 标定 |
| $P_{cth}$ | 待定：$2\gamma\cos\theta/r_{p,\rm outer}$；crust 用 $r_p$=3–4 nm | **B** | 与 $k_s$ 联合标定 |
| $r_p$（Kelvin） | 与 $P_{cth}$ 同一外层孔径 | A | 开关项 |
| $\tau_s$ | $\varepsilon_s^{-1/2}$ | A | 敏感性 |
| $\varepsilon_{s,\min}$（前沿密化下限） | 0.05 | A | 扫 0.02–0.10 |
| $\gamma\cos\theta$ | 0.075 N/m，$\theta$=0 | L | 固定 |
| $a_w(C_{eq})$ | 0.75 | L | 固定 |
| $\rho_g D_v$ | 3e−5 kg m$^{-1}$ s$^{-1}$（沿用 salt_lattice） | L | 外部对流用有效 $\delta$ 吸收 [V] |
| $\delta$（管口凹陷深度，无盖板）→ $\delta_{\rm eff}$ | 15 mm（Ben Amara）、2.5/7.5/15 mm（2014）+ 管口外等效长度 | 实验几何 + S1 标定 | $E-1$ 随 $\delta$ 变，须同报 [V] |
| $p_{\ell,0}$ | §6 | 实验几何 | 固定 |
| 种子 | 位置、数量 | A | 敏感性 |

需要标定的只有 $(k_s,P_{cth})$ 到 $d_{\rm sub}$ 的映射，其余全部有出处。这是模型唯一的自由度，须在论文中如此陈述。

## 12. 验证序列

1. **一维锚点**：平层推进 = `core.growth_rate`（1e−12）；$U_{\rm dem}$ 与式(33) 一致。
2. **守恒**：沉积盐质量 = $\int G F j\,\mathrm dt$ + $\int \dot S\rho_{cr}\,\mathrm dt$（前沿）；水质量 = $\int j\,\mathrm dt$ + 填孔；内部前沿的抽液 = $\int U_{\rm front}\mathrm dt$。
3. **对称**：镜像种子给镜像形貌。
4. **[2014 Fig 26–28] 定性序列**：固定 $P_{cth}$、扫 $k_s$，patchy → crusty → blocked；patchy 期 $E\approx$ 常数、crust 期 $E\to$ 无盐面积分数。
5. **[2014 Table 1]**：外部有效扩散长度用 $R_{\rm ext}$ 常数校核。
6. **Ben Amara #4**：patch 高度 12–13 mm 量级、早期 $\bar h(t)$；**#1**：壳高趋稳、$E\approx1/14$（需 §2 的干盐扩散 + Kelvin）。
   两者都是 drying，只比早期饱和阶段并注明。
7. **[W25]**：2D 箱计数维数与 1.71 比较（定性）。
8. **Norouzi Rad 2013**：`prat2026/fig08` 的质量速率在耦合模型中复现（$\eta$ 不再是输入而是输出）。

## 13. 数值实现要点

- 两个 Laplace 型方程用同一有限体积求解器：面扩散率取调和平均，Dirichlet 半格系数 2 [V，efflorescence3d]；2D 规模 ≤1e5 未知数，稀疏直接解即可。
- 干涸标记循环（§5）在每个事件步内收敛通常 1–2 次。
- 气相域随盐生长缩小、干盐域扩大，热启动需按空间场重映射 [V]。
- 二维会夸大自屏蔽（[V] 2D $E$=2.04 对 3D 1.10）。求解器自 S1 起维数无关，默认 2D，$E$ 的定量报告用 3D 算一次（2026-09-02 决定）。

## 14. 与 ../salt_grow 现有代码的关系

| 现有 | 复用 | 不复用 |
|---|---|---|
| `prat2026/core.py` | §3 全部 | — |
| `salt_lattice.py` / `salt_deterministic.py` | §8 填充规则、FV 求解器、`GrowthLaw` 接口思想 | 随机容量（等价于非均质 $\varepsilon_s$，V1 不要）、干格零通量（改为 §2 低扩散率区） |
| `growth_laws.capillary_rise_limit_decoupled` | §5 的一维闭式作检验 | 一维闭式本身（V1 解二维压力场） |
| `efflorescence3d` | 半格 Dirichlet、热启动、VOF 思想 | 三维、VOF（V1 用格子填充） |
