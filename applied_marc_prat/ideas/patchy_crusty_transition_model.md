# 多孔介质表面盐结晶 Patchy–Crusty 转变：机制假设、跨领域启发与数值建模方案

> 2026-09-02 修订：按 `plan/11_转变模型_可行性与实施阶梯.md` 的审查改了五处（§12 生长闭合、§13–16 水力判据与 d_b 缩放、§18/20 时间步、§19 二维贯通判据、§23 验证 regime），并把 §28 的开发顺序与 `saltv2/` 的实现状态对齐。改动处均标 **[修订]**。
> **Regime：全文第一版为 saturated evaporation-wicking**（§6.1）；drying 只出现在 §23 的数据讨论里。

## 1. 研究目标

我们希望建立一个能够**自发预测盐结晶形态**的模型，而不是预先规定某个条件下为 patchy efflorescence、另一个条件下为 crusty efflorescence。

核心科学问题可以表述为：

> **What determines whether initially discrete salt nuclei remain isolated or coalesce into a continuous crust?**

即：

> 初始离散盐核在什么条件下会保持为相互分离的 patch，什么条件下会发生横向扩展、合并并最终形成连续 crust？

我们的核心判断是：

\[
\boxed{\text{patchy–crusty transition} \neq \text{单一毛细作用或单一颗粒尺寸控制}}
\]

而更可能是以下过程之间竞争的结果：

\[
\boxed{
\text{nucleation}
+
\text{hydraulic-limited growth}
+
\text{evaporation screening}
+
\text{salt-solution pumping}
+
\text{coalescence/percolation}
}
\]

---

# 2. 为什么优先研究 patchy–crusty 转变

与一开始就精确测量盐壳孔隙率、渗透率相比，patchy–crusty 转变更基础，因为：

1. **形态决定后续输运路径**：
   - patchy：气相仍能直接接触大量裸露多孔介质表面；
   - crusty：连续盐壳可能大面积覆盖并重构蒸发边界条件。

2. **形态决定盐结构中的流动方式**：
   - patchy：液体可能被局部 pumping 到少量 efflorescence；
   - crusty：供液可能分布在大面积连续盐层中。

3. **形态决定蒸发反馈**：
   - patch 会在局部增强蒸发并屏蔽周围表面；
   - crust 更接近一个新的连续 porous layer。

4. 如果模型能够预测 transition，那么盐层孔隙率、渗透率等参数可以在第二阶段再加入，而无需一开始把所有细节都耦合进去。

---

# 3. 已有盐结晶研究给出的直接线索

## 3.1 Fine beads → crusty；coarse beads → patchy

已有实验在 drying 与 wicking 两类 regime 中都观察到：

- 细颗粒基质更容易形成 crusty efflorescence；
- 粗颗粒基质更容易形成 patchy efflorescence。

典型实验系统包括：

\[
1\text{–}50\ \mu m \quad \rightarrow \quad \text{crusty}
\]

\[
400\text{–}600\ \mu m \quad \rightarrow \quad \text{patchy}
\]

一个尤其有价值的 drying 实验还显示：

- 在 fine-bead substrate 顶部增加 1–5 mm coarse layer，能够推动 morphology 向 patchy 方向变化。

这些相关性说明：

> 真正控制 transition 的很可能是**近表面的输运条件、成核密度和新生盐结构本身的微结构/水力性质**，而不是整个样品的平均颗粒尺寸。珠径是实验控制量，不应直接当作盐结构渗透率的代理量。

---

## 3.2 初始成核位置不是完全随机

Prat/Veran-Tissoires 一系列工作表明：

- porous medium 内部的 pore-scale velocity fluctuations；
- meniscus distribution；
- evaporation flux heterogeneity；

会导致表面盐浓度出现空间 maxima，从而产生离散初始 nucleation sites。

因此不建议在模型里简单“随机撒盐核”，更合理的做法是让 substrate 中的 transport heterogeneity 自然产生 nucleation locations。

---

## 3.3 已有 efflorescence 会改变周围输运场

一个已经形成的 patch 会产生两个重要反馈：

### (1) Pumping effect

patch 表面持续蒸发，因此必须从下面的 porous substrate 获得液体补给：

\[
J_{\rm evap}^{patch}>0
\]

从而诱导：

\[
\mathbf u \rightarrow \text{existing patch}
\]

盐也随液体被输送到已有 patch：

\[
C\mathbf u \rightarrow \text{patch}
\]

因此已有结构会进一步得到更多盐供应。

### (2) Vapor screening

突出于表面的盐结构会改变 gas-phase vapor concentration field。

patch 顶部通常具有更高局部蒸发通量：

\[
j_e^{top}\uparrow
\]

而 patch 周围裸露表面会受到 screening：

\[
j_e^{bare}\downarrow
\]

于是这些区域的 upward liquid flow 和 salt accumulation 同时下降，从而抑制新核形成。

因此所谓 exclusion zone 不一定需要人为指定；它可以是 pumping + screening 的自然结果。

---

# 4. 跨领域可以借鉴的核心思想

## 4.1 Thin-film growth：island → coalescence → continuous film

这是目前最有启发性的类比之一。

经典 Volmer–Weber 型薄膜生长过程：

\[
\text{nucleation}
\rightarrow
\text{isolated islands}
\rightarrow
\text{island growth}
\rightarrow
\text{coalescence}
\rightarrow
\text{continuous film}
\]

与盐结晶可以对应为：

\[
\text{salt nucleation}
\rightarrow
\text{isolated efflorescence}
\rightarrow
\text{patch growth}
\rightarrow
\text{patch interaction}
\rightarrow
\text{crust}
\]

这启发我们：

> patchy 与 crusty 可能不是两种完全不同的晶体生长模式，而是同一 nucleation–growth 系统在“是否发生 lateral coalescence”上的不同结果。

---

## 4.2 KJMA / Avrami：用成核和生长描述覆盖率

定义单位面积成核速率：

\[
I(t) \qquad [\mathrm{m^{-2}s^{-1}}]
\]

若在 \(\tau\) 时刻产生的一个 patch 到 \(t\) 时刻面积为：

\[
A_p(t,\tau)
\]

则 extended coverage：

\[
X^*(t)=\int_0^t I(\tau)A_p(t,\tau)d\tau
\]

真实表面盐覆盖率可写成类似 KJMA 的形式：

\[
\theta_s(t)=1-\exp[-X^*(t)]
\]

如果 patch 近似圆形：

\[
A_p=\pi R_p^2
\]

则 patchy–crusty 可以从一个“离散 island 是否合并成连续覆盖”的问题重新定义：二维截面用 \(\theta\to1\)，三维表面才用 spanning / percolation 判据（§19）。

---

## 4.3 Capture zone / exclusion zone

薄膜 island growth 中已有岛会捕获周围扩散粒子，使附近新核形成受抑制。

盐系统中可以定义：

\[
L_{\rm ex}
\]

表示已有 efflorescence 对周围新成核的有效抑制长度。

平均 nuclei spacing：

\[
L_n\sim N_n^{-1/2}
\]

可定义一个形态参数：

\[
\boxed{\Lambda=\frac{L_{\rm ex}}{L_n}}
\]

物理解释：

- \(\Lambda\) 大：少量已有 patch 强烈抑制周围成核 → 更容易保持 isolated → patchy；
- \(\Lambda\) 小：大量 nuclei 可以并存 → lateral overlap/coalescence 概率高 → crusty。

但在最终模型中，\(L_{\rm ex}\) 更适合作为**模型输出**，而不是输入。

---

## 4.4 Reactive precipitation：Pe–Da morphology map

在反应沉淀、多孔介质堵塞等领域，经常用：

\[
Pe=\frac{UL}{D}
\]

描述 advection / diffusion 竞争，用：

\[
Da=\frac{\text{reaction/precipitation timescale}}{\text{transport timescale}}
\]

描述 precipitation / transport 的竞争。

这个思路可以直接迁移到盐结晶，但单纯 \(Pe-Da\) 可能还不足以描述 patchy–crusty，因为盐结构形成后会反馈改变 flow 和 evaporation field。

因此更适合最终形成：

\[
\boxed{
\mathcal M=F(Pe_p,Da_s,\Pi_h,\Lambda)
}
\]

其中：

- \(Pe_p\)：substrate 中盐输运 Pe；
- \(Da_s\)：盐结构内部 precipitation/growth Da；
- \(\Pi_h\)：盐结构 hydraulic accessibility；
- \(\Lambda\)：emergent exclusion / nuclei spacing ratio。

---

## 4.5 DLA / electrodeposition：diffusion screening

电沉积和 diffusion-limited aggregation 的共同特征是：

\[
\text{geometry}
\rightarrow
\text{transport field}
\rightarrow
\text{local flux}
\rightarrow
\text{growth}
\rightarrow
\text{new geometry}
\]

盐 efflorescence 同样具有：

\[
\boxed{
\text{salt geometry}
\leftrightarrow
\text{evaporation field}
}
\]

突出结构处 local flux 增强，相邻区域受到 screening，这会导致生长进一步 localization。

因此 patchy 可以被理解为一种：

> **transport-field localization instability**。

---

# 5. 推荐的最小可行模型（MVP）

建议第一版模型由三个区域组成：

\[
\boxed{\Omega_p+\Omega_s(t)+\Omega_g}
\]

其中：

- \(\Omega_p\)：原始 porous substrate；
- \(\Omega_s(t)\)：随时间生长的 porous salt structure；
- \(\Omega_g\)：外部气相。

二维截面示意：

```text
             dry air

      Ωg : vapor diffusion
 ───────────────────────────
             ███
            █████       ← salt patch Ωs
             ██
────────────────────────────  original surface
      Ωp : porous substrate

              ↑
        brine / reservoir
```

第一阶段完全可以做 2D；后续再升级为 surface \(x-y\) + height \(h(x,y,t)\) 的 2.5D 模型。

---

# 6. Substrate 中的控制方程

## 6.1 第一版假设 saturated evaporation-wicking

如果底部持续液体补给充足，可以暂时认为 substrate 始终 saturated。

因此第一版无需 Richards equation。

### Darcy flow

\[
\nabla\cdot\mathbf u=0
\]

\[
\mathbf u=-\frac{K(\mathbf x)}{\mu}
\left(\nabla p-\rho_l\mathbf g\right)
\]

其中 \(K(\mathbf x)\) 不建议完全均匀。

可以使用 correlated permeability field：

\[
K(\mathbf{x})=K_0\exp[\sigma_K\xi(\mathbf{x})]
\]

其中：

\[
\langle\xi\rangle=0
\]

相关长度可取：

\[
l_c\sim d_p
\]

用来在 continuum framework 中表示 pore-scale disorder。

---

# 7. 盐输运方程

substrate 中求解：

\[
\boxed{
\phi\frac{\partial C}{\partial t}
+
\nabla\cdot
\left(C\mathbf u-\phi D_{\rm eff}\nabla C\right)
=-R_{\rm prec}
}
\]

其中：

\[
D_{\rm eff}=D_s f(\phi,\tau)
\]

底部 reservoir：

\[
C=C_0
\]

NaCl 不挥发，因此在尚未结晶的 bare surface：

\[
\boxed{
\left(C\mathbf u-\phi D_{\rm eff}\nabla C\right)\cdot\mathbf n=0
}
\]

这使蒸发引起的 upward liquid flow 导致表面盐浓度不断积累。

---

# 8. 成核模型

第一版可以从最简单 criterion 开始：

\[
S(\mathbf x,t)=\frac{C}{C_{\rm sat}}
\]

如果：

\[
S>S_{\rm nuc}
\]

则允许 nucleation。

可以先设：

\[
S_{\rm nuc}\approx1
\]

再通过实验 induction time 校准。

更一般可以使用概率模型：

\[
P_{\rm nuc}
=1-\exp[-J(S)\Delta A\Delta t]
\]

其中 \(J(S)\) 是 supersaturation-dependent nucleation rate。

第一版不需要显式设置 exclusion radius，而让 pumping + vapor screening 自然产生 nucleation suppression。

---

# 9. Gas phase vapor diffusion

外部空气中第一版可采用准稳态 diffusion：

\[
\boxed{\nabla\cdot(D_v\nabla c_v)=0}
\]

如果 \(D_v\) constant：

\[
\nabla^2c_v=0
\]

远场：

\[
c_v=c_{v,\infty}
\]

所有湿润 surface（bare wet substrate + wet efflorescence surface）：

\[
c_v=a_w(C)c_{v,\rm sat}(T)
\]

局部蒸发通量：

\[
\boxed{
j_e=-D_v\nabla c_v\cdot\mathbf n
}
\]

由于 salt geometry 显式进入 gas domain，vapor screening 会自动出现。

---

# 10. Pumping 的实现

对于第 \(k\) 个 salt cluster，计算总蒸发量：

\[
Q_k=\int_{\Gamma_k}j_e\,dA
\]

对应液体补给流量：

\[
Q_{l,k}=\frac{Q_k}{\rho_l}
\]

将其作为 salt footprint \(F_k\) 上的液体 sink / boundary flux：

\[
\boxed{
\mathbf u\cdot\mathbf n
=
\frac{Q_{l,k}}{|F_k|}
}
\]

第一版可均匀分配。

更一般：

\[
u_n(\mathbf x)=\frac{Q_{l,k}}{|F_k|}w(\mathbf x)
\]

满足：

\[
\int_{F_k}w\,dA=1
\]

这样已有 patch 会主动重定向 substrate 中的液体和盐输运。

---

# 11. Exclusion zone 应作为 emergent quantity

已有 patch 周围同时发生：

### Vapor screening

\[
j_e^{bare}\downarrow
\]

### Flow redirection

\[
u^{bare}\downarrow
\]

因此：

\[
\frac{\partial C^{bare}}{\partial t}\downarrow
\]

如果：

\[
C^{bare}<C_{\rm nuc}
\]

新核自然无法形成。

于是模型中可以后处理得到：

\[
L_{\rm ex}^{sim}
\]

并与实验/已有 scaling 比较。

---

# 12. 盐结构的生长律

已有工作指出，多孔 salt structure 的局部生长速度可以与局部蒸发通量建立线性关系：

\[
V_n\propto j_e
\]

**[修订]** 原稿的工程闭合 \(V_n=\eta C_s j_e/[\rho_l\rho_{\rm NaCl}(1-\varepsilon_s)]\) 量纲不闭合：盐通量应为 \(J_{\rm salt}=\eta C_s j_e\)（kg m\(^{-2}\) s\(^{-1}\)），不含 \(1/\rho_l\)，多出的一项使 \(V_n\) 小 1200 倍。
从第一版起直接用 Prat (2026) 的严格形式（`prat2026/core.growth_rate`，已验证到 0.08 %）：

\[
\boxed{
V_n=G\,F(Da^*)\,j_e,\qquad
G=\frac{C_{\rm eq}}{\rho_{\rm cr}(1-\varepsilon_s)(1-C_{\rm eq})},\qquad
F=\frac{\sqrt{Da^*}}{1+\sqrt{Da^*}}
}
\]

其中 \(C_{\rm eq}/(1-C_{\rm eq})\) 是每千克蒸发水携带的盐，\(F\) 是在顶部局部析出的比例（其余经液相扩散回结构内部），\(Da^*\) 由盐结构的 \(\varepsilon_s\)、比表面积（经 \(d_c=0.14\,d_{\rm sub}\)）、\(D_{s,\rm eff}\)、\(k_r\) 给出。
供液需求同样取 Prat 式 (33)：\(U_{\rm dem}=j_e/\rho_l+[\varepsilon_s+(\rho_{\rm cr}/\rho_l)(1-\varepsilon_s)]V_n\)，第二项是新生孔隙的填充与将要结晶的盐。

---

# 13. 为什么必须加入 hydraulic accessibility

如果仅使用：

\[
V_n\propto j_e
\]

突出结构顶部由于 vapor flux 较大，很容易持续向上生长，模型可能天然偏向 patchy。

要允许 crust 出现，需要考虑：

> efflorescence 内部液体是否能持续输送到高处。

把 salt structure 看成 porous medium：

\[
\mathbf u_s
=-\frac{K_s}{\mu}
(\nabla p_s-\rho_l\mathbf g)
\]

盐结构 capillary pressure：

\[
p_{c,s}\sim\frac{2\gamma\cos\theta}{r_s}
\]

输液到高度 \(h\) 所需压力：

\[
\Delta p_{\rm req}
=
\Delta p_{\rm visc}+\rho_lgh
\]

其中：

\[
\Delta p_{\rm visc}
\sim
\frac{\mu q_s h}{K_s}
\]

定义：

\[
\boxed{
\Pi_h
=
\frac{\Delta p_{\rm req}}
{\Delta p_{\rm cap}}
}
\]

**[修订]** 可用毛细压力在两种 regime 下不同：

\[
\Delta p_{\rm cap}=
\begin{cases}
p_{c,s}-\big(p_{\rm atm}-p_{l,\rm surf}\big) & \text{wicking：基底饱和，表面液压由储液池深度与基底黏性压降给出，}p_{c,pm}\text{ 不进入}\\[4pt]
p_{c,s}-p_{c,pm} & \text{drying：基底失饱和后液压降到 }-p_{c,pm}\text{，与盐的吸力竞争（Eloukabi 2013 式 11）}
\end{cases}
\]

第一版是 wicking，用上式第一行。这意味着在 wicking 里**基底孔径不能以 \(p_{c,pm}\) 的形式进入 \(\Delta p_{\rm cap}\)**。基底渗透率仍可通过饱和基底内的黏性压降影响 \(p_{l,\rm surf}\)，珠径也可通过输运非均匀性影响成核；但不能把 drying 的“基底吸力与盐吸力竞争”搬到 wicking。2014 在 wicking 下仍观察到 5–50 µm 珠结壳，因此其 crust 必须从盐微结构、基底输运/成核密度等通道解释，而不能归因于 \(p_{c,pm}\)。
黏性项与重力项的量级（K§8）：patchy 型盐结构 \(\rho_l g\approx1.2\times10^4\) Pa/m 远大于 \(\mu q_s/K_s\approx 30\) Pa/m，上限由重力定（米级）；纳米孔 crust 黏性项 \(\sim10^{10}\) Pa/m，上限由黏性定（毫米）。

---

# 14. Hydraulic number 对 morphology 的意义

## Patchy 倾向

如果：

\[
\Pi_h\ll1
\]

顶部供液容易：

\[
j_e^{top}\uparrow
\]

\[
V_n^{top}\uparrow
\]

因此 vertical growth 强：

\[
\boxed{\text{isolated vertical structures → patchy}}
\]

---

## Crusty 倾向

如果：

\[
\Pi_h\to1
\]

顶部供液受限。

可定义 hydraulic accessibility：

\[
A_h
=
\frac{1}{1+\exp[\alpha(\Pi_h-1)]}
\]

然后：

\[
\boxed{V_n=A_h(\Pi_h)\,G\,F(Da^*)\,j_e}
\]

顶部随高度增加：

\[
A_h\downarrow
\]

vertical growth 被抑制。

但靠近 substrate 的 lateral surface 因 \(h\approx0\)，仍容易获得液体供应，因此：

\[
V_{\parallel}>V_{\perp}
\]

更容易发生 lateral coalescence：

\[
\boxed{\text{lateral spreading → coalescence → crust}}
\]

---

# 15. Fine beads → crusty 不能由单一孔径缩放解释

实验相关性 \(d_b\downarrow\Rightarrow\text{crusty}\) 是可靠现象，但下面这条简单因果链在本问题的毫米尺度上不成立：

\[
d_b\downarrow\Rightarrow r_s\downarrow,
\qquad K_s\propto r_s^2,
\qquad p_{c,s}\propto r_s^{-1}
\Rightarrow \Pi_h\uparrow .
\]

把这组缩放代入含重力的供液高度上限

\[
H_{\lim}
=
\frac{\Delta p_{\rm cap}}
{\rho_l g+\mu q_s/K_s}
\approx
\frac{p_{c,s}}
{\rho_l g+\mu q_s/K_s}
\]

得到：500 µm 基底对应约 0.7 m，25 µm 基底仍约 8 m，都远大于毫米级盐结构，所以 \(\Pi_h\ll1\)。要使 \(\Pi_h\sim1\)，\(K_s\) 需比继承孔径缩放给出的值再低 4–5 个数量级。这不是同一微结构的连续缩小，而是**盐层微结构发生 regime switch**；Veran-Tissoires & Prat 2014 §6.3 也指出 crust 的差异不只是孔径变化，Licsandru 2022 测得的 3–4 nm 孔与 \(\varepsilon_s=0.25\) 支持这一点。

因此保留两条彼此独立、可检验的候选通道：

1. **成核通道（V2）**：细基底改变表面盐输运非均匀性和成核密度，使平均核间距变小，促进横向合并；
2. **盐微结构通道（V4）**：某些条件生成低 \(K_s\) 的致密 crust，使 \(\Pi_h\) 接近 1，水力限制把竖向生长转成横向生长。

H3 检验的是第二条机制的条件效应，而不是证明 \(d_b\) 必然决定 \(K_s\)。在相图中，\(K_s\)（以及决定可用吸力的 \(p_{c,s}\) 或 \(r_s\)）必须是独立轴；在得到微结构测量或独立闭合前，不写 \(d_b\to K_s\) 的确定映射。

---

# 16. Coarse beads → patchy 的条件性解释

粗珠体系更常见 patchy 仍是需要模型解释的实验事实，不能简单由 \(d_b\uparrow\Rightarrow K_s\uparrow\) 推出。只有在**独立给定或测得**盐结构具有较高 \(K_s\)、因而 \(\Pi_h\ll1\) 时，下面的反馈链才成立。

已有 patch 可持续向上生长：

\[
V_{\perp}\uparrow
\]

突出结构进一步增强 tip evaporation：

\[
j_e^{tip}\uparrow
\]

并加快 vapor screening 建立：

\[
t_{scr}\downarrow
\]

从而抑制周围新成核和横向连接，最终维持 isolated patchy morphology。

与此同时，如果粗珠体系产生的初始核更稀疏，平均核间距增大也会推迟 lateral coalescence。V2 检验成核/间距通道，V4 在独立扫描 \(K_s\) 与 \(p_{c,s}\) 时检验水力通道；二者都通过后，才讨论它们与 \(d_b\) 的经验映射。

---

# 17. 一个更直接的 transition timescale criterion

定义四个时间尺度：

### 成核时间

\[
t_{nuc}
\]

### screening 建立时间

\[
t_{scr}
\]

### lateral coalescence 时间

\[
t_{coal}
\]

### vertical hydraulic limitation 时间

\[
t_h
\]

然后定义：

\[
\boxed{
\mathcal T
=
\frac{t_{scr}}{t_{coal}}
}
\]

可能存在：

## Patchy

\[
\boxed{\mathcal T<1}
\]

即：

\[
t_{scr}<t_{coal}
\]

已有 patch 在邻近 patch 合并之前就建立强 screening，因此周围新核和横向扩展被抑制。

## Crusty

\[
\boxed{\mathcal T>1}
\]

即：

\[
t_{coal}<t_{scr}
\]

盐岛在强 screening 建立之前已经发生 lateral coalescence，从而在二维形成近连续覆盖、在三维形成 spanning layer。

这个 criterion 可能比单独的 critical particle size 更具有普适性。

---

# 18. 数值几何实现：voxel / cellular growth

第一版不建议使用 phase-field。

可以采用简单 grid/voxel：

```text
0 = gas
1 = porous substrate
2 = porous/wet salt
```

candidate growth cells：

\[
\mathcal C
=
\{\text{gas cells adjacent to salt}\}
\]

对于每个 candidate face：

\[
\Delta V_i
=
V_{n,i}A_i\Delta t
\]

定义 occupancy：

\[
0\le f_i\le1
\]

更新：

\[
f_i^{n+1}
=
f_i^n
+
\frac{\Delta V_i}{V_{cell}}
\]

当：

\[
f_i\ge1
\]

该 cell 转为 salt cell。

优点：

- 易实现；
- 易做 connected-component analysis；
- 易显式改变 gas geometry；
- 很适合第一阶段验证 morphology mechanism。

**[修订] 时间步**：V2 以后的耦合主循环采用共同的**固定 \(\Delta t\)** 与占有率累积，不用“每步翻一格”的事件步，因为基底对流–扩散是瞬态方程。每次算例开始前按全程保守上界取

\[
\Delta t
\le
\min\!\left(
f_{\max}\frac{a}{V_{n,\max}},
\Delta t_{\rm CFL},
\Delta t_{\rm nuc}
\right),
\]

其中最快格填充尺度约 20 min、基底对流 CFL 尺度约 10 min，扩散项隐式；同一算例内 \(\Delta t\) 不随下一次体素翻转时刻改变。运行时只监测上述界限，若失效则该算例以更小的固定 \(\Delta t\) 重跑，而不是切回事件步。占有率 \(\ge1-\text{tol}\) 的格同时翻转，溢出量分给相邻候选格，体积守恒到舍入误差。V0 的 `GeometryEngine` 保留可选自适应步，仅用于独立几何测试与调试。
**实现状态**：`saltv2/geometry.py`（`GeometryEngine`）已实现，维数无关（2-D/3-D），`tests/test_V0_geometry.py` 8 项通过（2026-09-02）。

---

# 19. Patchy / crusty 的判据

不建议只用一个任意 coverage threshold。

计算：

## Surface coverage

\[
\theta(t)
=
\frac{A_{\rm salt,projected}}
{A_{\rm surface}}
\]

## Cluster number

\[
N_c(t)
\]

## Largest cluster fraction

\[
P_\infty
=
\frac{A_{max}}
{A_{salt}}
\]

## Height statistics

\[
H_{mean},\quad H_{rms}
\]

## 二维判据：覆盖率是否趋于 1

二维竖直截面的表面只有一条线，一维逾渗阈值为 1；只要保留一个裸露区间，就没有连续表面层。因此 spanning 在这里退化为完整覆盖，不能用 \(\theta_c<1\) 的“临界覆盖率”解释转变。

二维的操作性分类为：

\[
\boxed{
\theta(t)\text{ 在 }\theta_\infty<1\text{ 处饱和}
\Rightarrow\text{patchy}
}
\]

\[
\boxed{
\theta(t)\to1
\Rightarrow\text{crusty}
}
\]

数值上用 \(1-\theta<\varepsilon_\theta\) 识别“趋于 1”，但 \(\varepsilon_\theta\) 只是网格容差，必须做 \(a\to0\)、\(\varepsilon_\theta\to0\) 收敛检查，不能把它解释为物理阈值。Veran-Tissoires & Prat 2014 的二维格点模型采用的正是覆盖趋于连续这一思路。

## 三维判据：surface percolation

三维模型的表面是二维平面，才可定义某个盐簇跨越横向计算域的 onset，并结合 \(P_\infty\) 判断连续 crust：

\[
\boxed{\text{spanning cluster appears}\Rightarrow\text{crusty onset}}
\]

随机 continuum percolation 的 \(\theta_c\approx0.68\) 只可作为代码 sanity check 的数量级参照；盐核受输运场相关、尺寸也非独立同分布，不能把该数值直接当作模型阈值。

因此 H1–H4 先在二维用 \(\theta(t)\) 检验，形貌相图的真正逾渗版再用三维粗网格的 spanning 与 \(P_\infty\)。`saltv2/morphology.py` 已维数无关；当前 `spanning` 使用非周期切面的有限域判据，三维相图前还需用已知二维随机覆盖问题验证其有限尺寸行为。

三维只做粗网格代表点：例如 \(20\times20\times15\) mm、\(a=200\ \mu\mathrm m\) 时约有 \(7.5\times10^5\) 个未知量；按每次线性解 1–5 s、\(10^3\)–\(10^4\) 个时间步估算，单次为数小时到一天，因此不用于大矩阵扫描。

---

# 20. 每一个 timestep 的完整耦合流程

## Step 1：已知当前 salt geometry

\[
\Omega_s^n
\]

## Step 2：解 gas vapor diffusion

\[
\nabla^2 c_v=0
\]

得到：

\[
j_e(\mathbf x)
\]

## Step 3：V2 过渡闭合——计算每个 salt cluster 的总 pumping

\[
Q_k
=
\int_{\Gamma_k}j_e\,dA
\]

## Step 4：构造液相边界或联立条件

**[修订]** 足迹均匀汇只是 V2 的过渡做法；V3 起盐内 Darcy（S4）与基底 Darcy 联立，足迹上的通量由联立解自动给出，不再需要权函数 \(w(\mathbf x)\)。

bare surface：

\[
u_n=\frac{j_e}{\rho_l}
\]

salt footprint：

\[
u_n
=
\frac{Q_k}{\rho_lA_{F,k}}
\]

## Step 5：求 Darcy flow

\[
\nabla\cdot\mathbf u=0
\]

## Step 6：求 salt advection–diffusion

\[
\phi\frac{\partial C}{\partial t}
+
\nabla\cdot(C\mathbf u-\phi D\nabla C)
=-R_p
\]

## Step 7：bare surface nucleation

如果：

\[
C>C_{nuc}
\]

创建新的 salt seed。

## Step 8：计算已有 salt surface 的 hydraulic accessibility

\[
\Pi_h
\]

以及：

\[
A_h
\]

## Step 9：salt growth

\[
V_n=A_h(\Pi_h)\,G\,F(Da^*)\,j_e
\]

## Step 10：更新 voxel occupancy

\[
f_i^{n+1}
\]

## Step 11：cluster analysis

计算：

\[
N_c,
\theta,
P_\infty,
H_{mean},
H_{rms}
\]

## Step 12：按维数检测 transition

- 二维：记录 \(\theta\to1\) 的 crust formation time，并区分 \(\theta_\infty<1\) 的 patchy 饱和态；
- 三维：记录首次出现 spanning cluster 的 \(t_c\)，同时保存 \(P_\infty(t_c)\) 与 \(\theta(t_c)\)。

---

# 21. 伪代码

```text
initialize porous medium
initialize C = C0
initialize salt_geometry = empty

while t < tend:

    # 1. gas phase
    solve vapor diffusion
    compute local evaporation flux je

    # 2. liquid supply / pumping
    if version == V2:
        for each salt cluster:
            Qevap = integrate(je over salt surface)
            impose a uniform liquid sink on its substrate footprint
    else:  # V3+
        solve coupled substrate-salt Darcy problem
        enforce pressure and mass-flux continuity at the interface
        obtain footprint flux and pumping from the coupled solution

    # 3. porous substrate
    solve Darcy flow
    solve salt advection-diffusion

    # 4. nucleation
    for each bare surface cell:
        if Csurface > Cnuc:
            create salt seed

    # 5. growth
    for each exposed salt face:
        compute hydraulic number Pi_h
        Ah = hydraulic_accessibility(Pi_h)
        Vn = Ah * G * F(Da_star) * je
        add salt volume

    # 6. geometry update
    update salt geometry

    # 7. morphology metrics
    calculate:
        surface coverage
        cluster number
        nearest-neighbor distance
        mean height
        height standard deviation
        2-D: coverage convergence toward 1
        3-D: spanning cluster and P_inf

    t += dt
```

---

# 22. 推荐的参数标定顺序

原则：**不要同时 fit 所有参数。**

## Stage A：无盐 pure-water evaporation

标定：

\[
D_v,
K,
j_e
\]

目标：

- 总蒸发质量 \(m(t)\)；
- steady evaporation rate；
- gas-domain boundary condition。

---

## Stage B：只看 nucleation onset

标定：

\[
C_{nuc},
\sigma_K,
l_c
\]

目标：

- 首次析盐时间 \(t_{nuc}\)；
- 初始 nuclei 数 \(N_{nuc}\)；
- nuclei spatial distribution；
- nearest-neighbor spacing \(L_n\)。

---

## Stage C：isolated patch early growth

在 patch 尚未明显相互影响时，直接检验 Prat (2026) 闭合：

\[
V_n=G\,F(Da^*)\,j_e
\]

其中 \(\varepsilon_s,a_v,D_{s,\rm eff},k_r\) 优先取独立测量或文献值并传播不确定度，不再引入工程系数 \(\eta\) 吸收误差。

目标：

- patch height \(H(t)\)；
- lateral radius \(R(t)\)；
- total salt mass。

---

## Stage D：验证 screening / exclusion

从 simulation 后处理中得到：

\[
L_{ex}^{sim}(t)
\]

与实验或已有 scaling 比较。

如果 pumping + gas screening 足够，原则上不应再直接 fit 一个固定 \(L_{ex}\)。

---

## Stage E：最后标定 hydraulic limitation

利用：

\[
H(t)
\]

和：

\[
R_{\parallel}(t)
\]

独立约束或扫描：

\[
K_s,
\quad p_{c,s}\ \text{（或 }r_s\text{）}
\]

两者不能由同一个孔径通过 Carman–Kozeny 同时闭合。若数据不足，则把它们作为带先验范围的相图轴，不用形貌标签反标出唯一值。

---

# 23. 最有价值的 validation strategy

**[修订] regime 先行**：Ben Amara 2025 的四组实验是 **drying**（真空饱和、有限盐水）；分层基底 #2/#3 的关键过程是粗层先失饱和，需要基底两相模型，第一版（wicking）不能预测。
wicking 下 regime 匹配的转变数据是 **Veran-Tissoires & Prat 2014 第二组实验**：\(\delta=15\) mm，珠径 5–50 / 50–62 / 100–120 / 1000 µm 四档，其中 5–50 µm 结壳（Fig. 23），Table 2 给起始时间。
因此第一版的 calibration–prediction 改为：用 1000 µm（patchy）与 5–50 µm（crusty）两端约束跨算例共享的传质、成核与数值参数，冻结这些参数后预测 50–62 与 100–120 µm 两档。\(K_s\) 与 \(p_{c,s}\) 必须来自独立测量/先验范围或作为相图轴扫描，不能按已知的 patchy/crusty 标签分别回填，否则预测会变成循环论证。
Ben Amara 的数据只比早期饱和段（`plan/09` §3 判定规则），分层基底留到 drying 扩展。原稿以下内容作为 drying 阶段的验证保留：

近期 fine/coarse bead 实验非常适合作为 calibration–prediction 数据集（**drying 阶段**）。

## Calibration

只使用两个极端体系：

\[
1\text{–}50\ \mu m
\]

和：

\[
400\text{–}600\ \mu m
\]

分别校准 crusty 与 patchy 的增长和输运条件。

## Prediction

然后**不再调整参数**，直接模拟：

- fine substrate + 1 mm coarse top layer；
- fine substrate + 5 mm coarse top layer。

如果模型能够预测：

\[
\text{crusty}
\rightarrow
\text{mixed}
\rightarrow
\text{patchy}
\]

则说明模型抓住了 near-surface transport + salt-growth transition 的主要机制。

---

# 24. 第一版模型刻意不加入的机制

为了集中回答 patchy–crusty transition，第一版建议暂时不加入：

- phase-field crystal morphology；
- crystal lattice orientation；
- cubic NaCl anisotropy；
- full pore-network inside salt；
- crystallization pressure / mechanics；
- subflorescence damage；
- unsaturated drying；
- temperature coupling；
- full 3D micro-CT resolved geometry。

这些都可以在 MVP 成功后逐步加入。

---

# 25. 第一阶段最值得验证的 4 个假设

## H1：离散 nucleation 是否可以由 substrate transport heterogeneity 自然产生？

如果可以，则不需要人为随机撒核。

---

## H2：pumping + gas screening 是否能自然形成 exclusion zone？

如果可以，则：

\[
L_{ex}
\]

可以成为 emergent variable。

---

## H3：hydraulic limitation 是否能把 vertical growth 转化为 lateral growth？

即：

\[
\Pi_h\uparrow
\Rightarrow
V_{\parallel}/V_{\perp}\uparrow
\]

这是 crust formation 的关键假设。

---

## H4：patchy–crusty 是否可以由时间尺度竞争解释？

即：

\[
\mathcal T
=
\frac{t_{scr}}{t_{coal}}
\]

是否存在近似 transition：

\[
\mathcal T<1 \Rightarrow patchy
\]

\[
\mathcal T>1 \Rightarrow crusty
\]

如果成立，这会比“critical bead size”具有更强的一般性。

---

# 26. 最终希望得到的 morphology map

初步可以考虑：

\[
\boxed{
\mathcal M
=
F(Pe_p,Da_s,\Pi_h,\Lambda)
}
\]

其中：

\[
Pe_p=\frac{UL}{D_{eff}}
\]

控制 substrate 中 advection / diffusion；

\[
Da_s
\]

控制盐结构内部 precipitation/growth；

\[
\Pi_h
\]

控制 efflorescence 内部 liquid accessibility；

\[
\Lambda=\frac{L_{ex}}{L_n}
\]

描述已有 patch 对后续 nucleation/coalescence 的空间抑制。

最终 morphology 不再写成：

\[
d_p<d_c\Rightarrow crusty
\]

而是类似：

\[
(Pe_p,Da_s,\Pi_h,\Lambda)
\longrightarrow
\{patchy,\ mixed,\ crusty\}
\]

---

# 27. 可能形成的论文主线

一个比较有潜力的论文题目方向可以是：

> **Efflorescence as a nucleation–growth–percolation process: a mechanistic model for the patchy-to-crusty transition**

核心贡献可以包括：

1. 将 patchy/crusty 统一解释为 nucleation–growth–percolation transition；
2. 显式耦合 substrate salt transport、salt pumping 和 gas-phase vapor screening；
3. 引入 hydraulic accessibility 控制 vertical / lateral growth competition；
4. 让 exclusion zone 和 crust formation 成为模型自然输出；
5. 提出时间尺度 criterion：

\[
\mathcal T=\frac{t_{scr}}{t_{coal}}
\]

6. 第一阶段用 2014 wicking 的四档珠径做 calibration–prediction；coarse top-layer 数据留给后续 drying 扩展。

---

# 28. 当前最推荐的开发顺序

**[修订] 与 `saltv2/` 及子计划的对应**（详见 `plan/11` §4）：V0 已实现；V1 = 子计划 S1 + S2（气相 FV、Prat 生长律）；V2 基底域新建；V3 = S4 盐内 Darcy 与基底联立；V4 = S4 门控加软化参数 \(\alpha\)，\(\alpha\to\infty\) 回硬门控；V5 = 相图，轴取 §26。干盐蒸汽阻力、Kelvin、前沿密化（S5）与分级规则（S6）推迟到相图之后，因为转变是几何问题，蒸发抑制的 14× 是另一层。

## Version 0：几何算法测试（**已完成，2026-09-02**：`saltv2/grid.py, geometry.py, morphology.py`；8 项测试通过）

只做：

- prescribed nucleation；
- prescribed \(j_e\)；
- voxel growth；
- connected-component / percolation analysis。

确认 growth 和 morphology detection 代码工作正常。

当前八项回归覆盖：平层高度、含溢出的体积守恒、单核镜像对称、双核合并、周期边缘团簇去重、二维完整覆盖判据、高度/最近邻间距统计，以及三维生长与标记冒烟测试。

## Version 1：加入 gas diffusion

求：

\[
\nabla^2 c_v=0
\]

实现 geometry-dependent evaporation 和 vapor screening。

## Version 2：加入 substrate salt transport

求：

- Darcy flow；
- advection–diffusion；
- supersaturation-based nucleation。

让 nuclei locations 自发产生。

## Version 3：加入 pumping feedback

把盐内 Darcy（S4）与基底 Darcy 联立，在界面施加压力与质量通量连续；已有 salt cluster 对基底流场和盐供给的重定向由联立解自动产生，不再使用人为的 footprint 权函数。

## Version 4：加入 hydraulic accessibility

实现：

\[
V_n=A_h(\Pi_h)\,G\,F(Da^*)\,j_e
\]

测试 patchy ↔ crusty transition 是否出现。

## Version 5：参数扫描与 morphology map

物理输入扫描：

- 饱和基底输运：\(K_{pm}\)、相关非均匀场 \(\sigma_K,l_c\) 与 \(Pe_p\)；
- 盐微结构：\(K_s\) 与 \(p_{c,s}\) 独立扫描，再换算为 \(\Pi_h\)；
- 边界/成核：RH（或 evaporation strength）、initial salt concentration、\(S_{nuc}\)；
- substrate particle size 只用于设置有独立依据的基底参数与实验分组，不用于闭合 \(K_s(d_b)\)。

输出：

\[
\{patchy,mixed,crusty\}
\]

以及 transition map。

V0–V5 到相图的总工作量估计为 **15–17 人日**；其中 V0 已完成，余下约 14–16 人日。H1–H4 的二维机制算例先完成，三维只做粗网格 percolation 验证与少量代表点，不做大矩阵扫描。

---

# 29. 后续与 DuMuX 的衔接

如果 MVP 证明上述机制能够产生 transition，可以进一步把模型迁移到 DuMuX：

- porous substrate：现有 Darcy / compositional transport；
- gas domain：free-flow diffusion 或简化 vapor transport；
- moving salt region：独立 solid/porous domain 或 dynamic boundary layer；
- surface coupling：mass transfer + precipitation + pumping；
- morphology：先用 surface height field，后续再考虑显式 geometry evolution。

建议在迁移到 DuMuX 前，先用一个更容易快速迭代的 research code 验证机制。

---

# 30. 当前最核心的研究判断

我们目前最值得沿着下面这条路线继续：

\[
\boxed{
\text{substrate heterogeneity}
\rightarrow
\text{discrete nucleation}
}
\]

\[
\boxed{
\text{existing salt}
\rightarrow
\text{pumping + vapor screening}
\rightarrow
\text{exclusion/localization}
}
\]

\[
\boxed{
\text{salt hydraulic accessibility}
\rightarrow
\text{vertical vs lateral growth}
}
\]

\[
\boxed{
\text{coalescence/percolation}
\rightarrow
\text{crust formation}
}
\]

因此 patchy–crusty transition 可以被看作：

> **一个受多孔输运、蒸发场反馈和盐结构水力可达性共同控制的 nucleation–growth–percolation transition。**

这比单纯寻找一个“critical particle size”或经验 morphology criterion 更有机会形成具有普适性的理论框架。

---

# 31. 讨论中涉及的代表性文献方向

下面列出本轮讨论中重点涉及的文献线索，后续建议逐篇精读并整理参数表。

1. **Veran-Tissoires / Prat 系列**：NaCl efflorescence、surface nucleation、vapor screening、pumping、growth morphology。
2. **Eloukabi et al. (2013)**：fine/coarse porous media 上 crusty / patchy 转变。
3. **Roy et al. (2023)**：efflorescence localization / exclusion distance。
4. **Prat (2026)**：porous salt structure evaporation-driven growth law，局部生长与 evaporation flux 的关系。
5. **近期 fine/coarse bead structured-light scanning 实验（2025）**：crusty / patchy morphology 与 top-layer manipulation。
6. **Thin-film Volmer–Weber growth**：island nucleation–growth–coalescence。
7. **KJMA / Avrami theory**：nucleation + growth + impingement。
8. **DLA / electrodeposition**：transport screening 与 morphology instability。
9. **Reactive precipitation / clogging**：Pe–Da morphology competition。
10. **Continuum percolation**：仅在三维表面用 spanning / \(P_\infty\) 定义 continuous-layer onset；二维截面用 \(\theta\to1\)。

---

# 32. 下一步最具体的工作

建议下一步直接做两件事：

### A. 建立模型变量与方程表

把所有未知量、单位、方程、边界条件和参数整理成一张数学模型表，例如：

| 变量 | 含义 | 方程/来源 |
|---|---|---|
| \(p\) | substrate liquid pressure | Darcy |
| \(\mathbf u\) | liquid Darcy velocity | Darcy |
| \(C\) | dissolved NaCl concentration | advection–diffusion |
| \(c_v\) | vapor concentration | gas Laplace equation |
| \(j_e\) | local evaporation flux | vapor gradient |
| \(f_s\) | salt occupancy | geometry update |
| \(\Pi_h\) | hydraulic accessibility number | salt supply closure |
| \(V_n\) | salt interface growth velocity | growth law |

### B. 开发一个最小二维 prototype

先验证：

1. 单个 patch 是否会产生 vapor screening；
2. 两个 patch 是否存在 mutual screening；
3. 加入 pumping 后 exclusion 是否增强；
4. 改变 \(\Pi_h\) 后能否从 vertical patch growth 转变为 lateral coalescence；
5. 二维中 \(\theta(t)\) 是否从 \(<1\) 的饱和态转为趋近 1；随后在三维粗网格中是否出现相应的 spanning / \(P_\infty\) transition。

只要这五步能够跑通，就已经有足够基础进入系统参数扫描和实验对比。
