# 外部盐沉积模型开发：基于 Ben Amara et al. (2025) 的建模启示

## 1. 论文定位

Ben Amara et al. (2025) 研究的是**初始近饱和盐水充满多孔介质后的 drying（干燥）过程**，而不是具有外部盐水储液池持续供液的 evaporation–wicking 过程。

实验中，玻璃珠首先浸入接近饱和的 NaCl 溶液，并通过抽真空去除孔隙中的气泡。随后样品置于受控温湿度环境中，从顶部表面蒸发，并通过质量损失测量干燥动力学。

因此，其基本过程为：

\[
\text{initially saturated porous medium}
\rightarrow
\text{evaporation}
\rightarrow
\text{progressive drying}
\]

而不是：

\[
\text{reservoir}
\rightarrow
\text{continuous capillary supply}
\rightarrow
\text{evaporation}
\]

需要注意，论文中虽然反复讨论 **capillary pumping**，但这里的液体来源仍然是样品内部有限的初始盐溶液，而不是外部持续补给的 reservoir。

---

## 2. 这篇论文是否适合用于外部盐沉积模型开发？

结论是：

> **适合，而且非常有价值，但更适合作为物理图景和实验验证依据，而不是完整的数学闭合模型。**

尤其对于希望描述如下过程的模型：

\[
\text{surface nucleation}
\rightarrow
\text{external efflorescence growth}
\rightarrow
\text{液体进入盐结构}
\rightarrow
\text{顶部继续析盐}
\rightarrow
\text{反过来改变 evaporation}
\]

这篇论文能够提供非常重要的物理约束。

---

## 3. 最重要的建模启示：外部盐不能简单视为不可渗透固体

论文明确指出，patchy efflorescence 本身具有明显的多孔结构。

因此，外部析出的盐不能简单建模为：

\[
\text{impermeable solid}
\]

更合理的描述是：

\[
\boxed{
\text{external porous salt domain}
}
\]

即：

\[
\text{original porous medium}
\rightarrow
\text{porous salt structure}
\rightarrow
\text{air}
\]

盐结构内部存在连通孔隙，底部仍可与原始多孔介质中的盐溶液保持 hydraulic connection。

溶液可以通过盐结构内部的毛细作用继续向上输运，而蒸发主要发生在盐结构外表面，尤其是突出部位。

---

## 4. Active evaporation surface 会随盐析结构生长而移动

在没有外部盐沉积时，蒸发表面为原始多孔介质顶部：

\[
\Gamma_{\mathrm{evap}}=\Gamma_{\mathrm{PM}}
\]

但随着 efflorescence 生长，真正的活跃蒸发表面应逐步变成：

\[
\boxed{
\Gamma_{\mathrm{evap}}
=
\text{outer surface of efflorescence}
}
\]

因此，外部盐沉积不仅改变固体质量，还会改变：

- 蒸发表面积；
- 局部蒸发通量；
- 气相扩散路径；
- 毛细供液路径；
- 盐离子输运路径。

这意味着 geometry evolution 是整个过程的重要反馈环节。

---

## 5. 推荐的三域建模框架

一个较合理的连续介质框架可以写成：

\[
\boxed{
\Omega_{\rm PM}
+
\Omega_{\rm salt}
+
\Omega_{\rm gas}
}
\]

### 5.1 原始多孔介质区域 \(\Omega_{\rm PM}\)

主要描述：

- 水/盐水流动；
- 盐离子对流；
- 盐离子扩散；
- 饱和度变化；
- 向表面的盐输运。

基本控制方程可包括：

\[
\nabla\cdot \mathbf V=0
\]

\[
\mathbf V
=
-\frac{k}{\mu}
\left(
\nabla P_l-\rho_l\mathbf g
\right)
\]

盐输运：

\[
\frac{\partial(\rho\epsilon SC)}{\partial t}
+
\nabla\cdot
(\rho\epsilon USC)
-
\nabla\cdot
(\rho\epsilon SD_{\rm eff}\nabla C)
=0
\]

并通过 Péclet 数衡量 advection 和 diffusion 的相对作用：

\[
Pe=
\frac{U_{\rm ref}L}
{D_{\rm eff}}
\]

---

### 5.2 外部盐析层区域 \(\Omega_{\rm salt}\)

当表面达到析盐条件，例如：

\[
C\geq C_{\rm sat}
\]

或：

\[
S_{\rm supersat}
=
\frac{C}{C_{\rm eq}}
>1
\]

则形成新的外部盐析区域：

\[
\Omega_{\rm salt}
\]

该区域应具有自己的有效物性参数：

\[
\epsilon_s,\qquad
k_s,\qquad
p_{c,s},\qquad
D_{v,s}^{\rm eff}
\]

其中：

- \(\epsilon_s\)：盐析层孔隙率；
- \(k_s\)：盐析层渗透率；
- \(p_{c,s}\)：盐析层毛细压力；
- \(D_{v,s}^{\rm eff}\)：水蒸气在盐结构内部的有效扩散系数。

这篇论文的重要贡献之一，就是明确支持将盐析层看作**新的多孔介质**。

---

## 6. 盐析层内部的毛细供液

论文引用 Veran-Tissoires & Prat 的模型，认为盐析结构内部的液体输运由 capillary pumping 驱动。

可写成压力平衡：

\[
P_c
\sim
\rho gH+
\Delta P_{\rm visc}
\]

其中 viscous pressure loss 近似为：

\[
\Delta P_{\rm visc}
\sim
\frac{\mu UH}{k_s}
\]

因此：

\[
\boxed{
P_c
\sim
\rho gH+
\frac{\mu UH}{k_s}
}
\]

这里：

- \(P_c\)：盐析层内部毛细压力；
- \(H\)：盐析结构高度；
- \(k_s\)：盐析层渗透率；
- \(U\)：盐析层中的液体速度；
- \(\mu\)：液体黏度。

该关系可以作为外部盐结构液体输运模型的重要基础。

---

## 7. Crusty–patchy 转变的核心物理机制

论文实验显示：

- 细玻璃珠（约 1–50 μm）更容易形成连续 crust；
- 粗玻璃珠（约 400–600 μm）更容易形成 patchy efflorescence；
- 细颗粒基质上加入 1 mm 粗颗粒表层，会出现 crust + patchy 共存；
- 加入 5 mm 粗颗粒层后，patchy 特征更加明显。

因此，表层孔径和孔径分布是重要控制因素。

### 7.1 粗孔结构：有利于 patchy

粗颗粒对应较大的盐析层孔隙：

\[
r_s\uparrow
\]

通常带来更高渗透率：

\[
k_s\uparrow
\]

因此液体更容易持续输送到已有盐结构顶部：

\[
\text{solution}
\rightarrow
\text{capillary pumping}
\rightarrow
\text{tip}
\]

有利于：

\[
\boxed{
\text{vertical / branched growth}
}
\]

最终形成离散的 patchy efflorescence。

### 7.2 细孔结构：有利于 crust

对于较细的孔隙：

\[
r_s\downarrow
\]

虽然毛细压力会增加：

\[
P_c\sim\frac{2\gamma\cos\theta}{r_s}
\]

但渗透率大致满足：

\[
k_s\sim r_s^2
\]

因此 viscous resistance：

\[
\Delta P_{\rm visc}
\sim
\frac{\mu UH}{k_s}
\]

会随孔径减小快速增大。

于是，当盐析层孔隙非常小时，纵向供液能力受到限制：

\[
\boxed{
\text{vertical growth 被抑制}
}
\]

盐结构更倾向于：

\[
\text{lateral growth}
\rightarrow
\text{coalescence}
\rightarrow
\text{continuous crust}
\]

---

## 8. 一个有用的尺度估计

若暂时忽略重力：

\[
P_c
\sim
\frac{\mu UH}{k_s}
\]

则可以得到一个最大可持续供液高度：

\[
H_c
\sim
\frac{k_sP_c}{\mu U}
\]

进一步利用：

\[
k_s\sim r_s^2
\]

和：

\[
P_c\sim\frac{\gamma}{r_s}
\]

得到：

\[
\boxed{
H_c
\sim
\frac{\gamma r_s}{\mu U}
}
\]

因此：

\[
H_c\propto r_s
\]

即：

- 盐析层孔径越大，越容易维持纵向生长；
- 盐析层孔径越小，越容易转向横向铺展形成 crust。

这一尺度关系可以作为未来建立 crusty–patchy transition criterion 的一个出发点。

需要强调的是，这不是 Ben Amara et al. (2025) 直接给出的正式临界判据，而是基于其引用的机理模型进行的简化推导。

---

## 9. 蒸发通量同样参与控制形貌

如果：

\[
U\sim\frac{j}{\rho}
\]

则上式可写为：

\[
H_c
\sim
\frac{\gamma\rho r_s}{\mu j}
\]

因此：

\[
j\uparrow
\Rightarrow
H_c\downarrow
\]

即更强蒸发意味着更高的液体补给需求，也意味着更大的 viscous pressure loss。

所以在相同孔隙结构下：

\[
\boxed{
\text{更高 evaporation flux 可能促进 crust formation}
}
\]

这说明形貌不是单由孔径决定，而是由：

\[
r_s,\quad k_s,\quad P_c,\quad j
\]

等因素共同控制。

---

## 10. Patchy 为什么不会最终铺满整个表面？

论文给出两个关键机制：

### 10.1 Evaporation screening

随着已有盐结构向空气中突出，其尖端蒸发增强。

气相中若扩散占主导，可近似满足：

\[
\nabla^2 c_v=0
\]

蒸发通量：

\[
j=-D_v\nabla c_v\cdot n
\]

突出部位通常有：

\[
j_{\rm tip}>j_{\rm flat}
\]

因此已有盐结构会“屏蔽”周围区域，使两个 patch 之间的局部蒸发减弱，新的表面浓缩和成核因此受到抑制。

### 10.2 Ion redirection

已有 patchy efflorescence 与下方盐水保持液压连通，相当于一个毛细泵。

液体和盐离子被优先输送至已有晶体结构顶部：

\[
\text{substrate}
\rightarrow
\text{salt structure}
\rightarrow
\text{tip}
\]

于是形成正反馈：

\[
\text{existing crystal}
\rightarrow
\text{higher evaporation}
\rightarrow
\text{stronger ion supply}
\rightarrow
\text{faster growth}
\]

同时周围裸露区域获得的盐离子供应减少。

论文将这一机制概括为：

\[
\boxed{
\text{screening + ion redirection}
}
\]

这可以解释为什么 patchy morphology 能够长期保持离散，而不是最终演化成均匀 crust。

---

## 11. Laplacian growth 对 patchy 形貌的重要性

在气相扩散控制条件下：

\[
\nabla^2 c_v=0
\]

突出部位由于浓度梯度更大而蒸发更强：

\[
j_{\rm tip}>j_{\rm valley}
\]

从而：

\[
\text{tip grows faster}
\]

形成类似 diffusion-limited growth 的形貌不稳定性：

\[
\text{tip growth}
\rightarrow
\text{higher evaporation}
\rightarrow
\text{more solution supply}
\rightarrow
\text{more salt deposition}
\]

因此 patchy efflorescence 可以视为一种：

\[
\boxed{
\text{evaporation-controlled Laplacian growth}
}
\]

这为未来的 geometry evolution 或 phase-field / moving-boundary 模型提供了很好的理论方向。

---

## 12. Crust 对蒸发的反馈

对于 crusty efflorescence，论文观察到蒸发通量显著下降。

其主要原因包括：

1. crust 的有效渗透率降低；
2. 水蒸气穿过 crust 的扩散阻力增加；
3. 随 crust 逐渐变干、致密和连续，其 vapor permeability 进一步下降；
4. 很小的孔隙可能引入 Kelvin effect。

论文采用 Kelvin equation：

\[
a_w
=
\frac{P_{vs}}{P_{v\infty}}
=
\exp
\left[
-\frac{M_v}{RT}
\frac{2\gamma\cos\theta}{r}
\right]
\]

说明当孔径达到纳米尺度时，water activity 会受到明显影响。

因此，crust 一旦形成，会对后续 evaporation 产生显著抑制。

---

## 13. 这篇论文目前不能提供的关键闭合关系

虽然物理图景非常有价值，但论文没有给出真正闭合的外部盐沉积模型。

目前仍缺少：

\[
v_n=f(C,j,k_s,P_c,\ldots)
\]

即外部盐结构界面的生长速度定律。

此外也缺少：

\[
\epsilon_s(C)
\]

\[
k_s(C)
\]

\[
r_{p,s}(d_{\rm substrate})
\]

等定量 constitutive relations。

论文也没有给出明确的 crusty–patchy 临界判据，例如：

\[
r_{\rm critical}
\]

或：

\[
Pe_{\rm critical}
\]

或某个：

\[
\Pi_{\rm critical}=1
\]

因此，其主要作用应被理解为：

\[
\boxed{
\text{experiment + mechanistic picture}
}
\]

而不是：

\[
\boxed{
\text{fully predictive transition theory}
}
\]

---

## 14. 一个值得发展的无量纲 transition criterion

可以在 Veran-Tissoires–Prat 的物理图景基础上定义：

\[
\Pi
=
\frac{\text{available capillary pressure}}
{\text{pressure required to feed evaporation}}
\]

即：

\[
\boxed{
\Pi=
\frac{P_c}
{\rho gH+\mu UH/k_s}
}
\]

可以提出如下假设：

当：

\[
\Pi\gg1
\]

毛细供液足够：

\[
\boxed{
\text{vertical growth}\rightarrow\text{patchy}
}
\]

当：

\[
\Pi\lesssim1
\]

供液能力不足：

\[
\boxed{
\text{vertical growth suppressed}
\rightarrow
\text{lateral growth}
\rightarrow
\text{crust}
}
\]

这不是论文已经建立好的理论，而是非常适合作为后续模型开发的方向。

---

## 15. 对我们的模型，最核心的反馈回路

最终模型应尽量体现：

\[
\boxed{
\text{evaporation}
\rightarrow
\text{liquid flow}
\rightarrow
\text{ion transport}
\rightarrow
\text{precipitation}
\rightarrow
\text{geometry change}
\rightarrow
\text{evaporation change}
}
\]

其中 geometry change 不应只改变盐质量，还需要改变：

- 活跃蒸发表面位置；
- 有效蒸发表面积；
- 局部气相浓度梯度；
- 毛细供液路径；
- 盐离子输运路径；
- 盐层有效渗透率和扩散系数。

---

## 16. 这篇论文适合作为哪些模型部分的依据？

| 模型问题 | 适用程度 |
|---|---|
| 外部盐是否应视为 porous domain | 很适合 |
| 外部盐内是否存在液体输运 | 很适合 |
| capillary pumping 是否重要 | 很适合 |
| evaporation surface 是否随盐结构移动 | 很适合 |
| efflorescence 是否反馈 evaporation | 很适合 |
| pore size 是否影响 morphology | 很适合 |
| patchy/crusty morphology validation | 很适合 |
| salt-layer \(k_s,\epsilon_s\) | 仅提供概念和估算 |
| precipitation kinetic law | 不提供 |
| nucleation law | 不提供 |
| moving-boundary growth law | 不提供 |
| crusty–patchy 临界判据 | 尚无成熟模型 |
| 长期 wicking 条件下盐结构生长 | 不适合直接验证 |

---

## 17. 对 evaporation–wicking 模型的局限

这篇论文属于 drying。

因此：

\[
M_{\rm water}(t)\downarrow
\]

内部盐水储量有限。

后期 efflorescence 停止生长可能来自：

- hydraulic connection 丧失；
- 可蒸发水耗尽；
- 盐离子可输运量逐渐减少。

而在 evaporation–wicking 情况下：

\[
\text{reservoir}
\rightarrow
\text{continuous saline supply}
\rightarrow
\text{efflorescence}
\]

盐和水可以持续补给。

因此需要进一步回答：

\[
H_s(t)\ ?
\]

例如：

\[
H_s\rightarrow H_{\max}
\]

还是：

\[
H_s\propto t^\alpha
\]

甚至在理想条件下持续增长。

这些长期行为无法直接由这篇 drying 论文验证。

---

## 18. 推荐的模型开发路线

### V1：先建立可生长的 porous salt layer

第一阶段建议不直接追求完整 morphology prediction，而是先保证：

- water mass conservation；
- salt mass conservation；
- precipitation；
- porous salt layer 形成；
- 盐层内部液体输运；
- 盐层内部 vapor diffusion；
- evaporation surface 随盐层生长更新。

目标是先建立：

\[
\boxed{
\text{conservative external deposition model}
}
\]

### V2：加入 morphology selection

第二阶段再进一步加入：

- crusty vs patchy；
- 局部 geometry evolution；
- evaporation screening；
- ion redirection；
- Laplacian growth；
- pore-size-dependent growth direction；
- transition criterion。

可考虑建立：

\[
v_n
=
f(J_{\rm salt},C-C_{\rm sat},j_{\rm evap},k_s,P_c)
\]

并引入无量纲判据：

\[
\Pi=
\frac{P_c}
{\rho gH+\mu UH/k_s}
\]

作为 morphology selection 的候选变量。

---

## 19. 总结

这篇论文非常适合作为 external salt deposition model 的核心参考文献之一。

其最重要的建模启示可以概括为：

\[
\boxed{
d_{\rm grain}
\rightarrow
r_{\rm pore}
\rightarrow
r_{\rm salt}
\rightarrow
(k_s,P_c)
\rightarrow
\text{capillary supply}
\leftrightarrow
\text{evaporation demand}
}
\]

最终决定：

\[
\begin{cases}
\text{vertical growth} & \rightarrow \text{patchy}\\
\text{lateral growth} & \rightarrow \text{crusty}
\end{cases}
\]

一旦形成 patchy：

\[
\boxed{
\text{Laplacian evaporation}
+
\text{screening}
+
\text{ion redirection}
}
\]

会进一步稳定离散生长。

一旦形成 crust：

\[
\boxed{
k_s\downarrow
+
\text{vapor diffusion resistance}
+
\text{Kelvin effect}
}
\]

会显著抑制 evaporation。

因此，这篇论文最适合用于：

\[
\boxed{
\text{physical picture + experimental validation}
}
\]

而真正的 governing equations 和 constitutive closure，还需要结合 Veran-Tissoires & Prat (2014)、Sghaier et al. (2014) 等工作进一步发展。

---

## 参考文献

Ben Amara, M. E. A., Sghaier, N., & Kharaghani, A. (2025). *Salt crystallization on the surface of porous media: Morphological evolution of salt crystals and growth kinetics*. Drying Technology, 43(13), 1992–2007.

Veran-Tissoires, S., & Prat, M. (2014). *Evaporation of a sodium chloride solution from a saturated porous medium with efflorescence formation*. Journal of Fluid Mechanics, 749, 701–749.

Sghaier, N., Geoffroy, S., Prat, M., Eloukabi, H., & Ben Nasrallah, S. (2014). *Evaporation-driven growth of large crystallized salt structures in a porous medium*. Physical Review E, 90, 042402.
