# Patchy–Crusty Transition：统一数值模型框架

可以，而且我觉得现在已经能把这个想法收敛成一个**真正可实现、可验证、最后能够预测 patchy/crusty 的模型框架**了。

我建议不要一上来做一个特别复杂的 3D 全耦合晶体生长模型。最合理的路线是：

> **先做一个“2D/2.5D substrate + 显式盐结构 + gas diffusion”的最小模型，让 patchy/crusty 自发出现；随后再逐步把盐结构内部的多孔流动加进去。**

这里最关键的是把已有文献里的三块拼起来：

\[
\boxed{
\text{surface nucleation}
+
\text{pumping/screening feedback}
+
\text{capillary-controlled salt growth}
}
\]

2012/2014 的 Prat 工作已经证明，初始离散成核可以来自多孔介质内部速度场的随机性，而已有 efflorescence 会同时产生液体/离子 pumping 和邻近区域的 vapor screening；2012 年的模型甚至已经采用 gas-phase Laplace 方程并令盐结构局部生长速度正比于局部蒸发通量。

2026 年 Prat 又进一步给出了理论上的局部 growth law：盐结构局部生长率与局部蒸发通量线性相关，并由包含盐结构孔隙率、比表面积、反应系数和有效扩散系数的 Damköhler 数修正。

我认为我们现在缺的，其实就是把这些机制组织成一个**可以跨越 patchy–crusty transition 的统一数值模型**。

---

## 一、先把物理图像明确下来

我现在最倾向于下面这个机制。

初始时：

\[
\mathrm{brine}
\rightarrow
\text{evaporation-induced upward flow}
\rightarrow
\text{surface salt accumulation}.
\]

但是表面不是均匀的。

因为：

\[
K(\mathbf x),\quad
u(\mathbf x),\quad
j_e(\mathbf x)
\]

存在 pore-scale spatial fluctuation，所以表面首先达到析盐条件的位置不同。

2012 PRL 的 pore-network simulation 就表明，由随机孔隙结构造成的 velocity fluctuations 能够产生稳定的表面盐浓度 maxima，而这些 maxima 与实验中的初始 efflorescence locations 对应。

所以：

\[
\boxed{\text{第一阶段：形成离散 nuclei}}
\]

然后一个 nucleus 一旦长出来，它就改变整个问题。

盐本身是 porous：

\[
\text{substrate}
\rightarrow
\text{liquid enters salt structure}
\rightarrow
\text{evaporation from outer salt surface}.
\]

于是发生两个 feedback。

### Pumping

已有 patch 蒸发：

\[
J_{\rm patch}\uparrow
\]

导致下面的液体：

\[
\mathbf u\rightarrow\text{patch}.
\]

同时：

\[
\mathrm{NaCl}\rightarrow\text{patch}.
\]

### Screening

盐结构突出到空气中以后改变 vapor diffusion field：

\[
\nabla^2 c_v=0
\]

导致 patch 顶部：

\[
j_e\uparrow
\]

而 patch 周围裸露 substrate：

\[
j_e\downarrow.
\]

2012 PRL 的二维模型正是这么做的：gas 区域求准稳态纯扩散 \(\nabla^2c_v=0\)，盐结构边界和湿表面保持饱和蒸汽浓度，顶部给低蒸汽浓度，然后由

\[
J_{\rm edge}=-\nabla c_v\cdot{\bf n}
\]

算局部蒸发，盐结构按照局部 \(J_{\rm edge}\) 生长。结果自然出现 salt structure 顶部增强蒸发、底部和相邻表面受到 screening。

于是形成：

\[
\boxed{
\text{existing patch}
\rightarrow
\text{more salt supply}
\rightarrow
\text{more growth}
}
\]

同时：

\[
\boxed{
\text{existing patch}
\rightarrow
\text{neighbor screening}
\rightarrow
\text{less new nucleation}
}
\]

这就是 patchy 的自增强机制。

---

# 二、那 crusty 怎么出现？

这里才是整个模型最有意思的部分。

我觉得不能简单写：

\[
d_b<d_c
\Rightarrow
\mathrm{crusty}.
\]

真正应该模拟的是：

\[
\boxed{
\text{vertical growth}
\quad \text{vs}\quad
\text{lateral growth/coalescence}.
}
\]

Eloukabi 2013 和最新 2025 实验都支持一个重要假设：

> substrate pore size 会影响生成 efflorescence 的孔结构；当 efflorescence 孔尺度较小时，其内部 viscous resistance 增大，向上 capillary supply 会受到限制，使盐结构更倾向 lateral development，最终相互 coalesce 成 crust。

2025 年实验仍然观察到：

\[
1-50~\mu m
\rightarrow
\text{crusty}
\]

而：

\[
400-600~\mu m
\rightarrow
\text{patchy}.
\]

更重要的是，在细颗粒基质上只加 **1–5 mm 粗颗粒表层**，就能够推动 crusty→patchy 的转变，这说明真正控制 morphology 的很可能是**近表面的输运和盐结构生长条件**，而不是整个样品的宏观性质。

---

# 三、我建议的计算域

第一版不要做整个 crystal 微结构。

做三个 domain：

\[
\boxed{
\Omega_p+\Omega_s(t)+\Omega_g
}
\]

分别是：

- \(\Omega_p\)：原始 porous substrate；
- \(\Omega_s(t)\)：随时间增长的 porous salt；
- \(\Omega_g\)：外部气相。

二维截面就够测试理论：

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

以后再改成：

\[
x-y
\]

二维 substrate surface 加：

\[
h(x,y,t)
\]

的 2.5D 模型。

---

# 四、substrate 中求什么？

对于你现在 evaporation-wicking 的实验，我甚至建议第一版直接假设 substrate 始终 saturated。

那么不需要 Richards equation。

只需要：

## 1. Darcy flow

\[
\nabla\cdot\mathbf u=0
\]

\[
\mathbf u
=
-\frac{K(\mathbf x)}{\mu}
\left(
\nabla p-\rho_l\mathbf g
\right).
\]

关键是：

\[
K(\mathbf x)
\]

不能完全均匀。

因为如果你设：

\[
K=\mathrm{constant}
\]

同时：

\[
j_e=\mathrm{uniform},
\]

那么所有 surface cells 的盐浓度几乎同时达到 saturation，你必然直接得到 uniform crust。

这反而把最核心的物理机制消掉了。

---

# 五、怎么产生 realistic 的初始 nuclei？

这里我强烈建议**不要直接随机撒 nuclei**。

第一版本可以生成一个 correlated permeability field：

\[
K(\mathbf{x})
=
K_0\exp[\sigma_K\xi(\mathbf{x})]
\]

其中：

\[
\langle\xi\rangle=0
\]

并给一个 correlation length：

\[
l_c\sim d_p
\]

或者：

\[
l_c\sim d_b.
\]

这样：

\[
K(\mathbf x)
\rightarrow
u(\mathbf x)
\rightarrow
c_s(\mathbf x)
\]

自然产生 local maxima。

这实际上是在 continuum framework 中模拟 PRL 2012 的 pore-network disorder。

---

# 六、盐输运

substrate 中：

\[
\boxed{
\phi
\frac{\partial C}{\partial t}
+
\nabla\cdot
\left(
C\mathbf u
-
\phi D_{\mathrm{eff}}\nabla C
\right)
=
-R_{\mathrm{prec}}
}
\]

其中：

\[
D_{\mathrm{eff}}
=
D_s f(\phi,\tau).
\]

底部 reservoir：

\[
C=C_0.
\]

---

# 七、结晶之前的 surface boundary condition 很重要

NaCl 不挥发。

所以水离开：

\[
j_e>0
\]

但是盐不能离开液相。

因此在 bare surface 上：

\[
\boxed{
\left(
C\mathbf u-
\phi D_{\mathrm{eff}}\nabla C
\right)\cdot\mathbf n=0
}
\]

这正是 evaporation-induced concentration buildup 的来源。

2012 PRL 的 pore-network model 在 crystallization onset 前同样采用了 surface meniscus 上的 zero salt flux condition。

于是：

\[
C_s(\mathbf x,t)
\]

逐渐升高。

---

# 八、nucleation criterion

最简单先用：

\[
C_s(\mathbf x,t)
>
C_{\rm nuc}.
\]

写成 supersaturation：

\[
S(\mathbf x,t)
=
\frac{C}{C_{\rm sat}},
\]

如果：

\[
S>S_{\rm nuc}
\]

则允许 nucleation。

第一版甚至可以：

\[
S_{\rm nuc}=1.
\]

因为经典 Prat 的 NaCl 实验近似把 supersaturation effect 忽略掉；但这点之后最好作为待校准量，而不是永远固定。

也就是说我们得到：

\[
\boxed{
\chi_s(\mathbf x)=
\begin{cases}
1,&C>C_{\rm nuc}\\
0,&\text{otherwise}.
\end{cases}
}
\]

---

# 九、但是仅有这个还不会得到 patchy

关键来了。

假设：

```text
surface

x x x x X x x x x
```

中间 \(X\) 首先析盐。

如果之后你仍然保持原来的：

\[
j_e(x)
\]

和：

\[
u(x),
\]

那么周围：

```text
x x X X X x x
```

很快全部达到 saturation。

最后还是 crust。

所以一定要加入：

\[
\boxed{\text{salt → transport feedback}}
\]

---

# 十、先显式求 gas vapor diffusion

这一步非常值得做，因为计算成本其实不高。

gas domain：

\[
\boxed{
\nabla\cdot
(D_v\nabla c_v)=0
}
\]

如果 \(D_v\) constant：

\[
\nabla^2c_v=0.
\]

远场：

\[
c_v=c_{v,\infty}.
\]

所有 wet surfaces：

\[
c_v
=
a_w(C)
c_{v,\mathrm{sat}}(T).
\]

包括：

- bare wet substrate；
- wet efflorescence outer surface。

然后：

\[
\boxed{
j_e
=
-D_v
\nabla c_v\cdot\mathbf n
}
\]

得到真实的 local evaporation flux。

这一步实际上就是直接继承 Veran-Tissoires et al. 2012 的方法。

---

# 十一、这样 screening 会自己出现

假设有两个 patch：

```text
        ███         ███
       █████       █████
─────────────x x──────────
```

解 Laplace equation 后自然会出现：

顶部：

\[
j_e^{\rm top}\uparrow
\]

而中间：

\[
j_e^{\rm substrate}\downarrow.
\]

所以：

\[
u_{\rm upward}
=
\frac{j_e}{\rho_l}
\]

下降。

进一步：

\[
C_{\rm surface}
\]

增长变慢。

于是：

\[
P_{\rm nuc}\downarrow.
\]

**这就是 exclusion zone。**

所以我会把：

\[
L_{\rm ex}
\]

从“输入参数”降级成“diagnostic quantity”。

最后从 simulation 中测：

\[
L_{\rm ex}^{sim}
\]

然后跟 Roy et al. 2023 的实验/scaling 比较。

---

# 十二、第二个 feedback：pumping

仅仅 gas screening 还不完整。

已有 salt patch 的 outer surface 会蒸发：

\[
Q_k
=
\int_{\Gamma_k}
j_e\,dA.
\]

那么该 patch 必须从 substrate 吸取：

\[
Q_{l,k}
=
\frac{Q_k}{\rho_l}.
\]

所以在它与 substrate 的 contact footprint \(F_k\) 上：

\[
\boxed{
\mathbf u\cdot\mathbf n
=
\frac{Q_{l,k}}{|F_k|}
}
\]

第一版可以均匀分布。

更好一点：

\[
u_n(\mathbf x)
=
\frac{Q_{l,k}}{|F_k|}
w(\mathbf x)
\]

其中：

\[
\int_{F_k}w\,dA=1.
\]

这就产生：

\[
\mathbf u
\rightarrow
\text{existing patch}.
\]

于是盐也会：

\[
C\mathbf u
\rightarrow
\text{existing patch}.
\]

这就是 PRL 中所谓的 **pumping effect**。

---

# 十三、这样 exclusion 就有两个来源

现在附近的裸表面同时经历：

### vapor screening

\[
j_e^{bare}\downarrow
\]

以及：

### flow redirection

\[
u^{bare}\downarrow.
\]

所以：

\[
\frac{\partial C^{bare}}{\partial t}\downarrow.
\]

最终：

\[
C^{bare}<C_{\rm nuc}.
\]

于是自然得到：

```text
       SALT

         █
        ███
       █████
───────xxxxx────────
       exclusion
```

这个机制比直接写：

\[
r<L_{\rm ex}\Rightarrow P_{\rm nuc}=0
\]

高级很多，而且可预测。

---

# 十四、接下来真正决定 patchy/crusty 的，是盐怎么长

这里可以直接利用 Prat 2026：

\[
\boxed{
V_n\propto j_e
}
\]

并由盐结构内部的 Damköhler number 修正。

第一版我们甚至没必要立即复现他的完整 analytical expression。

可以先写一个质量守恒形式：

\[
J_{\rm salt}
=
\eta(Da)
C_s
\frac{j_e}{\rho_l}.
\]

这里：

\[
\frac{j_e}{\rho_l}
\]

是供液速度，

所以：

\[
C_s\frac{j_e}{\rho_l}
\]

是盐 mass flux。

然后 solid volume growth：

\[
V_n
=
\frac{J_{\rm salt}}
{\rho_{\rm NaCl}(1-\varepsilon_s)}.
\]

即：

\[
\boxed{
V_n
=
\eta(Da)
\frac{C_s}
{\rho_l\rho_{\rm NaCl}(1-\varepsilon_s)}
j_e
}
\]

注意：**这是我建议第一版采用的质量守恒 closure，不是我声称它等同于 Prat 2026 的原始公式。**

后面拿到完整 equation 后直接替换：

\[
V_n=F(Da)j_e.
\]

这样模型结构完全不用变。

---

# 十五、用 voxel/cellular growth 实现最简单

比如 grid：

```text
0 = gas
1 = porous substrate
2 = liquid-filled salt
3 = solid/porous salt cell
```

盐的 candidate cells 是：

\[
\mathcal C=
\{\text{gas cells adjacent to salt}\}.
\]

对于每个 candidate face：

\[
\Delta V_i
=
V_{n,i}A_i\Delta t.
\]

定义 occupancy：

\[
0\le f_i\le1.
\]

更新：

\[
f_i^{n+1}
=
f_i^n+
\frac{\Delta V_i}{V_{\rm cell}}.
\]

当：

\[
f_i\ge1
\]

该 cell 转成 salt cell。

这和 phase-field 相比容易实现很多，而且特别适合我们现在验证 morphology mechanism。

---

# 十六、真正核心：加入 salt 内部 hydraulic limitation

如果只采用：

\[
V_n\propto j_e
\]

模型很可能会非常倾向向上长：

```text
      ██
      ██
      ██
──────██─────
```

也就是 patchy。

怎么得到：

```text
██████████████
──────────────
```

这种 crust？

我认为必须加入：

\[
\boxed{
\text{efflorescence hydraulic accessibility}.
}
\]

盐结构看作一个 porous medium：

\[
\mathbf u_s
=
-\frac{K_s}{\mu}
(\nabla p_s-\rho_l\mathbf g).
\]

内部 capillary suction：

\[
p_{c,s}
\sim
\frac{2\gamma\cos\theta}{r_s}.
\]

而把液体输送到高度 \(h\) 所需要的压力：

\[
\Delta p_{\rm req}
=
\Delta p_{\rm visc}
+
\rho_lgh.
\]

其中：

\[
\Delta p_{\rm visc}
\sim
\frac{\mu q_s h}{K_s}.
\]

定义一个非常重要的新无量纲数：

\[
\boxed{
\Pi_h
=
\frac{\Delta p_{\rm req}}
{\Delta p_{\rm cap}}
}
\]

其中：

\[
\Delta p_{\rm cap}
=
p_{c,s}-p_{c,pm}.
\]

---

# 十七、这个 \(\Pi_h\) 很可能就是 morphology 的关键之一

如果：

\[
\Pi_h\ll1
\]

说明 liquid supply 很容易到达 patch top。

于是：

\[
j_e^{top}\text{ high}
\]

\[
V_n^{top}\text{ high}.
\]

得到：

\[
\boxed{\text{vertical growth → patchy}}
\]

---

如果：

\[
\Pi_h\rightarrow1
\]

则 patch 顶部 supply 不够。

可以定义 wetness：

\[
A_h
=
\frac{1}
{1+\exp[\alpha(\Pi_h-1)]}.
\]

然后：

\[
V_n
=
A_hF(Da)j_e.
\]

于是顶部：

\[
A_h\rightarrow0.
\]

vertical growth 被压制。

但靠近 substrate 的 lateral surface：

\[
h\approx0
\]

所以：

\[
\Pi_h\ll1.
\]

仍然可以长。

因此：

```text
initial

       ██
───────██──────

later

     ██████
─────██████─────

later

 ██████████████
────────────────
```

最终邻近 islands：

\[
\boxed{\text{coalescence → crust}}.
\]

这就把 Eloukabi/Veran-Tissoires 的 capillary-rise hypothesis 真正嵌入 growth model 里面了，而不是只用它做事后解释。

---

# 十八、这样 patchy/crusty 就不用人为定义

让模型自己长。

然后做 connected-component analysis。

表面 salt coverage：

\[
\boxed{
\theta(t)
=
\frac{A_{\rm salt,projected}}
{A_{\rm surface}}
}
\]

cluster 数：

\[
N_c(t)
\]

largest connected cluster：

\[
A_{\max}(t)
\]

定义：

\[
P_\infty
=
\frac{A_{\max}}
{A_{\rm salt}}.
\]

---

# 十九、crust 的最好判据不是一个任意 coverage

我建议：

\[
\boxed{
\text{spanning cluster criterion}
}
\]

如果某一个 salt cluster：

\[
\text{left}\leftrightarrow\text{right}
\]

或者在二维表面：

\[
x\text{ direction or }y\text{ direction spans system},
\]

就认为发生 surface percolation。

随机重叠圆盘的经典 continuum percolation threshold 是大约：

\[
\theta_c\approx0.676.
\]

但是**我不会直接给盐使用 0.676**。

因为我们的 nuclei：

- 不随机；
- 存在 exclusion；
- patch 形状不是圆；
- growth strongly correlated。

所以 0.676 只能作为 sanity check。

真正 criterion 是：

\[
\boxed{\text{出现 spanning connected salt cluster}}
\]

而不是：

\[
\theta>0.676.
\]

---

# 二十、我觉得甚至可以得到一个非常漂亮的 transition theory

这里是我认为这项工作的核心潜力。

其实 patchy/crusty 可以转化成几个时间尺度的竞争。

定义：

### 新核形成时间

\[
t_{\rm nuc}
\]

### 一个 patch 建立明显 screening 的时间

\[
t_{\rm scr}
\]

### 邻近 patches lateral coalescence 的时间

\[
t_{\rm coal}
\]

### vertical growth 达到 hydraulic limitation 的时间

\[
t_h.
\]

那么：

## Patchy

如果：

\[
\boxed{
t_{\rm scr}<t_{\rm coal}
}
\]

一个 patch 在碰到邻居以前就已经把周围区域 screen 掉。

于是：

\[
\text{new nucleation suppressed}
\]

\[
\Downarrow
\]

isolated salt structures survive。

---

## Crusty

相反，如果：

\[
\boxed{
t_{\rm coal}<t_{\rm scr}
}
\]

surface islands 在 strong screening 建立之前就已经 lateral connect。

于是：

\[
\boxed{\text{crust}}.
\]

---

# 二十一、为什么 fine beads 可能 crusty？

这时候整个机制链就完整了：

\[
d_b\downarrow
\]

可能导致：

\[
r_s\downarrow,
\qquad
K_s\downarrow
\]

而通常：

\[
K_s\sim r_s^2.
\]

所以虽然：

\[
p_c\sim\frac1{r_s}
\]

增大，但 hydraulic resistance：

\[
\frac{1}{K_s}\sim\frac1{r_s^2}
\]

增长更快。

于是：

\[
\Pi_h\uparrow.
\]

vertical growth 更容易受到阻碍。

所以：

\[
V_\parallel/V_\perp\uparrow
\]

导致：

\[
t_{\rm coal}\downarrow.
\]

同时小颗粒实验中 efflorescence 确实出现得更早、更容易在表面 spread。

因此：

\[
\boxed{
d_b\downarrow
\rightarrow
\text{dense nucleation}
+
\text{lateral growth}
\rightarrow
t_{\rm coal}<t_{\rm scr}
\rightarrow
crusty
}
\]

这是一个非常漂亮、可验证的假说。

---

# 二十二、coarse beads 则反过来

\[
d_b\uparrow
\]

\[
K_s\uparrow
\]

vertical hydraulic transport 容易：

\[
\Pi_h\downarrow
\]

所以：

\[
V_\perp\uparrow.
\]

已有 patch 更快突出 surface。

突出之后：

\[
j_e^{tip}\uparrow
\]

gas screening 迅速增强。

于是：

\[
t_{\rm scr}\downarrow.
\]

最终：

\[
t_{\rm scr}<t_{\rm coal}.
\]

因此：

\[
\boxed{\text{patchy}}.
\]

这与 PRL 中“patch 越长大，顶部 evaporation 越强，同时相邻裸表面 flux 越被 screen”的机制直接吻合。

---

# 二十三、这样我们真正应该画的 morphology map 不是单纯 \(Pe-Da\)

我上一条说 \(Pe-Da\)，现在我觉得可以进一步升级。

最核心可能是三个量：

\[
\boxed{
Pe_p,\qquad Da_s,\qquad\Pi_h
}
\]

再加一个 emergent quantity：

\[
\boxed{
\Lambda
=
\frac{L_{\rm ex}}
{L_n}
}
\]

其中 \(L_n\) 是 nuclei spacing。

---

可以最终得到：

\[
\boxed{
\mathcal M
=
F(
Pe_p,
Da_s,
\Pi_h,
\Lambda
)
}
\]

然后：

```text
               Πh
               ↑

        crusty │ mixed
               │
               │
───────────────┼──────────→ Λ
               │
        crusty │ patchy
```

具体边界让 simulation 给。

---

# 二十四、我建议第一版暂时不要显式求 salt 内部 Darcy flow

否则一下子太复杂。

先定义：

\[
\Pi_h
=
\frac{\mu j_e h_s}
{\rho_l K_s\Delta p_{\rm cap}}
\]

然后：

\[
A_h
=
\frac{1}
{1+\exp[\alpha(\Pi_h-1)]}.
\]

盐生长：

\[
\boxed{
V_n
=
A_h
F(Da)
j_e
}
\]

这已经足够测试：

> hydraulic limitation + screening + pumping 是否能够产生 crusty/patchy transition。

如果可以，再把 \(A_h\) 替换成真正求解：

\[
\nabla\cdot\mathbf u_s=0
\]

和：

\[
\mathbf u_s=-K_s\nabla p/\mu.
\]

---

# 二十五、完整的 timestep 实际上非常清楚

每一个：

\[
t^n\rightarrow t^{n+1}
\]

按照下面顺序。

### ① 更新 salt geometry

已知：

\[
\Omega_s^n.
\]

### ② 求 gas vapor diffusion

\[
\nabla^2c_v=0.
\]

得到：

\[
j_e(\mathbf x).
\]

### ③ 计算每个 salt cluster 总 pumping

\[
Q_k
=
\int_{\Gamma_k}j_e\,dA.
\]

### ④ 把 evaporation/pumping 转成 substrate boundary flux

裸表面：

\[
u_n=\frac{j_e}{\rho_l}.
\]

salt footprint：

\[
u_n
=
\frac{Q_k}{\rho_lA_{F,k}}.
\]

### ⑤ 求 Darcy flow

\[
\nabla\cdot\mathbf u=0.
\]

### ⑥ 求 salt advection–diffusion

\[
\phi
\frac{\partial C}{\partial t}
+
\nabla\cdot
(C\mathbf u-\phi D\nabla C)
=
-R_p.
\]

### ⑦ 检查 bare surface nucleation

如果：

\[
C>C_{\rm nuc}
\]

创建：

\[
\text{new salt seed}.
\]

### ⑧ 对每个 existing salt surface 算 hydraulic accessibility

\[
\Pi_h
\]

以及：

\[
A_h.
\]

### ⑨ 生长 salt

\[
V_n
=
A_h F(Da)j_e.
\]

### ⑩ 更新 voxel occupancy

\[
f_i^{n+1}.
\]

### ⑪ 做 cluster analysis

输出：

\[
N_c,
\quad
\theta,
\quad
P_\infty,
\quad
H_{\rm mean},
\quad
H_{\rm rms}.
\]

### ⑫ 判断是否 percolate

如果存在 spanning cluster：

\[
\boxed{\text{crusty onset}}
\]

并记录：

\[
t_c.
\]

---

# 二十六、伪代码基本就是这样

```text
initialize porous medium
initialize C = C0
salt_geometry = empty

while t < tend:

    # gas phase
    solve vapor diffusion
    compute local evaporation flux je

    # salt pumping
    for each salt cluster:
        Qevap = integrate(je over salt surface)
        impose corresponding liquid sink on footprint

    # porous substrate
    solve Darcy flow
    solve salt advection-diffusion

    # nucleation
    for each bare surface cell:
        if Csurface > Cnuc:
            create salt seed

    # growth
    for each exposed salt face:

        compute hydraulic number Pi_h
        Ah = hydraulic_accessibility(Pi_h)

        Vn = Ah * growth_law(Da, je)

        add salt volume

    update geometry

    # morphology
    calculate:
        surface coverage
        cluster number
        nearest-neighbor distance
        mean height
        height standard deviation
        spanning cluster

    t += dt
```

这其实已经完全可以编码了。

---

# 二十七、实验参数怎么标定——这是最重要的设计

我不建议所有参数一起 fit。

要**分阶段校准**。

## Stage A：没有盐

用 pure-water experiment。

标定：

\[
D_v,
\quad
K,
\quad
j_e.
\]

首先把：

\[
m(t)
\]

搞对。

---

## Stage B：只看 nucleation onset

看：

\[
t_{\rm nuc}
\]

和第一批 salt spot locations。

拟合：

\[
C_{\rm nuc},
\quad
\sigma_K,
\quad
l_c.
\]

目标包括：

\[
N_{\rm nuc}
\]

和：

\[
L_n.
\]

---

## Stage C：只看 isolated patch 的 early growth

在 patch 尚未互相影响的时候：

\[
H(t),
\quad
R(t)
\]

拟合：

\[
F(Da)
\]

或者简化的：

\[
\eta.
\]

Prat 2026 的 growth law 就特别适合这一步，因为他的核心结论是 local growth rate 与 local evaporation flux 线性相关。

---

## Stage D：看 exclusion zone

测：

\[
L_{\rm ex}(t).
\]

对比：

- simulation；
- Roy 2023；
- PRL 2012。

如果模型正确：

\[
L_{\rm ex}
\]

应该不需要直接 fit，而能够从 pumping + screening emergently 得到。

---

## Stage E：最后才调 hydraulic limitation

用：

\[
H(t)
\]

和：

\[
R_{\parallel}(t)
\]

反演：

\[
K_s,
\quad
r_s
\]

或者：

\[
\Pi_h.
\]

---

# 二十八、2025 那篇论文其实给了我们一个极好的 validation case

它有：

\[
1-50~\mu m
\]

以及：

\[
400-600~\mu m
\]

两个极端。

还有：

- 纯 fine beads；
- fine + 1 mm coarse surface layer；
- fine + 5 mm coarse surface layer。

而且用了 structured-light surface scanning，所以不仅有 top-view morphology，还有三维高度演化。

这意味着一个很强的 validation strategy：

### Calibration

只使用：

\[
1-50~\mu m
\]

和：

\[
400-600~\mu m.
\]

### Prediction

**不再调参数**，直接模拟：

\[
1~\mathrm{mm}
\]

和：

\[
5~\mathrm{mm}
\]

coarse surface layer。

如果模型能够预测：

\[
\text{crusty}
\rightarrow
\text{mixed}
\rightarrow
\text{patchy},
\]

这会是非常有说服力的结果。

---

# 二十九、如果我们最终想形成论文，我会把核心科学问题定义成这个

不是：

> Why do small beads produce crusty efflorescence?

而是：

> **What determines whether initially discrete salt nuclei remain isolated or coalesce into a continuous crust?**

然后提出：

\[
\boxed{
\text{patchy/crusty transition is a competition between}
}
\]

\[
\boxed{
\text{nucleation}
+
\text{hydraulic-limited lateral growth}
+
\text{evaporation screening}
+
\text{salt-solution pumping}.
}
\]

更进一步，用 timescale 表达：

\[
\boxed{
\mathcal T
=
\frac{t_{\rm screen}}
{t_{\rm coalescence}}
}
\]

那么可能：

\[
\mathcal T<1
\Rightarrow
\text{patchy}
\]

而：

\[
\mathcal T>1
\Rightarrow
\text{crusty}.
\]

我觉得这个指标甚至比单纯 \(Pe-Da\) 更有潜力，因为它直接描述了**为什么 island 最终保持孤立，还是连接成 crust**。

---

## 我建议我们现在的开发顺序

第一版我会**故意保持很简单**：

1. **2D substrate Darcy + salt transport**；
2. **2D gas Laplace vapor diffusion**；
3. surface nucleation；
4. voxel salt growth；
5. \(V_n\propto j_e\)；
6. pumping feedback；
7. simplified hydraulic factor \(A_h(\Pi_h)\)；
8. connected-component / percolation analysis。

暂时**不要**加：

- phase field；
- crystal lattice orientation；
- cubic crystal anisotropy；
- full salt pore network；
- mechanical stress；
- subflorescence；
- unsaturated drying。

因为这些都不是回答 **patchy–crusty transition** 所必须的。

如果这个最小模型就能从同一套方程，在改变 \(d_b\)、\(K\)、\(K_s\) 后分别出现：

\[
\boxed{\text{isolated vertical patches}}
\]

和

\[
\boxed{\text{laterally coalesced crust}},
\]

那么我们就已经抓住了真正值得深入研究的机制。之后再把它嵌入 DuMuX 或你目前的 salt-growth framework，会比从一开始做全耦合容易很多。
