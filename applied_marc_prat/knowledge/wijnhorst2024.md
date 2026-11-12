# Wijnhorst, Van der Sloot, Pel, Shahidzadeh (2024) —— 蒸发面积对 efflorescence / subflorescence 的影响；盐壳自抬升（盐腿）

*Effect of evaporative surface area on salt efflorescence and subflorescence formation in a given porous material.*
Phys. Rev. Applied 21, 064055. 文件：`literature/Wijnhorst et al - 2024 - Effect of evaporative surface area….pdf`

**regime：drying**（截锥形 Prague 砂岩，底与侧面 Parafilm 封闭，只从顶面干燥；有限盐水，无补给）。

## 实验条件

| 项 | 值 |
|---|---|
| 样品 | 截锥砂岩，长 5 cm，孔隙率 30 %，平均孔径 30 µm；大底 / 小底（A_0 = 0.79 cm²）朝上两种 |
| 溶液 | 10 wt % NaCl（S ≈ 0.3） |
| 环境 | RH 45 ± 5 %，21 °C，氮气吹扫 |
| 测量 | 质量、顶面图像、NMR 水与 Na 剖面、SEM、轮廓仪；Pe = UL/D*，L = 5 cm，D* = 1.5e−5 cm²/s |

## 主要观测

- **大底朝上**（初始 Pe = 1.5，随干燥降到 <1）：6 h 出现 fairy ring，向心蔓延；efflorescence 是**覆盖颗粒的薄膜，孔口保持开放、不阻塞**，恒速期持续到 S = 0.14；60 h 后盐霜变白（干）、不再演化；壳与石头强粘附；57 % 的盐留在内部成 subflorescence。
- **小底朝上**（初始 Pe = 2.7，全程 >1.5）：30 min 出现 ring；cauliflower 状簇很早出现并长大；约 40 h（S = 0.8）簇盖满整个小底 → 阻塞 → 干燥速率指数衰减；31 % subflorescence。
- **盐壳自抬升**：50 h 后壳顶不再有结晶（顶面失去水力连通），结晶转到壳下方砂岩表面的液膜里，形成**盐柱/盐腿**，靠结晶压把整块壳顶起，壳与基底粘附极弱、可整块揭下（Fig 5、6c、7、8）。
- 干燥过程中观察到表面析出晶体的**湿/干循环**（Video 1）：efflorescence 本身是孔径小于石头的多孔介质，靠毛细泵吸把溶液抽到外层。
- 顶面浓度经 NMR 证实在早期即达 6 M（溶解度）峰。

## 对模型的用法

| 结论 | 用途 | 限制 |
|---|---|---|
| 同一材料、同一环境，Pe 1.5 → 非阻塞薄膜，Pe 2.7 → cauliflower 合并成阻塞壳 | **阻塞型 crust 可由 cauliflower 簇合并 + 顶部干涸产生**，不必是纳米孔薄膜；支持 T 文档的成核密度/合并通道与 S4 的门控 | drying，Pe 随时间变 |
| 顶部失去水力连通后结晶转到壳底、盐腿抬壳 | 审查 §6 的"失连不等于横向铺展"有实证；对应本项目推迟的 CONTACT/GAP 状态（记忆 `salt-v2-no-detachment-state`），机制是壳下液膜结晶压 | drying、高 Pe、基底表面失饱和；wicking 里基底保持连通，未见报道 |
| 表面晶体湿/干循环 | 审查 §7 的回湿有实证；模型不应默认"干格永久干" | drying |
| efflorescence 孔径小于基底、可毛细泵吸 | 与 2013/2014/2025/2026 一致 | — |

**注**：Prat 2026 引用本文（ref 20）作为 d_c/d_sub ≈ 0.14 的出处，但本文正文未找到该比值的明确表述，可能来自 SEM 图像；待核。
