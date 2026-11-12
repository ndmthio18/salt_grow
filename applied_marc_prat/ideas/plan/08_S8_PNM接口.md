# S8 PNM 接口冻结

**Regime：接口按 wicking 设计，但字段留 drying 钩子。** V1 传常数 p_l,0、C_0；drying 阶段 PNM 传 p_l,0(x,t)、C_0(x,t) 与失饱和标记，本模型方程不变。

## 1. 目的

把 F§6 的基底边界表变成代码契约，使本模型可作为 PNM 的表面模块调用。不加物理。

## 2. 依据

F§6：基底 z = 0 上两个输入场 p_l,0(x)、C_0(x)，两个输出场 j(x)（裸露格）、U_dem(x)（被盖格）。
反馈含义：(a) 盐结构改变裸露基底的 j（屏蔽）；(b) 被覆盖孔从蒸发边界变为抽液边界；(c) 抽液量含填孔项。

## 3. 前置

S4 通过即可开始接口设计；S7 通过后冻结。

## 4. 实现

- `saltv2/interface.py`：
  - `SubstrateBC(p_l0: array[surface], C0: array[surface], saturated: bool[surface])` 输入；V1 常数且 `saturated` 全真；drying 阶段 PNM 每步更新三者（regime 钩子）。数组按表面体素索引，2D/3D 通用。
  - `SurfaceExchange(j: array[nx], U_dem: array[nx], covered: bool[nx], state_bottom)` 输出。
  - `SaltModel.step_with_bc(bc) -> SurfaceExchange`：一个事件步。
- 调用顺序文档：PNM 步（给 p_l,0）→ 本模型事件步（返回 j、U_dem）→ PNM 把 j 作蒸发边界、U_dem 作抽液边界 → 下一步。
  时间步由本模型的事件步决定，PNM 侧按其 Δt 子循环或插值，写明。

## 5. 验收（`tests/test_S8_interface.py`）

| 测试 | 判据 |
|---|---|
| mock PNM | p_l,0(x) 取随 x 线性变化的场，干涸首先出现在吸力最大处 |
| 守恒 | Σ j·a + Σ U_dem ρ_l a = 本模型账本的总蒸发 + 总抽液 |
| 常数回退 | p_l,0 常数时结果与 S4–S7 逐格一致 |
| 覆盖切换 | 某格被盖后其输出从 j 切到 U_dem，且 U_dem 含填孔项 |

## 6. 完成后

可以声称：本模型可作为表面模块嵌入任何给出表面孔液压的 PNM。
不能声称：基底内部任何东西；成核（V2）。

## 7. 工作量

1–2 人日。
