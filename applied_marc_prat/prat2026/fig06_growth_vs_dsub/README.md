# Figure 6 — 生长速度与结构内 Darcy 速度

`../.venv/bin/python run.py` → `fig06.png`

## 对应论文内容

主图，式(32)：

$$\frac{dh/dt}{U_{evap}}=\frac{\rho C_{eq}}{\rho_{cr}(1-\varepsilon_s)(1-C_{eq})}
\cdot\frac{\sqrt{Da^*}}{1+\sqrt{Da^*}}$$

inset，式(33)：

$$\frac{U}{U_{evap}}=1+\frac{(\rho\varepsilon_s+\rho_{cr}(1-\varepsilon_s))C_{eq}}
{\rho_{cr}(1-\varepsilon_s)(1-C_{eq})}\cdot\frac{\sqrt{Da^*}}{1+\sqrt{Da^*}}$$

## 复现结果

| $d_{sub}$ | $(dh/dt)/U_{evap}$ | $U/U_{evap}$ |
|---|---|---|
| 1 µm | 0.1221 | 1.178 |
| 10 µm | 0.2196 | 1.321 |
| 25 µm | 0.2540 | 1.371 |
| 50 µm | 0.2759 | 1.403 |
| 100 µm | 0.2937 | 1.429 |

与论文一致：比值落在 0.1–0.3 区间、随 $d_{sub}$ 单调增（因 $Da^*\propto d_{sub}$），
inset 的 $U/U_{evap}$ 始终 > 1 且在 1.15–1.45 之间。

## 这张图承担的特殊职责：标定 ρ

论文**从未给出溶液密度 ρ**，但式(32) 需要它。本项目把 ρ = 1200 kg/m³
（饱和 NaCl 盐水）作为反标定值，其**全部依据就是这张图的两个点**：

- $d_{sub}$ = 25 µm → 我方 **0.254**，论文正文写 "0.25 cm/day for $U_{evap}$ = 1 cm/day"
- $d_{sub}$ = 10 µm → 我方 **0.2196**，论文 Fig 6 读数 ≈ 0.21–0.22

两个不同粒径同时对上，说明这不是单点凑合。`run.py` 末尾有一条 assert 守住它，
`tests/test_anchors.py` 里有对应的两条测试。

**任何改动 `params.RHO` 的操作都必须重跑这两条检查。**

## 为什么 U 一定大于 U_evap

由式(6) 的界面质量平衡：向上输送的溶液不仅要补偿蒸发掉的水，还要
(i) 填满新生盐结构的孔隙空间，(ii) 送来变成晶体的那部分盐。
`tests/test_consistency.py::test_darcy_velocity_matches_the_interface_mass_balance`
直接从式(6) 重算了一遍 inset 的值，两者逐位相同。

## 参数

$\varepsilon_s$=0.43，$k_r$=2.3×10⁻³ m/s，$D_s$=1.3×10⁻⁹ m²/s，$\beta$=2，$d_{sub}$ 扫 1–100 µm。

**与原图的已知差异**：额外标注了 25 µm 处的标定锚点（诊断用途）。
