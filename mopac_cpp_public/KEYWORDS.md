# MOPAC_C++ 输入关键词使用说明

本程序是 MOPAC2016（Fortran）的 C++17 整库翻译版，**输入文件格式与关键词体系与 MOPAC2016 官方版本一致**（默认路径与官方参考程序逐位一致）。本文按类别列出程序实际解析的全部输入关键词及其功能。

---

## 1. 输入文件格式（.mop）

```
<关键词行>            # 第一行：关键词，空格分隔，不区分大小写，可用 MOPAC 官方缩写前缀
<注释行 1>            # 第二行：标题/注释
<注释行 2>            # 第三行：注释（可选）

<原子坐标行>          # 每行一个原子：元素  x  y  z  [净电荷]
...
                     # 末尾空行结束
```

坐标支持 **Z 矩阵**（默认）与 **XYZ 笛卡尔坐标**（加 `XYZ` 关键词）；原子坐标行格式：

- Z 矩阵：`C 0 1.54 1 109.5 2 120.0 0`（原子、连接、键长、角度、二面角）
- 笛卡尔：`C 0.0000 0.0000 0.0000`（需 `XYZ` 关键词）

**运行**：`mopac.exe <输入名>.mop` → 输出写入 `<输入名>.out`。

---

## 2. 计算方法（哈密顿量）

| 关键词 | 功能 |
|---|---|
| `AM1` | Austin Model 1 半经验哈密顿量 |
| `PM3` | Parametric Method 3 |
| `PM6` | Parametric Method 6 |
| `PM6-D` | PM6 + D3 色散修正 |
| `PM6-H` | PM6 + 氢键修正 |
| `PM6-DH+` | PM6 + D3 色散 + 氢键修正 |
| `PM7` | Parametric Method 7（默认推荐，精度最高） |
| `PM7-TS` | PM7 过渡态专用参数 |
| `D3H4` | Grimme D3(BJ) 色散 + H4 氢键修正项 |
| `MOZ` | MOZYME 局域化近似（大分子，数千原子） |
| `STO3G` | STO-3G 极小基组（用于 ESP 静电势拟合等） |

---

## 3. 任务类型

| 关键词 | 功能 |
|---|---|
| `1SCF` | 单点计算（不优化几何） |
| `EF` | 特征向量跟随法几何优化（默认优化器） |
| `BFGS` | BFGS 拟牛顿优化器 |
| `DFP` | DFP 变尺度优化器 |
| `LBFGS` | 有限内存 BFGS 优化器 |
| `FLEPO` | Fletcher-Powell 优化器 |
| `NLLSQ` | 非线性最小二乘优化（适合过渡态搜索） |
| `POWSQ` | 幂次法优化器 |
| `TS` | 过渡态搜索（一阶鞍点） |
| `SADDLE` | 鞍点搜索 |
| `FORCETS` | 力常数法过渡态搜索 |
| `SIGMA` | Sigma 方法过渡态搜索 |
| `LOCATE-TS` | 双端点 QST 过渡态定位 |
| `REFINE-TS` | 精化已有的过渡态猜测 |
| `DRC` | 动态反应坐标（动力学轨迹） |
| `IRC` | 内禀反应坐标 |
| `MOLDYN` | 分子动力学模拟 |
| `THERMO` | 热力学性质计算（熵、热容等） |
| `FORCE` | 力常数矩阵（振动频率）计算 |
| `DFORCE` | 解析力常数 |
| `POLAR` | 极化率计算 |
| `RAPID` | 快速模式（配合 MOZYME 使用） |
| `CAMP-KING` | Camp-King 过渡态搜索方法 |

---

## 4. SCF 收敛控制

| 关键词 | 功能 | 参数 |
|---|---|---|
| `PULAY` | Pulay DIIS 外推加速 SCF 收敛 | — |
| `SHIFT=` | Hückel 能级移动（稳定 SCF，-80 默认） | 数值 |
| `FILL=` | 设置占据轨道数（如开壳层） | 数值 |
| `ITRY=` | SCF 最大迭代次数 | 数值 |
| `SCFCRT=` | SCF 收敛判据（能量，默认 1.D-10） | 数值 |
| `RELSCF` | 相对论修正 | — |
| `OLDENS` | 用旧密度矩阵热启动 SCF | — |
| `DEBUGPULAY` | PULAY 调试输出 | — |
| `DEBUG` | 调试输出（配合 ITER/FOCK 等） | — |

---

## 5. 几何优化控制

| 关键词 | 功能 | 参数 |
|---|---|---|
| `GNORM=` | 梯度收敛标准（默认 0.01 kcal/mol·Å） | 数值 |
| `CYCLES=` | 最大优化循环数 | 数值 |
| `BIGCYCLES=` | DRC/MD 大循环数 | 数值 |
| `PREC` | 提高精度（更严收敛判据） | — |
| `PRECISE` | 高精度模式（SCFCRT 缩紧 100 倍） | — |
| `LET` | 允许不太合理的几何（跳过部分检查） | — |
| `GEO-OK` | 跳过几何合理性检查 | — |
| `NOUPD` | EF 优化中不更新 Hessian | — |
| `NONR` | EF 中不重新计算 Hessian | — |
| `RSCAL` | EF Hessian 缩放 | — |
| `GRAD` / `GRADI` | 输出梯度 | — |
| `STEP1` / `STEP2` | 输出优化步进细节 | — |
| `INTERP` | 几何内插 | — |
| `LINMIN` | 线搜索最小化 | — |
| `NOOPT` | 禁止几何优化（与 OPT 相反） | — |
| `OPT` | 强制几何优化 | — |

---

## 6. 性质计算与输出

| 关键词 | 功能 |
|---|---|
| `POP` | 布居分析（Mulliken 电荷输出） |
| `MULLIK` / `MULL` | Mulliken 布居 |
| `CHARGES` | 输出原子净电荷 |
| `PRTCHAR` | 输出电荷 |
| `ESP` | 静电势拟合点电荷 |
| `DIPOLE` | 偶极矩 |
| `PMEP` | 分子静电势（MEP）计算 |
| `PRTMEP` | 输出 MEP 网格 |
| `MINMEP` | 最小静电势搜索 |
| `ENPART` | 能量分解（双原子贡献） |
| `BONDS` | 键级/键序输出 |
| `PI` | π 键级 |
| `SPIN` | 自旋密度 |
| `EIGS` | 输出轨道本征值 |
| `VEC` / `VECTORS` | 输出本征向量 |
| `ALLVEC` / `ALLV` | 输出全部本征向量 |
| `DENS` | 输出密度矩阵 |
| `FOCK` | 输出 Fock 矩阵 |
| `1ELE` | 输出单电子矩阵 |
| `HCORE` | 输出 Hcore 矩阵 |
| `HESSIAN` | 输出 Hessian |
| `COMPFG` | compfg 详细输出 |
| `THERM` | 热力学量输出 |
| `GRAPH` / `GRAPHF` | 分子图/键序图输出 |
| `TIDY` | 整理并输出规范化几何 |
| `SILENT` | 安静模式（减少输出） |

---

## 7. 组态相互作用（CI）

| 关键词 | 功能 |
|---|---|
| `C.I.` / `MECI` | 多电子组态相互作用 |
| `PECI` | 完全组态相互作用（价电子全激发） |
| `MICROS` | 微态（指定占据数）CI |
| `CIOSCI` | CI-OSCI 方法 |
| `EXCI` | 指定激发数 |
| `ROOT=` | 计算的第 n 个根 | 
| `OPEN` | 开壳层 CI |
| `BIRAD` | 双自由基 |
| `SING` / `DOUB` / `TRIP` | 单/双/三重态 |
| `QUAR` / `QUIN` / `SEXT` / `SEPT` / `OCTE` | 四/五/六/七/八重态 |
| `ESR` | 电子自旋共振参数（超精细耦合） |
| `NONE` | 关闭 CI |
| `NOANCI` | 关闭解析 CI 导数 |

---

## 8. 输入输出与坐标

| 关键词 | 功能 |
|---|---|
| `XYZ` | 笛卡尔坐标输入 |
| `ANG` | 键长单位用埃（默认） |
| `COORD` | 输出笛卡尔坐标 |
| `PDBOUT` | 输出 PDB 文件 |
| `RESEQ` | 原子重新编号 |
| `RESID` | 保留残基标识 |
| `NORES` | 不输出残基标识 |
| `NORESEQ` | 不重排残基 |
| `ADD-H` | 自动加氢 |
| `SITE` | 指定活性位点 |
| `DENOUTF` | 密度矩阵格式化输出文件 |
| `AIGIN` / `AIGOUT` | Archive 文件读入/写出 |
| `MOLDAT` | 输出 moldat 文件 |
| `AUX` | 输出辅助文件（Aux） |
| `HTML` | 输出 HTML 格式 |
| `LOG` | 输出日志文件 |
| `OUTPUT` | 指定输出文件名 |
| `CHECK` | 输入检查模式 |
| `THREADS=` | 设置并行线程数 | 
| `CHAINS` | PDB 多链处理 |
| `OLDGEO` | 从旧几何启动 |
| `RESTART` | 从断点重启 |
| `NSPA` | 特殊原子处理（Sparkle 体系） |
| `PINOUT` | 输出 PIN 文件 |
| `SETPI` | 设置 π 键体系 |
| `LEWIS` | Lewis 结构检查 |
| `PARASOK` | 参数检查（Sokalski） |

---

## 9. 对称性

| 关键词 | 功能 |
|---|---|
| `SYMM` | 启用分子对称性（DRC 中） |
| `NOSYM` | 禁用对称性 |
| `NOREOR` | 不重定向分子 |
| `MOLSYM` | 分子对称性分析 |
| `SYMOIR` | 对称性约化 |
| `SYMTRZ` | 对称性变换 |
| `CARTAB` | 笛卡尔表输出 |
| `ALLBOND` | 输出全部键 |

---

## 10. 电场、溶剂化与特殊模型

| 关键词 | 功能 |
|---|---|
| `FIELD` | 外加电场 |
| `EXTERNAL` | 外部势/电荷 |
| `QMMM` | QM/MM 混合计算 |
| `SPARKLE` | Sparkle 模型（镧系/锕系） |
| `WILLIAMS` | Williams 静电势拟合（ESP 权重） |
| `SYMAVG` | 对称平均（ESP） |
| `ESPGRID` | ESP 网格输出 |
| `ESPRST` | ESP 重启 |
| `POTWRT` | 静电势写出 |
| `CUBE` | 输出 cube 文件（密度/势） |
| `DISP` | 色散校正输出 |
| `NODIRECT` | 禁用直接积分（整库存储） |
| `SEMIDIRECT` | 半直接积分 |
| `DAMP` | 阻尼因子（SCF 收敛辅助） |
| `DIAGG2` | 对角化调试输出 |

---

## 11. DRC / 分子动力学专用

| 关键词 | 功能 | 参数 |
|---|---|---|
| `VELO=` | 指定初始速度（速度标度） | 数值 |
| `KINE=` | 指定初始动能 | 数值 |
| `X-PRIO` | 优先轴设置 | — |
| `EPS=` | 能量判据 | 数值 |
| `SMOOTH` | 势能面平滑 | — |
| `POINT` / `POINT1` / `POINT2` | 势能面扫描点控制 | 数值 |

---

## 12. 自定义扩展（本程序新增，MOPAC2016 官方无此关键词）

| 关键词 | 功能 |
|---|---|
| `FULLDIAG` | **可选加速对角化模式**：SCF 改用完整 LAPACK `dsyevd` 对角化（默认走官方伪对角化路径，与 Fortran 参考逐位一致）。适合大体系（1SCF 约 2–3.6× 提速）与近简并体系（SCF 收敛更稳健，如 C40H82 案例）。与设置环境变量 `MOPAC_FULL_DIAG=1` 等价。优化终点与默认路径数学等效（NIST 验证）。 |

---

## 13. 环境变量（运行时开关，全部默认关闭）

| 变量 | 功能 |
|---|---|
| `MOPAC_FULL_DIAG` | 等效于 `FULLDIAG` 关键词（见上） |
| `MOPAC_DSPAR` / `MOPAC_DSPAR_NB` | 确定性并行 Fock 归约（FBDR，NB=桶数） |
| `MOPAC_GAB` | EF 优化器中简并基的确定性选择 |
| `EFPROBE` / `HESSPROBE` | 非侵入式诊断转储 |
| `MOPAC_ULPPERT` | ±1 ULP 扰动注入（研究用） |

所有变量不设置时，程序与官方 Fortran 路径逐位一致。

---

## 14. 示例

**单点计算**（水分子，PM7，笛卡尔坐标）：
```
PM7 1SCF XYZ
water single point

O 0.0000 0.0000 0.1173
H 0.0000 0.7572 -0.4692
H 0.0000 -0.7572 -0.4692

```

**几何优化**（PM7，梯度收敛 0.01，Pulay 加速）：
```
PM7 GNORM=0.01 PULAY EF
ethanol optimisation

C 0.0000 0.0000 0.0000
...

```

**过渡态搜索 + FULLDIAG 加速**：
```
PM7 TS FULLDIAG
transition state search

...

```

---

*说明：本清单由程序源码关键词解析逻辑提取整理；未列出的 MOPAC2016 官方关键词在本程序中同样按官方语义处理（关键词前缀缩写规则与官方一致）。详细数学定义见 MOPAC2016 官方手册。*
