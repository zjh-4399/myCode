# 2D Bellhop 射线追踪代码优化记录

## 背景

基于参考代码（`otherCode/`）的加速技巧，对 `myCode/myCode/` 中的 MATLAB 2D Bellhop 风格高斯波束追踪代码进行逐步性能优化。

**优化原则**：
- 每次一个优化，最小化 diff，不改变功能
- 不破坏已有接口和兼容性
- 不以牺牲可读性为代价过度优化

---

## 代码架构概述

| 文件 | 职责 |
|------|------|
| `run_solver_2d.m` | 主入口，分发 ray/coherent_tl 输出类型 |
| `trace_beam_auxiliary.m` | 单条波束追踪主循环（积分 + 事件检测 + 在线累积） |
| `+utils/advance_beam_step.m` | 积分器分发（RK4 / bellhop2） |
| `+utils/beam_rhs.m` | ODE 右端项（调用 `sound_speed_derivatives`） |
| `+utils/step_beam_bellhop2.m` | 二阶中点法积分步 |
| `+utils/rk4_step_beam.m` | 经典 RK4 积分步 |
| `sound_speed_derivatives.m` | 声速剖面及其导数（热路径核心） |
| `accumulate_segment_to_targets_2d.m` | 将一段波束贡献累积到接收点（最重要的热路径） |
| `local_find_bottom_hit.m` | 预测当前步与底部边界的交点 |
| `reflect_beam_surface.m` / `reflect_beam_bottom.m` | 海面 / 海底反射 |
| `attach_bottom_normal.m` | 底部命中点的法向/切向/曲率 |
| `bottom_reflection_coeff_fluid.m` | 平面波底部反射系数 |
| `+utils/apply_ssp_interface_jump_2d.m` | SSP 层界面处的声速跳变修正 |

**波束状态向量**：
```
y = [r; z; pr; pz; Re(P); Im(P); Re(Q); Im(Q); tau; Re(G); Im(G)]
```

**支持的波束类型**（`meta.familyCode`）：
- `1` — paraxial（复数 P/Q，高斯包络）
- `2` — geometric_hat（实数 q = Re(Q)，帽形包络）
- `3` — geometric_gaussian（实数 q，高斯包络 + stent 宽度下限）

---

## 优化详情

### Opt #1–4：多文件 persistent 缓存（dispatch/config 字段访问）

**涉及文件**：`+utils/advance_beam_step.m`、`trace_beam_auxiliary.m`（`local_predict_event` 内）

**问题**：热路径函数在每次调用时重复执行 `isfield`、`lower`、`strcmpi`、`char(string(...))` 等字符串操作来读取 config 中的固定配置项。

**优化**：用 `persistent` 变量在首次调用时缓存这些标志和字符串，后续调用直接使用缓存值。

**示例（`advance_beam_step.m`）**：
```matlab
persistent cachedMethod
if isempty(cachedMethod)
    cachedMethod = local_get_integrator_name(config);
end
switch cachedMethod
    case {'rk4', 'rk4_fixed'}   ...
    case {'bellhop2', ...}       ...
end
```

**注意**：切换 config 时需调用 `clear <函数名>` 重置缓存。

---

### Opt #5：`sound_speed_derivatives` 单点位置缓存

**文件**：`sound_speed_derivatives.m`

**问题**：每个积分步中，`local_predict_event`（预测阶段）和 `beam_rhs`（积分器首步）会在完全相同的 `(r, z)` 处重复计算声速及其导数。

**调用链分析**：
```
local_predict_event → sound_speed_derivatives(y(1), y(2))   ← 计算并缓存
advance_beam_step
  └─ beam_rhs(y)     → sound_speed_derivatives(y(1), y(2))   ← 缓存命中！
  └─ beam_rhs(yMid)  → sound_speed_derivatives(yMid...)       ← cache miss，填充缓存
local_point_from_state(y_next)
  └─ sound_speed_derivatives(y_next...)                        ← cache miss
下一步 local_predict_event → sound_speed_derivatives(y_next)  ← 缓存命中！
```

**优化**：在 dispatch 缓存块之后插入标量位置缓存：
```matlab
persistent posR posZ posC posCr posCz posCrr posCrz posCzz
if isscalar(r) && isscalar(z) && ~isempty(posR) && r == posR && z == posZ
    c = posC; cr = posCr; cz = posCz;
    crr = posCrr; crz = posCrz; czz = posCzz;
    return;
end
```
函数末尾更新缓存（仅标量查询）。

**效果**：每步减少约 2 次完整 SSP 求值（bellhop2 模式下 4 次调用中命中 2 次）。

---

### Opt #6：`accumulate_segment_to_targets_2d` 向量化索引区间预计算

**文件**：`accumulate_segment_to_targets_2d.m`

**问题**：原代码在深度循环内逐深度调用 `local_receiver_index_span`，对均匀网格重复做相同的除法和 ceil/floor。

**优化**：对均匀单调接收距离网格，一次性批量计算所有活跃深度的 `[i1, i2]` 索引区间：
```matlab
if target.gridInfo.isUniform && target.gridInfo.isIncreasing
    r0_g = target.gridInfo.r0;  dr_g = target.gridInfo.dr;
    x1v = (rMinAll(depthIdx) - r0_g) / dr_g + 1.0;
    x2v = (rMaxAll(depthIdx) - r0_g) / dr_g + 1.0;
    i1v = max(1, ceil(x1v - 1.0e-12));
    i2v = min(nR, floor(x2v + 1.0e-12));
% else（非均匀网格回退，已注释保留）
end
```

同时将 `drab < 1e-12` 的零长度步检查合并进 `maskOverlap`，消除深度循环内的 `continue` 分支。

---

### Opt #7：`local_find_bottom_hit` persistent 缓存启用标志

**文件**：`local_find_bottom_hit.m`

**问题**：函数在每步调用时重复执行 7 次 `isfield` 检查来判断底部边界是否启用。

**优化**：
```matlab
persistent bottomEnabled cachedBottom cachedNSeg
if isempty(bottomEnabled)
    bottomEnabled = isfield(config, 'boundary') && ...（7个检查合并为一次）
    if bottomEnabled
        cachedBottom = config.bottom;
        cachedNSeg   = numel(config.bottom.r) - 1;
    end
end
if ~bottomEnabled; return; end
bottom = cachedBottom;  nSeg = cachedNSeg;
```

---

### Opt #8：`+utils/apply_ssp_interface_jump_2d` persistent 缓存

**文件**：`+utils/apply_ssp_interface_jump_2d.m`

**问题**：SSP 层界面跳变修正函数每次被调用时重复检查 4 个 guard 条件，并拷贝 `zw`/`cw` 数组。

**优化**：
```matlab
persistent jumpApplicable cachedZw cachedCw cachedNZw
if isempty(jumpApplicable)
    jumpApplicable = local_jump_enabled(config) && ...
        strcmpi(config.env.profile.type, 'tabular') && ...
        strcmpi(local_get_tabular_mode(config), 'piecewise_linear');
    if jumpApplicable
        cachedZw = config.env.profile.zw(:);
        cachedCw = config.env.profile.cw(:);
        cachedNZw = numel(cachedZw);
    end
end
if ~jumpApplicable; return; end
```

---

### Opt #9：`accumulate_segment_to_targets_2d` lerp + delta 预计算

**文件**：`accumulate_segment_to_targets_2d.m`

**问题**：深度循环内每次用 `(1-w)*A + w*B` 形式插值，每行需要 2 次 nRange 长向量乘法；且 `pointA`/`pointB` 字段访问在循环内重复。

**优化**：在深度循环外预计算所有差值（delta），循环内改用 lerp 形式 `A + w*(B-A)`：
```matlab
gainA = pointA.gain;  gainDelta = pointB.gain - pointA.gain;
PA    = pointA.P;     PDelta    = pointB.P    - pointA.P;
QA    = pointA.Q;     QDelta    = pointB.Q    - pointA.Q;
tauA  = pointA.tau;   tauDelta  = pointB.tau  - pointA.tau;
cA    = pointA.c;     cDelta    = pointB.c    - pointA.c;
rBeamA = pointA.r;    rBeamDelta = pointB.r   - pointA.r;
omega  = config.omega;
```
循环内：`G = gainA + w .* gainDelta;`（一次乘法，而非两次）

---

### Opt #10：`accumulate_segment_to_targets_2d` switch 提到深度循环外

**文件**：`accumulate_segment_to_targets_2d.m`

**问题**：`switch meta.familyCode` 原在深度循环内执行，每次都判断波束类型。

**优化**：将 `switch` 提到深度循环外，为 case 1/2/3 分别创建专用循环。case 2/3 不需要 `P`/`PDelta`，跳过其计算：
```matlab
switch meta.familyCode
    case 1
        for ii = 1:numel(depthIdx)  % 含 P/PDelta
            ...
        end
    case 2
        for ii = 1:numel(depthIdx)  % 不含 P/PDelta
            ...
        end
    case 3
        for ii = 1:numel(depthIdx)  % 不含 P/PDelta
            ...
        end
end
```

---

### Opt #11：`accumulate_segment_to_targets_2d` 均匀网格等比数列相位

**文件**：`accumulate_segment_to_targets_2d.m`（case 2/3 的深度循环内）

**原理**：对均匀接收距离网格，`tau = tauA + w*tauDelta` 关于接收点下标是等差数列，因此
`exp(-i·ω·tau)` 构成等比数列，可用 `cumprod` 代替逐点 `exp`：
- 原始：`nRange` 次 `exp` 调用（复数指数，开销较大）
- 优化后：2 次 `exp` + `nRange` 次复数乘法（`cumprod`）

**实现**：
```matlab
if isUniform_g
    nR_loc = numel(idx);
    phaseStart = exp(-1i * omega * tau(1));
    if nR_loc > 1
        phaseStep = exp(-1i * omega * (dr_g / drab * tauDelta));
        phase_g = cumprod([phaseStart, repmat(phaseStep, 1, nR_loc - 1)]);
    else
        phase_g = phaseStart;
    end
else
    phase_g = exp(-1i * omega * tau);
end
```

**注意**：gaussian（case 3）中 `tau` 还用于 stent 宽度下限计算，故保留 `tau` 参数，另外新增 `phase` 参数。

---

### Opt #12：reflect/coeff 文件 persistent 缓存

**涉及文件**：
- `reflect_beam_surface.m`
- `reflect_beam_bottom.m`
- `attach_bottom_normal.m`
- `bottom_reflection_coeff_fluid.m`

**reflect_beam_surface**：缓存 `config.beam.model.curvature_factor`。

**reflect_beam_bottom**：缓存 curvature_factor + plane_wave 模式标志：
```matlab
persistent cachedCurvFactor cachedPlaneWave
if isempty(cachedCurvFactor)
    cachedCurvFactor = config.beam.model.curvature_factor;
    cachedPlaneWave  = isfield(config.boundary.bottom, 'reflection') && ...
        strcmpi(config.boundary.bottom.reflection, 'plane_wave');
end
```

**attach_bottom_normal**：缓存 `config.bottom`、插值模式字符串、`use_curvature_on_facet` 标志（消除每次底部命中的 struct 拷贝、`lower()`、`isfield()` 调用）。

**bottom_reflection_coeff_fluid**：一次性预计算复数底部声速 `cbx`（从 dB/wavelength 衰减系数转换），缓存 `rho2` 和 `cbx`：
```matlab
persistent cachedRho2 cachedCbx
if isempty(cachedRho2)
    cachedRho2 = config.boundary.bottom.material.rho;
    cachedCbx  = complex_bottom_speed_from_dbwl(...);
end
```

---

### Opt #13：`local_find_bottom_hit` 向量化分段交叉检测

**文件**：`local_find_bottom_hit.m`

**问题**：候选分段上的 `for segId = segStart:segEnd` 循环逐段做 bbox 筛选和 Cramer 法求交，每段为标量操作。

**优化**：将整个循环替换为批量向量化版本：
1. 一次性提取所有候选分段端点 `r1v, z1v, r2v, z2v`
2. 批量 bbox 筛选得到 `bboxOk` 掩码
3. 对通过筛选的分段同时求解 Cramer 公式，得到 `tv`/`uv` 向量
4. `hitMask` 筛选有效命中；非命中项设为 `Inf`；`min(tv)` 取最早交点

```matlab
segRange = segStart:segEnd;
r1v = bottom.r(segRange);  z1v = bottom.z(segRange);
r2v = bottom.r(segRange+1); z2v = bottom.z(segRange+1);
bboxOk = ~((ray_rmax < min(r1v,r2v)-tolBox) | ...);
if any(bboxOk)
    ...（向量化 Cramer 求解）...
    hitMask = nonPar & (tv > tolParam) & ...;
    if any(hitMask)
        tv(~hitMask) = Inf;
        [tRaw, bestLocal] = min(tv);
        ...
    end
end
```

---

### Opt #14：`local_find_bottom_hit` persistent 空命中模板

**文件**：`local_find_bottom_hit.m`

**问题**：函数在每次调用时用 `struct(...)` 从头构造包含 8 个字段的命中结构体，该结构体在无命中时（最常见情况）直接被丢弃。

**优化**：将空命中结构体的构造合并到 `persistent` 初始化块，后续调用仅做一次 struct 拷贝（MATLAB copy-on-write，开销极低）：
```matlab
persistent bottomEnabled cachedBottom cachedNSeg emptyHit
if isempty(bottomEnabled)
    ...（原有缓存初始化）...
    emptyHit = struct('exists', false, 'rHit', nan, ...);
end
hit = emptyHit;
if ~bottomEnabled; return; end
```

---

### Opt #15：`trace_beam_auxiliary` persistent 空命中模板

**文件**：`trace_beam_auxiliary.m`（`local_empty_bottom_hit` / `local_empty_ssp_hit`）

**问题**：`local_predict_event` 在每个积分步中各调用一次 `local_empty_bottom_hit()` 和 `local_empty_ssp_hit()`，创建固定内容的结构体，共 `nBeams × nSteps` 次。

**优化**：两个辅助函数均改用 `persistent` 模板：
```matlab
function hit = local_empty_bottom_hit()
persistent tmpl
if isempty(tmpl)
    tmpl = struct('exists', false, 'rHit', nan, ...);
end
hit = tmpl;
end
```

---

## 优化效果汇总

| # | 文件 | 优化手段 | 主要收益 |
|---|------|---------|---------|
| 1–4 | 多文件 | persistent 缓存 config 字段访问 | 消除热路径中重复的 isfield/strcmpi/lower |
| 5 | `sound_speed_derivatives` | 单点位置缓存 | bellhop2 模式每步减少 ~2 次 SSP 求值 |
| 6 | `accumulate_segment_to_targets_2d` | 批量向量化索引区间 | 消除逐深度的 ceil/floor 循环调用 |
| 7 | `local_find_bottom_hit` | persistent 启用标志缓存 | 消除每步 7 次 isfield 检查 |
| 8 | `apply_ssp_interface_jump_2d` | persistent 缓存 | 消除重复 guard 检查和数组拷贝 |
| 9 | `accumulate_segment_to_targets_2d` | lerp + delta 预计算 | 每行减少 1 次 nRange 向量乘法；12 次字段访问移出循环 |
| 10 | `accumulate_segment_to_targets_2d` | switch 提到深度循环外 | 消除逐深度 switch 判断；case 2/3 跳过 P 计算 |
| 11 | `accumulate_segment_to_targets_2d` | 等比数列相位（cumprod） | nRange 次 exp → 2 次 exp + nRange 次复数乘法 |
| 12 | reflect/coeff 文件（4 个） | persistent 缓存各类参数 | 消除每次反射的多级字段访问和重复计算 |
| 13 | `local_find_bottom_hit` | 向量化分段交叉检测 | 标量循环 → 批量 Cramer 求解 + mask 筛选 |
| 14 | `local_find_bottom_hit` | persistent 空命中模板 | struct 构造 → 极低开销的 persistent 拷贝 |
| 15 | `trace_beam_auxiliary` | persistent 空命中模板 | 每步 2 次 struct 构造 → persistent 拷贝 |

---

## 重要注意事项

### persistent 缓存重置
切换 `config`（不同声速剖面、不同边界条件）时，需手动重置各文件的 persistent 变量：
```matlab
clear sound_speed_derivatives
clear local_find_bottom_hit
clear attach_bottom_normal
clear reflect_beam_surface
clear reflect_beam_bottom
clear bottom_reflection_coeff_fluid
clear utils.advance_beam_step
clear utils.apply_ssp_interface_jump_2d
```

### 非均匀接收网格
Opt #6 中非均匀网格的回退代码已**注释保留**（未删除），以便将来需要时恢复：
```matlab
% else
%     i1v = zeros(nActive, 1);
%     ...
end
```

### gaussian stent 宽度下限
Opt #11 中 gaussian（case 3）的 `tau` 参数不能用 `phase` 替换，因为
`Wfloor = min(0.2*freq*max(real(tau),0), pi*lambda)` 依赖实际旅行时。
函数签名改为同时传入 `tau`（用于 stent）和 `phase`（用于指数项）。
