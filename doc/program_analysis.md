# COLA Halo OMP 项目分析

## 主程序确认

**主程序文件**: `main.c` (Makefile第30行 `OBJS := main.o` 是第一个目标文件)

## 程序执行流程

### 1. 初始化阶段 (main函数开始)

| 行号 | 函数调用 | 所在文件 | 功能 |
|------|----------|----------|------|
| 135 | `mpi_init()` | main.c:608 | MPI初始化 |
| 136 | `msg_init()` | msg.c | 消息系统初始化 |
| 143 | `read_parameters()` | read_param_lua.c | 读取参数 (若定义USE_LUA) |
| 161 | `confirm_parameters()` | confirm_param.c | 确认/打印参数 |
| 162 | `gb_solve_D_E_Qs()` | solve_growth.c | 求解增长因子 |
| 164 | `fft_init()` | main.c:638 | FFTW初始化 |
| 165 | `comm_init()` | comm.c | MPI通信初始化 |
| 200 | `power_init()` | power.c | 功率谱初始化 |
| 205 | `allocate_shared_memory()` | mem.c | 内存分配 |
| 206 | `lpt_init()` | lpt.c | LPT初始化 |
| 209 | `allocate_particles()` | mem.c | 分配粒子结构 |
| 220 | `pm_init()` | pm.c | PM方法初始化 |
| 222 | `fof_init()` | fof.c | FOF晕寻找初始化 |
| 223 | `subsample_init()` | subsample.c | 子采样初始化 |

### 2. 多重实现循环 (irealization循环, main.c:237)

| 行号 | 函数调用 | 所在文件 | 功能 |
|------|----------|----------|------|
| 248 | `lpt_set_displacement()` | lpt.c | 设置初始网格和2LPT位移 |
| 257 | `write_snapshot()` | write.c | 可选：写入初始条件 |

### 3. 时间演化循环 (istep循环, main.c:321)

| 行号 | 函数调用 | 所在文件 | 功能 |
|------|----------|----------|------|
| 326 | `move_particles2()` | move.c | 粒子移动/通信 |
| 330 | `pm_calculate_forces()` | pm.c | PM计算力 |
| 355 | `cola_kick()` | cola.c | Leap-frog kick (速度更新) |
| 365 | `cola_drift()` | cola.c | Leap-frog drift (位置更新) |

### 4. 输出阶段 (在时间演化循环中)

| 行号 | 函数调用 | 所在文件 | 功能 |
|------|----------|----------|------|
| 348/359 | `snapshot_time()` | main.c:713 | 快照输出入口函数 |
| 735 | `cola_set_snapshot()` | cola.c | 设置快照数据 |
| 743 | `fof_find_halos()` | fof.c | FOF晕寻找 |
| 745 | `fof_write_halos()` | fof.c | 写入晕数据 |
| 753 | `write_snapshot()` | write.c | 写入快照文件 |
| 764 | `write_random_sabsample()` | write.c | 写入子采样 |
| 770 | `coarse_grid2()` | coarse_grid.c | 写入粗网格 |

### 5. 光锥输出 (xiaodong添加的功能)

- 在时间循环中处理光锥输出 (main.c:387-564)
- 包括边界穿越检查逻辑

### 6. 结束阶段

| 行号 | 函数调用 | 所在文件 | 功能 |
|------|----------|----------|------|
| 589 | `timer_print()` | timer.c | 打印计时信息 |
| 604 | `MPI_Finalize()` | - | MPI结束 |

## 文件依赖关系汇总

按首次调用顺序排列：

| 序号 | 文件 | 主要功能 | 对应头文件 |
|------|------|----------|------------|
| 1 | `read_param_lua.c` | 参数读取(LUA格式) | parameters.h |
| 2 | `confirm_param.c` | 参数确认 | parameters.h |
| 3 | `solve_growth.c` | 增长因子求解 | solve_growth.h |
| 4 | `msg.c` | 消息输出系统 | msg.h |
| 5 | `comm.c` | MPI通信 | comm.h |
| 6 | `power.c` | 功率谱处理 | power.h |
| 7 | `mem.c` | 内存管理 | mem.h |
| 8 | `lpt.c` | 拉格朗日扰动理论位移 | lpt.h |
| 9 | `pm.c` | 粒子网格法力计算 | pm.h |
| 10 | `cola.c` | COLA时间积分 | cola.h |
| 11 | `move.c` | 粒子移动/边界处理 | move.h |
| 12 | `move_min.c` | 粒子移动(最小数据版) | move_min.h (自动生成) |
| 13 | `fof.c` | Friends-of-Friends晕寻找 | fof.h |
| 14 | `write.c` | 数据输出 | write.h |
| 15 | `subsample.c` | 子采样功能 | subsample.h |
| 16 | `coarse_grid.c` | 粗网格生成 | coarse_grid.h |
| 17 | `timer.c` | 计时统计 | timer.h |

## 关键数据结构

```c
// 粒子结构 (particle.h)
typedef struct {
    float x[3];      // 位置
    float dx1[3];    // ZA位移
    float dx2[3];    // 2LPT位移
    float v[3];      // 速度
    long long id;    // ID
} Particle;

// 粒子集合 (particle.h)
typedef struct {
    Particle* p;
    float3* force;
    float a_x, a_v;
    int np_local, np_allocated;
    long long np_total;
    float np_average;
} Particles;

// 快照结构 (particle.h)
typedef struct {
    ParticleMinimum* p;
    int nc, np_local;
    long long np_total;
    float boxsize, omega_m, h;
    int seed;
    char* filename;
} Snapshot;
```

## 编译信息

- 编译器: `mpicc -std=c99`
- 并行: MPI + OpenMP混合并行
- 依赖库: FFTW3 (MPI + OMP), GSL, Lua5.4
- 可执行文件: `cola_halo`

## 程序类型

**宇宙学N体模拟程序** - 使用COLA(COmoving Lagrangian Acceleration)方法模拟暗物质晕的形成和演化。

## 待分析项目

- [x] PM方法内存分配分析
- [ ] 各模块内部函数实现细节
- [ ] 2LPT位移计算方法
- [ ] PM力计算算法
- [ ] COLA积分方案
- [ ] FOF晕寻找算法
- [ ] 光锥输出逻辑细节

---

## 内存分析

- 详细内存分析见: [memory_analysis.md](memory_analysis.md)

---

## LPT模块深度分析 (lpt.c)

### 模块概述

`lpt.c` 实现了**拉格朗日扰动理论(Lagrangian Perturbation Theory)**的位移场计算,包括一阶(Zel'dovich Approximation, ZA)和二阶(2LPT)位移场。这是COLA模拟的初始条件生成模块。

**代码来源**: 基于Roman Scoccimaro等人的2LPT代码 (http://cosmo.nyu.edu/roman/2LPT/)

### 核心数据结构

```c
// 位移场存储 (全局静态变量)
fftwf_complex *(cdisp[3]), *(cdisp2[3]);  // ZA和2LPT位移(傅里叶空间)
float         *(disp[3]), *(disp2[3]);     // ZA和2LPT位移(实空间)
fftwf_complex *(cdigrad[6]);                // ZA位移梯度的6个独立分量
float         *(digrad[6]);

// MPI并行参数
static int Local_nx, Local_x_start;         // 本进程负责的x方向网格范围
static int Nmesh;                            // 每维网格数
```

### 主要函数功能

#### 1. `lpt_init()` - 初始化函数

**功能**: 分配内存并创建FFT计划

**关键步骤**:
- 计算MPI并行下的本地数据大小 (`fftwf_mpi_local_size_3d`)
- 分配12个复数数组:
  - `cdisp[3]`, `cdisp2[3]`: ZA和2LPT位移场
  - `cdigrad[6]`: ZA位移梯度的6个独立分量 (∂ψ_i/∂q_j, i≤j)
- 创建10个FFT计划:
  - 6个逆变换计划 (复数→实数)
  - 1个正变换计划 (实数→复数)
  - 3个位移场逆变换计划
- 生成随机种子表 (确保可重复性)

**内存需求**: 约 `12 × Nmesh³ × sizeof(fftwf_complex)` 字节

#### 2. `lpt_set_displacement()` - 核心计算函数

**功能**: 生成随机高斯位移场并设置粒子初始位置

**参数**:
- `InitTime`: 初始时刻(尺度因子a)
- `omega_m`: 物质密度参数
- `Seed`: 随机数种子
- `Box`: 模拟盒子大小
- `particles`: 粒子结构体指针

**算法流程**:

##### 步骤1: 设置宇宙学参数和速度前置因子 (行145-163)

```c
// Hubble常数和速度前置因子计算
const double hubble_a = Hubble * sqrt(Omega/a³ + (1-Omega-ΩΛ)/a² + ΩΛ);
double vel_prefac  = a * hubble_a * F_Omega(a);   // 一阶速度因子
double vel_prefac2 = a * hubble_a * F2_Omega(a);  // 二阶速度因子
```

其中 `F_Omega(a)` 和 `F2_Omega(a)` 是增长因子的宇宙学修正函数。

##### 步骤2: 生成随机种子表 (行169-204)

使用复杂的对称性填充方式生成 `Nmesh × Nmesh` 的种子表,确保:
- 不同k模式使用不同随机种子
- 满足实数场的傅里叶共轭对称性
- 可重复性(相同Seed产生相同结果)

##### 步骤3: 生成ZA位移场 (行206-348)

**核心算法**: 在傅里叶空间生成高斯随机场

```c
for each k-mode (i,j,k):
    // 1. 生成随机相位和振幅
    phase = 2π × random()
    ampl = -log(random())
    
    // 2. 计算功率谱
    p_of_k = PowerSpec(|k|)
    
    // 3. 构造位移场(傅里叶空间)
    // ψ(k) = -i k/k² × δ(k)
    // 其中 δ(k) = √P(k) × √(-log(ampl)) × e^{i×phase}
    cdisp[axes][k] = -k[axes]/|k|² × δ(k) × sin/cos(phase)
```

**关键技术点**:
- 避免Nyquist频率 (i,j,k = Nmesh/2)
- 避免k=0模式 (零模式)
- 处理k=0平面的共轭对称性 (行286-343)
- MPI进程间的共轭模式可能分布在不同进程

##### 步骤4: 计算2LPT位移场 (行352-491)

**理论基础**: 二阶位移场由一阶位移的梯度张量构造

```
ψ⁽²⁾_i = ∇⁻² × S_i
其中 S_i = (∂ψ⁽¹⁾_j/∂q_i)(∂ψ⁽¹⁾_i/∂q_j) - (∂ψ⁽¹⁾_i/∂q_j)(∂ψ⁽¹⁾_j/∂q_k)
```

**计算步骤**:

1. **计算ZA位移梯度** (行355-395):
   ```c
   // 在傅里叶空间: ∂ψ_i/∂q_j → i×k_j×ψ_i(k)
   cdigrad[0] = i×k[0]×cdisp[0]  // ∂ψ_0/∂q_0
   cdigrad[1] = i×k[1]×cdisp[0]  // ∂ψ_0/∂q_1
   cdigrad[2] = i×k[2]×cdisp[0]  // ∂ψ_0/∂q_2
   cdigrad[3] = i×k[1]×cdisp[1]  // ∂ψ_1/∂q_1
   cdigrad[4] = i×k[2]×cdisp[1]  // ∂ψ_1/∂q_2
   cdigrad[5] = i×k[2]×cdisp[2]  // ∂ψ_2/∂q_2
   ```

2. **逆变换到实空间** (行397-402):
   ```c
   fftwf_mpi_execute_dft_c2r(Inverse_plan[i], cdigrad[i], digrad[i]);
   ```

3. **构造二阶源项** (行404-419):
   ```c
   // 实空间计算
   source = ∂ψ_0/∂q_0 × (∂ψ_1/∂q_1 + ∂ψ_2/∂q_2) 
          + ∂ψ_1/∂q_1 × ∂ψ_2/∂q_2
          - (∂ψ_0/∂q_1)² - (∂ψ_0/∂q_2)² - (∂ψ_1/∂q_2)²
   ```

4. **正变换回傅里叶空间** (行421-424):
   ```c
   fftwf_mpi_execute_dft_r2c(Forward_plan, digrad[3], cdigrad[3]);
   ```

5. **求解泊松方程得到2LPT位移** (行435-491):
   ```c
   // ψ⁽²⁾(k) = i×k/|k|² × source(k)
   cdisp2[axes][k] = i×k[axes]/|k|² × source(k)
   ```

##### 步骤5: 变换到实空间 (行495-500)

```c
for each axis:
    fftwf_mpi_execute_dft_c2r(Disp_plan[axis], cdisp[axis], disp[axis]);
    fftwf_mpi_execute_dft_c2r(Disp2_plan[axis], cdisp2[axis], disp2[axis]);
```

##### 步骤6: 设置粒子位置和速度 (行502-558)

**关键公式**:

```c
// 初始位置
x[axes] = q[axes] + D⁺(a)×ψ⁽¹⁾(q) - (3/7)×D₂(a)×ψ⁽²⁾(q)

// 存储外推到a=1的位移 (用于COLA方法)
p->dx1[axes] = ψ⁽¹⁾(q)                    // 一阶位移
p->dx2[axes] = -(3/7)×D₂₀×ψ⁽²⁾(q)       // 二阶位移

// 初始速度设为0 (COLA会在后续步骤中计算)
p->v[axes] = 0.0f
```

其中:
- `D⁺(a)` = 增长因子 (从InitTime到a=1)
- `D₂(a)` = 二阶增长因子,使用 `D₂ ∝ D⁺² × (Ω(a)/Ω_m)^{-1/143}` 近似
- `D₂₀` = a=1时的二阶增长因子归一化
- `-3/7` = 二阶LPT的理论系数

**粒子ID分配**: 按照网格顺序分配全局唯一ID

```c
id = Local_x_start × Nmesh² + 1  // 起始ID
```

#### 3. 辅助函数

##### `F_Omega(a)` - 一阶增长因子修正 (行591-597)

```c
F_Omega(a) = [Ω(a)]^0.6
其中 Ω(a) = Ω_m / (Ω_m + a×(1-Ω_m-Ω_Λ) + a³×Ω_Λ)
```

这是Peebles的增长因子近似公式。

##### `F2_Omega(a)` - 二阶增长因子修正 (行600-606)

```c
F2_Omega(a) = 2 × [Ω(a)]^{4/7}
```

用于二阶扰动场的增长因子修正。

### 关键技术特点

#### 1. MPI并行策略

- **数据分布**: x方向网格分解 (`Local_nx`, `Local_x_start`)
- **傅里叶变换**: 使用 `fftwf_mpi_plan_dft_*` 自动处理数据分布
- **共轭对称性**: 处理跨进程的复数共轭模式 (行312-342)

#### 2. 傅里叶空间处理

- **存储格式**: FFTW的实数输入/复数输出格式
  - 实空间: `Nmesh × Nmesh × Nmesh`
  - 傅里叶空间: `Nmesh × Nmesh × (Nmesh/2+1)` (利用实数场的共轭对称性)
- **内存对齐**: 使用 `2×(Nmesh/2+1)` 以支持in-place变换

#### 3. 随机数生成

- **算法**: GSL库的 `gsl_rng_ranlxd1` (Luxury random numbers)
- **种子表**: 预生成 `Nmesh²` 个种子,确保:
  - 不同k模式独立
  - 相同种子产生相同结果
  - 满足傅里叶共轭对称性

#### 4. 功率谱截断

```c
#ifdef SPHEREMODE
    if(|k| > Nsample/2) continue;  // 球形截断
#else
    if(|k_i| > Nsample/2) continue; // 立方体截断
#endif
```

#### 5. CIC修正 (可选)

```c
#ifdef CORRECT_CIC
    // 去卷积CIC插值的平滑效应
    smth = 1 / [sinc(k_x/2) × sinc(k_y/2) × sinc(k_z/2)]²
    cdisp[axes][k] *= smth;
#endif
```

### 内存使用分析

对于 `Nmesh = nc`:

| 数组 | 数量 | 大小(每个) | 总大小 |
|------|------|-----------|--------|
| `cdisp` | 3 | `nc³ × sizeof(fftwf_complex)` | 6×nc³ 字节 |
| `cdisp2` | 3 | `nc³ × sizeof(fftwf_complex)` | 6×nc³ 字节 |
| `cdigrad` | 6 | `nc³ × sizeof(fftwf_complex)` | 12×nc³ 字节 |
| **总计** | 12 | - | **24×nc³ 字节** |

**示例**: nc=512 → 约 3.2 GB

### 与COLA方法的集成

LPT模块为COLA时间积分提供关键数据:

1. **初始位置**: `p->x = q + D⁺×ψ⁽¹⁾ - (3/7)×D₂×ψ⁽²⁾`
2. **LPT位移轨迹**: 
   - `p->dx1` = ZA位移 (外推到a=1)
   - `p->dx2` = 2LPT位移 (外推到a=1)
3. **COLA公式**: `x(a) = x_LPT(a) + x_residual(a)`
   - LPT部分: `x_LPT(a) = q + D⁺(a)×dx1 + D₂(a)×dx2`
   - 残余部分: 由PM方法计算

---

## PM模块深度分析 (pm.c)

### 模块概述

`pm.c` 实现了**粒子网格(Particle Mesh, PM)方法**的引力计算,这是COLA模拟中计算残余引力的核心模块。PM方法通过将粒子质量分配到网格上,在傅里叶空间求解泊松方程,高效地计算长程引力。

**代码来源**: 基于Svetlin Tassev的COLAcode (Harvard/Princeton),采用GNU GPL v3许可

### 核心数据结构

```c
// PM网格参数
static int Ngrid;                    // 每维网格数
static int PM_factor;                // PM网格细化因子
static float BoxSize;                // 模拟盒子大小

// MPI并行参数
static int Local_nx, Local_x_start;  // x方向分解
static int Local_ny_td, Local_y_start_td;  // y方向分解(转置后)

// 傅里叶变换数据
static fftwf_complex* fftdata;       // 密度/力场(实空间+傅里叶空间)
static fftwf_complex* density_k;     // 密度场(傅里叶空间,转置存储)

// FFTW计划
static fftwf_plan p0;                // 正变换: r2c (实数→复数)
static fftwf_plan p11;               // 逆变换: c2r (复数→实数)

// 边界粒子缓冲区
static BufferVec3 BufPos;            // 位置缓冲区(用于MPI通信)
```

### 主要函数功能

#### 1. `pm_init()` - 初始化函数 (行87-198)

**功能**: 设置PM网格、分配内存、创建FFT计划

**关键步骤**:

1. **计算MPI并行分解** (行98-106):
   ```c
   // 使用转置输出优化性能
   fftwf_mpi_local_size_3d_transposed(Ngrid, Ngrid, Ngrid/2+1, 
       &local_nx, &local_x_start,    // x方向分解(实空间)
       &local_ny, &local_y_start);   // y方向分解(傅里叶空间转置)
   ```

2. **分配内存** (行110-127):
   - `fftdata`: 主数据数组,用于密度场和力场
   - `density_k`: 存储傅里叶空间的密度场(转置格式)

3. **创建FFT计划** (行141-147):
   ```c
   // 正变换: 实空间密度 → 傅里叶空间密度(转置输出)
   p0 = fftwf_mpi_plan_dft_r2c_3d(..., FFTW_MPI_TRANSPOSED_OUT);
   
   // 逆变换: 傅里叶空间力 → 实空间力(转置输入)
   p11 = fftwf_mpi_plan_dft_c2r_3d(..., FFTW_MPI_TRANSPOSED_IN);
   ```

4. **确定相邻进程** (行149-169):
   ```c
   // 找到左右邻居进程(周期性边界)
   LeftNode = 左邻居进程ID;
   RightNode = 右邻居进程ID;
   ```

5. **分配粒子缓冲区** (行176-194):
   ```c
   // 缓冲区大小: 平均粒子数 + 5σ涨落
   nbuf = np_alloc_factor×ncp² + 5×sqrt(np_alloc_factor×ncp²);
   ```

**内存需求**: 
- `fftdata`: `Ngrid³ × sizeof(fftwf_complex)`
- `density_k`: `(Ngrid/2+1) × Ngrid × Local_ny_td × sizeof(fftwf_complex)`
- 总计约 `2 × Ngrid³ × sizeof(fftwf_complex)` 字节

#### 2. `pm_calculate_forces()` - 主计算函数 (行565-596)

**功能**: 计算粒子受到的引力

**算法流程**:

```
1. send_buffer_positions()    // 交换边界粒子
2. PtoMesh()                  // 粒子→网格 (CIC质量分配)
3. check_total_density()      // 检查质量守恒
4. compute_density_k()        // 傅里叶变换
5. compute_force_mesh()       // 计算力场 (3个分量)
6. force_at_particle_locations() // 网格→粒子 (CIC插值)
7. add_buffer_forces()        // 合并边界粒子力
```

#### 3. `PtoMesh()` - 粒子到网格的质量分配 (行209-279)

**功能**: 使用CIC(Cloud-In-Cell)插值将粒子质量分配到网格

**CIC插值原理**:

每个粒子被看作一个均匀密度的立方体"云",其质量分配到最近的8个网格点:

```
权重 = W(x) × W(y) × W(z)
其中 W(t) = max(0, 1 - |t|)  (三角形函数)
```

**实现细节**:

```c
// 1. 计算网格索引和距离
int iI = floor(X);  // 左网格索引
float D1 = X - iI;  // 到左网格的距离
float T1 = 1.0 - D1; // 到右网格的距离

// 2. 周期性边界处理 (y,z方向)
if(J >= Ngrid) J = 0;
if(K >= Ngrid) K = 0;

// 3. 分配质量到8个网格点
WRtPlus(density, iI, J,  K,  T3*T1*T2W);  // 左-下-后
WRtPlus(density, iI, J,  K1, D3*T1*T2W);  // 左-下-前
WRtPlus(density, iI, J1, K,  T3*T1*D2W);  // 左-上-后
WRtPlus(density, iI, J1, K1, D3*T1*D2W);  // 左-上-前
// ... 右侧4个点类似
```

**关键点**:
- `WPAR = PM_factor³`: 质量权重(考虑PM网格细化)
- x方向**不**做周期性处理,而是通过缓冲粒子处理边界
- 使用 `#pragma omp atomic` 保证线程安全

#### 4. `compute_density_k()` - 傅里叶变换 (行282-301)

**功能**: 将实空间密度场变换到傅里叶空间

```c
// 正变换: ρ(x) → ρ(k)
fftwf_mpi_execute_dft_r2c(p0, (float*) fftdata, fftdata);

// 复制到density_k (转置格式)
for(Jl = 0; Jl < Local_ny_td; Jl++)
  for(iI = 0; iI < Ngrid; iI++)
    for(K = 0; K < Ngrid/2+1; K++)
      density_k[index] = fftdata[index];
```

**转置优化**: 使用 `FFTW_MPI_TRANSPOSED_OUT`,输出数据按y方向分解,避免额外的通信开销。

#### 5. `compute_force_mesh()` - 力场计算 (行304-355)

**功能**: 在傅里叶空间计算引力场

**理论基础**: 泊松方程在傅里叶空间的解

```
∇²Φ = 4πGρ  →  Φ(k) = -4πG ρ(k)/|k|²
F_i = -∂Φ/∂x_i  →  F_i(k) = -i k_i Φ(k) = i k_i × 4πG ρ(k)/|k|²
```

**实现**:

```c
// 归一化因子
f1 = -1.0 / Ngrid³ / (2π/BoxSize);

for each k-mode (I0, J0, K):
    // 计算力场: F_axes(k) = i × k_axes/|k|² × ρ(k)
    f2 = f1 / (K² + I0² + J0²) × k[axes];
    
    // 复数乘法: i × ρ(k) = -Im[ρ(k)] + i×Re[ρ(k)]
    FN11[index][0] = -f2 × P3D[index][1];  // 实部
    FN11[index][1] =  f2 × P3D[index][0];  // 虚部
```

**关键点**:
- 跳过 k=0 模式(零力)
- 使用转置格式存储 (y方向分解)
- 计算完成后立即逆变换回实空间

#### 6. `force_at_particle_locations()` - 网格到粒子的力插值 (行359-412)

**功能**: 使用CIC插值从网格力场得到粒子受力

**插值公式**:

```c
// 与质量分配相同的权重
f[i][axes] = Σ mesh_force[j] × W_j(x_i)
```

其中权重 `W_j` 与 `PtoMesh()` 中相同,保证动量守恒。

**实现**:

```c
// 读取8个网格点的力值,用CIC权重插值
f[i][axes] = 
    REd(fmesh, iI, J,  K ) × T3×T1×T2 +
    REd(fmesh, iI, J,  K1) × D3×T1×T2 +
    // ... 其他6个点
```

#### 7. `send_buffer_positions()` - 边界粒子交换 (行416-502)

**功能**: 交换边界区域的粒子位置,用于CIC质量分配

**为什么需要缓冲粒子?**

- x方向按进程分解,不做周期性处理
- CIC插值需要粒子周围±1个网格的信息
- 边界粒子的CIC分配需要邻居进程的网格点

**实现流程**:

```c
// 1. 选择右边界粒子 (x > x_right)
for(i = 0; i < np_local; i++) {
    if(p[i].x[0] > x_right) {
        pos[nsend] = p[i].x;  // 记录位置
        index[nsend] = i;     // 记录索引
        nsend++;
    }
}

// 2. 与左右邻居交换粒子数
MPI_Sendrecv(&nsend, RightNode, &nrecv, LeftNode);

// 3. 交换粒子位置
MPI_Sendrecv(pos, RightNode, pos_recv, LeftNode);

// 4. 将接收的粒子添加到本地粒子数组末尾
for(i = 0; i < nrecv; i++)
    p[np_local + i].x = pos_recv[i];
```

**周期性处理**: 如果是盒子的最右边界,位置需要减去 `BoxSize`

#### 8. `add_buffer_forces()` - 合并边界粒子力 (行505-532)

**功能**: 将缓冲粒子的力发送回原进程

```c
// 1. 发送缓冲粒子的力到左邻居,接收右邻居发来的力
MPI_Sendrecv(fsend, LeftNode, frecv, RightNode);

// 2. 将接收的力加到对应粒子上
for(i = 0; i < nrecv; i++)
    force[index[i]] += frecv[i];
```

### 关键技术特点

#### 1. MPI并行策略

**数据分解**:
- **实空间**: x方向分解 (`Local_nx`, `Local_x_start`)
- **傅里叶空间**: y方向分解 (`Local_ny_td`, `Local_y_start_td`)

**转置优化**:
```c
FFTW_MPI_TRANSPOSED_OUT  // 正变换输出转置
FFTW_MPI_TRANSPOSED_IN   // 逆变换输入转置
```

优点:
- 避免额外的全局转置通信
- 傅里叶空间计算按y分解,缓存友好

**边界通信**:
- 只交换x方向边界粒子
- 使用左右邻居进程 (`LeftNode`, `RightNode`)

#### 2. CIC插值的对称性

**质量分配和力插值使用相同权重**,保证:
- 动量守恒
- 能量误差最小
- 力的平滑性

#### 3. OpenMP并行

```c
#pragma omp parallel for
for(i = 0; i < np; i++) {
    // 粒子循环
}

#pragma omp atomic
d[index] += f;  // 网格更新的原子操作
```

#### 4. 内存优化

**In-place FFT**: `fftdata` 同时存储实空间和傅里叶空间数据

**内存布局**:
```c
// 实空间: Ngrid × Ngrid × Ngrid
// 傅里叶空间: Ngrid × Ngrid × (Ngrid/2+1)
// 实际分配: Ngrid × Ngrid × 2×(Ngrid/2+1)  (对齐)
```

### PM方法的物理原理

#### 泊松方程求解

在宇宙学模拟中,引力势满足泊松方程:

```
∇²Φ = 4πGρ
```

在傅里叶空间:

```
Φ(k) = -4πG ρ(k)/|k|²
```

#### 力的计算

引力是势能的梯度:

```
F = -∇Φ
```

在傅里叶空间:

```
F_i(k) = -i k_i Φ(k) = i k_i × 4πG ρ(k)/|k|²
```

#### PM方法的优缺点

**优点**:
- 计算复杂度: O(N log N) (FFT)
- 长程力准确
- 易于并行化

**缺点**:
- 空间分辨率受网格限制
- 短程力不准确(软化长度 ~ 网格间距)
- 无法处理小尺度结构

### 与COLA方法的集成

PM模块为COLA时间积分提供**残余引力**:

```
F_total = F_LPT + F_PM_residual

其中:
- F_LPT: LPT预测的引力(由位移场计算)
- F_PM_residual: PM计算的残余引力(小尺度非线性效应)
```

**COLA时间积分** (在 `cola.c` 中):
```c
// Kick-Drift-Kick (KDK) 积分
cola_kick();   // 速度更新: v += F_PM × Δt/2
cola_drift();  // 位置更新: x += v × Δt
cola_kick();   // 速度更新: v += F_PM × Δt/2
```

### 性能分析

#### 计时统计 (使用timer模块)

```c
timer_start(assign);    // 质量分配
timer_start(fft);       // 傅里叶变换
timer_start(force_mesh); // 力场计算
timer_start(pforce);    // 力插值
timer_start(comm);      // MPI通信
```

#### 内存使用

对于 `Ngrid = nc_pm`:

| 数组 | 大小 |
|------|------|
| `fftdata` | `2 × nc_pm³` 字节 |
| `density_k` | `~2 × nc_pm³` 字节 |
| 缓冲区 | `~6 × nbuf` 字节 |
| **总计** | **~4 × nc_pm³ 字节** |

**示例**: nc_pm=1024 → 约 4 GB

### 待分析项目更新

- [x] PM方法内存分配分析
- [x] **2LPT位移计算方法** ✓
- [x] **PM力计算算法** ✓
- [x] **COLA积分方案** ✓
- [ ] 各模块内部函数实现细节
- [ ] FOF晕寻找算法
- [ ] 光锥输出逻辑细节

---

## COLA模块深度分析 (cola.c)

### 模块概述

`cola.c` 实现了**COLA(COmoving Lagrangian Acceleration)方法的时间积分**,这是整个模拟的核心动力学演化模块。COLA方法通过结合LPT(拉格朗日扰动理论)的长程引力和PM(粒子网格)方法的短程残余引力,在粗网格上实现高效的N体模拟。

**代码来源**: 基于Svetlin Tassev的原始COLA代码 (Harvard/Princeton),采用GNU GPL v3许可

### 核心思想

COLA方法的基本公式:

```
x(a) = x_LPT(a) + x_residual(a)

其中:
- x_LPT(a) = q + D⁺(a)×ψ⁽¹⁾(q) + D₂(a)×ψ⁽²⁾(q)  // LPT预测轨迹
- x_residual(a): PM方法计算的残余位移
```

**关键优势**:
- LPT在长波极限下精确,可在大步长下使用
- PM只计算短程残余力,可在粗网格上运行
- 总体计算量大幅减少,同时保持小尺度精度

### 核心数据结构

```c
// 全局参数
static float Om = -1.0f;                    // 物质密度参数(运行时设置)
static const float subtractLPT = 1.0f;      // LPT扣除开关(1=启用)
static const float nLPT = -2.5f;            // LPT时间依赖指数
static const int fullT = 1;                 // 速度增长模型选择

const int gb_use_numerical_q1q2 = 1;        // 使用数值计算的q1,q2系数
```

### 主要函数功能

#### 1. `cola_kick()` - 速度更新 (行77-141)

**功能**: Leap-frog积分的Kick步骤,更新粒子速度

**参数**:
- `particles`: 粒子集合
- `Omega_m`: 物质密度参数
- `avel1`: 新的速度尺度因子 (t + 0.5*dt)

**算法流程**:

```c
// 1. 设置时间参数
AI = particles->a_v;   // 当前速度尺度因子 (t - 0.5*dt)
A  = particles->a_x;   // 当前位置尺度因子 (t)
AF = avel1;            // 新的速度尺度因子 (t + 0.5*dt)

// 2. 计算宇宙学修正因子
Om143 = [Ω_m(A) / (Ω_m + (1-Ω_m)×A³×A^{-3(1+w)})]^{1/143}

// 3. 计算LPT系数 (两种模式)
if (gb_use_numerical_q1q2) {
    // 使用数值求解的增长因子
    q1 = 1.5×Om×growthD(A)
    q2 = 1.5×Om×growthD(A)²×(1 + 7/3×Om143)
    // 通过gb_Tsquare_DE进行数值修正
} else {
    // 使用解析近似
    q1 = 1.5×Om×growthD(A)
    q2 = 1.5×Om×growthD(A)²×(1 + 7/3×Om143)
}

// 4. 计算时间积分因子
dda = Sphi(AI, AF, A)  // ∫_{AI}^{AF} gpQ(a)/Q(a) da

// 5. 更新每个粒子的速度
for each particle i:
    // 总加速度 = PM力 + LPT修正力
    a_x = -1.5×Om×f[i][0] - subtractLPT×(P[i].dx1[0]×q1 + P[i].dx2[0]×q2)
    a_y = -1.5×Om×f[i][1] - subtractLPT×(P[i].dx1[1]×q1 + P[i].dx2[1]×q2)
    a_z = -1.5×Om×f[i][2] - subtractLPT×(P[i].dx1[2]×q1 + P[i].dx2[2]×q2)
    
    // Kick: v(t+0.5dt) = v(t-0.5dt) + a × dda
    P[i].v[0] += a_x × dda
    P[i].v[1] += a_y × dda
    P[i].v[2] += a_z × dda

// 6. 更新速度尺度因子
particles->a_v = AF
```

**物理含义**:

- **PM力项** (`-1.5×Om×f[i]`): 来自PM计算的残余引力
- **LPT力项** (`-subtractLPT×(dx1×q1 + dx2×q2)`): LPT预测的长程引力
- **符号**: `subtractLPT=1.0` 表示从总力中扣除LPT部分,只保留残余力

**时间积分因子 `dda`**:

```
dda = Sphi(AI, AF, A) = [gpQ(AF) - gpQ(AI)] × A / Q(A) / DERgpQ(A)

其中:
- gpQ(a) = a^{nLPT} = a^{-2.5}  (LPT时间依赖)
- Q(a) = a³H(a)/H₀  (无量纲Hubble参数)
- DERgpQ(a) = nLPT × a^{nLPT-1}  (gpQ的导数)
```

这是一个**近似的时间积分**,基于LPT的时间依赖形式。

#### 2. `cola_drift()` - 位置更新 (行143-177)

**功能**: Leap-frog积分的Drift步骤,更新粒子位置

**参数**:
- `particles`: 粒子集合
- `Omega_m`: 物质密度参数
- `apos1`: 新的位置尺度因子 (t + dt)

**算法流程**:

```c
// 1. 设置时间参数
A  = particles->a_x;   // 当前位置尺度因子 (t)
AC = particles->a_v;   // 当前速度尺度因子 (t + 0.5dt)
AF = apos1;            // 新的位置尺度因子 (t + dt)

// 2. 计算时间积分因子
dyyy = Sq(A, AF, AC)   // ∫_{A}^{AF} gpQ(a)/Q(a) da / gpQ(AC)

// 3. 计算增长因子的变化
da1 = growthD(AF) - growthD(A)    // ΔD⁺ (一阶增长因子变化)
da2 = growthD2(AF) - growthD2(A)  // ΔD₂ (二阶增长因子变化)

// 4. 更新每个粒子的位置
for each particle i:
    // 位置更新 = 漂移项 + LPT修正项
    x[0] += v[0]×dyyy + subtractLPT×(dx1[0]×da1 + dx2[0]×da2)
    x[1] += v[1]×dyyy + subtractLPT×(dx1[1]×da1 + dx2[1]×da2)
    x[2] += v[2]×dyyy + subtractLPT×(dx1[2]×da1 + dx2[2]×da2)

// 5. 更新位置尺度因子
particles->a_x = AF
```

**物理含义**:

- **漂移项** (`v×dyyy`): 粒子由于速度的自由流动
- **LPT修正项** (`dx1×da1 + dx2×da2`): LPT预测的位移增量

**时间积分因子 `dyyy`**:

```
dyyy = Sq(A, AF, AC) = ∫_{A}^{AF} [gpQ(a)/Q(a)] da / gpQ(AC)

这个积分通过GSL数值积分计算 (见下文)
```

#### 3. 增长因子函数

##### `growthD(a)` - 一阶增长因子 (行220-229)

**功能**: 计算线性增长因子 D⁺(a)

**两种模式**:

```c
if (gb_use_solve_growth) {
    // 使用预计算的数值解 (来自solve_growth.c)
    intpl_D_E_Qs(1/a - 1, &Da, &Ea, &QdDda, &QdEda);
    return Da;
} else {
    // 使用解析近似 (growthDtemp)
    return growthDtemp(a) / growthDtemp(1.0);
}
```

**`growthDtemp(a)`** (行180-218): 使用超几何函数的解析表达式

```c
// 对于ΛCDM宇宙学,增长因子可以解析表达为:
// D⁺(a) ∝ H(a) × ∫_0^a da' / (a'×H(a'))³

// 通过超几何函数 _2F_1 计算
x = -Ω_m / [(Ω_m - 1) × a³]

if (|x-1| < 0.001):
    // 使用泰勒展开近似 (在x≈1附近)
    hyperP = 0.8596 - 0.1017×(x-1) + 0.0258×(x-1)² - ...
    hyperM = 1.1765 + 0.1585×(x-1) - 0.0142×(x-1)² + ...
else if (x < 1):
    // 使用超几何函数
    hyperP = _2F_1(1/2, 2/3, 5/3, -x)
    hyperM = _2F_1(-1/2, 2/3, 5/3, -x)
else:
    // 变量变换后的表达式
    ...

// 最终结果
if (a > 0.2):
    return √[1 + (a^{-3}-1)×Ω_m] × 
           [3.449×(1/Ω_m-1)^{2/3} + (hyperP×(4a³(Ω_m-1)-Ω_m) 
           - 7a³×hyperM×(Ω_m-1)) / (a⁵(Ω_m-1) - a²Ω_m)]
else:
    // 早期宇宙的级数展开
    return (a×(1-Ω_m)^{1.5} × (...)) / (1.5625e10 × Ω_m⁵)
```

##### `growthD2(a)` - 二阶增长因子 (行247-254)

**功能**: 计算二阶增长因子 D₂(a)

```c
// 近似公式: D₂(a) ∝ D⁺(a)² × Ω(a)^{-1/143}
// 其中 Ω(a) = Ω_m / [Ω_m + (1-Ω_m)×a³×a^{-3(1+w)}]

if (gb_use_solve_growth):
    return Ea;  // 来自solve_growth.c的数值解
else:
    return growthD2temp(a) / growthD2temp(1.0)
```

##### `growthD2v(a)` - 二阶增长因子的速度项 (行257-268)

**功能**: 计算 dD₂/dy,用于快照的速度插值

```c
// 理论关系: dD₂/dy = Q(a) × (D₂/a) × 2 × Ω(a)^{6/11}

if (gb_use_solve_growth):
    return QdEda;  // 来自solve_growth.c的数值导数
else:
    return Qfactor(a) × (d2/a) × 2 × Ω(a)^{6/11}
```

##### `decayD(a)` - 衰减模式 (行270-273)

**功能**: 计算衰减模式 D₋(a)

```c
D₋(a) = √[Ω_m/a³ + (1-Ω_m)×a^{-3(1+w)}]
```

在标准COLA中不使用衰减模式,但定义用于完整性。

#### 4. 时间积分辅助函数

##### `Qfactor(a)` - 无量纲Hubble参数 (行232-235)

```c
Q(a) = a³ × H(a) / H₀
     = √[Ω_m/a³ + (1-Ω_m)×a^{-3(1+w)}] × a³
```

其中 `w` 是暗能量状态方程参数 (通过全局变量 `gb_w` 设置)。

##### `Sq(ai, af, aRef)` - 漂移时间积分 (行315-334)

**功能**: 计算位置更新的时间积分因子

```c
Sq(ai, af, aRef) = ∫_{ai}^{af} [gpQ(a)/Q(a)] da / gpQ(aRef)

其中:
- gpQ(a) = a^{nLPT} = a^{-2.5}
- Q(a) = Qfactor(a)

使用GSL自适应积分计算:
- 工作空间: 5000个点
- 相对精度: 1e-5
- 积分规则: Gauss-Kronrod (rule 6)
```

**为什么需要这个积分?**

在COLA方法中,时间坐标变换为:

```
dy = da / [a² × H(a)]

粒子的运动方程变为:
dx/dy = v
dv/dy = -∇Φ × [a² × H(a)]⁻¹

使用LPT的时间依赖 gpQ(a) = a^{nLPT},可以将积分分解为:
∫ dy = ∫ [gpQ(a)/Q(a)] da / gpQ(a)
```

##### `Sphi(ai, af, aRef)` - Kick时间积分 (行340-345)

**功能**: 计算速度更新的时间积分因子

```c
Sphi(ai, af, aRef) = [gpQ(af) - gpQ(ai)] × aRef / Q(aRef) / DERgpQ(aRef)

其中 DERgpQ(a) = nLPT × a^{nLPT-1}

这是基于LPT时间依赖的近似积分。
```

##### `DprimeQ(a, nGrowth)` - 增长因子的导数 (行275-287)

**功能**: 计算 Q(a) × d(D⁺^n)/da

```c
DprimeQ(a, nGrowth) = Q(a) × d[D⁺(a)^n] / da

用于快照的速度插值,其中 nGrowth=1.0。
```

#### 5. 快照插值函数

##### `set_noncola_initial()` - 初始快照设置 (行350-391)

**功能**: 在非COLA运行中(如纯PM),从粒子数据生成快照

```c
// 1. 获取当前尺度因子
aout = particles->a_x

// 2. 计算速度因子
vfac = 100.0 / aout  // 转换为 km/s, H₀=100 km/s/Mpc

// 3. 计算速度增长因子
Dv = DprimeQ(aout, 1.0)   // dD⁺/dy
Dv2 = growthD2v(aout)     // dD₂/dy

// 4. 对每个粒子插值速度
for each particle:
    v[i] = vfac × (dx1×Dv + dx2×Dv2)
    x[i] = p[i].x
    id[i] = p[i].id
```

**注意**: 这个函数**不使用PM力**,只使用LPT位移。适用于纯LPT运行。

##### `cola_set_snapshot()` - COLA快照设置 (行394-485)

**功能**: 在COLA运行中,从当前粒子状态插值到指定输出时刻的快照

**为什么需要插值?**

COLA的时间积分是离散的(Kick-Drift-Kick),输出时刻可能不在积分步上。需要在两个积分步之间插值。

**算法流程**:

```c
// 1. 设置时间参数
AI = particles->a_v;   // 上次Kick的时间 (t - 0.5dt)
A  = particles->a_x;   // 上次Drift的时间 (t)
AF = aout;             // 输出时刻 (t + δt, 0 ≤ δt ≤ dt)
AC = particles->a_v;   // 当前速度时间 (t + 0.5dt)

// 2. 计算LPT系数 (与cola_kick相同)
q1 = 1.5×Om×growthD(A)
q2 = 1.5×Om×growthD(A)²×(1 + 7/3×Om143)
// 可选: 使用数值修正

// 3. 计算时间积分因子
dda = Sphi(AI, AF, A)   // 从AI到AF的Kick积分
dyyy = Sq(A, AF, AC)    // 从A到AF的Drift积分

// 4. 计算增长因子变化
da1 = growthD(AF) - growthD(A)
da2 = growthD2(AF) - growthD2(A)

// 5. 计算速度增长因子
Dv = DprimeQ(aout, 1.0)
Dv2 = growthD2v(aout)

// 6. 对每个粒子进行插值
for each particle:
    // 计算PM加速度
    ax = -1.5×Om×f[0] - (dx1[0]×q1 + dx2[0]×q2)
    ay = -1.5×Om×f[1] - (dx1[1]×q1 + dx2[1]×q2)
    az = -1.5×Om×f[2] - (dx1[2]×q1 + dx2[2]×q2)
    
    // 速度 = 当前速度 + Kick + LPT速度项
    v_out[0] = vfac × [v[0] + ax×dda + (dx1[0]×Dv + dx2[0]×Dv2)]
    v_out[1] = vfac × [v[1] + ay×dda + (dx1[1]×Dv + dx2[1]×Dv2)]
    v_out[2] = vfac × [v[2] + az×dda + (dx1[2]×Dv + dx2[2]×Dv2)]
    
    // 位置 = 当前位置 + Drift + LPT位移项
    x_out[0] = x[0] + v[0]×dyyy + (dx1[0]×da1 + dx2[0]×da2)
    x_out[1] = x[1] + v[1]×dyyy + (dx1[1]×da1 + dx2[1]×da2)
    x_out[2] = x[2] + v[2]×dyyy + (dx1[2]×da1 + dx2[2]×da2)
    
    id_out = id
```

**关键点**:

- **两步插值**: 先从 (t-0.5dt) Kick到输出时刻,再从 (t) Drift到输出时刻
- **LPT修正**: 始终加上LPT的位移和速度贡献
- **单位转换**: 速度转换为 km/s (H₀=100 km/s/Mpc)

### 关键技术特点

#### 1. Leap-frog KDK (Kick-Drift-Kick) 积分

COLA使用标准的Leap-frog KDK格式:

```
Kick:    v(t+0.5dt) = v(t-0.5dt) + a(t) × dda
Drift:   x(t+dt)   = x(t) + v(t+0.5dt) × dyyy + LPT_drift
Kick:    v(t+1.5dt) = v(t+0.5dt) + a(t+dt) × dda'
```

**优点**:
- 辛积分器,能量守恒性好
- 二阶精度
- 时间可逆

#### 2. LPT扣除机制

COLA的核心思想是**从总力中扣除LPT预测的部分**,只计算残余力:

```
F_total = F_PM + F_LPT
F_residual = F_total - F_LPT = F_PM

其中:
- F_PM: PM方法计算的总引力
- F_LPT: LPT预测的长程引力 = -(dx1×q1 + dx2×q2)
```

**为什么有效?**

- LPT在长波极限下精确,可以大步长积分
- PM只计算短程残余力,可以在粗网格上运行
- 总体计算量 ~ O(N log N),而非O(N²)

#### 3. 数值 vs 解析增长因子

代码提供两种模式:

```c
if (gb_use_solve_growth) {
    // 使用solve_growth.c的数值解
    // 优点: 精确,适用于任意宇宙学
    // 缺点: 需要预计算
} else {
    // 使用growthDtemp的解析近似
    // 优点: 无需预计算
    // 缺点: 仅适用于ΛCDM
}
```

类似地,`gb_use_numerical_q1q2` 控制是否使用数值计算的q1,q2系数。

#### 4. 暗能量状态方程

代码支持一般暗能量状态方程:

```c
// 默认: w = 0 (宇宙学常数)
// 通过全局变量 gb_w 设置

Ω(a) = Ω_m / [Ω_m + (1-Ω_m)×a³×a^{-3(1+w)}]

当 w = -1: ΛCDM
当 w = 0: 无暗能量
当 w ≠ -1: 动力学暗能量
```

#### 5. 超几何函数计算

`growthDtemp()` 使用GSL的超几何函数 `_2F_1`:

```c
// 对于ΛCDM,增长因子可以解析表达为:
// D⁺(a) ∝ H(a) × ∫_0^a da' / (a'×H(a'))³
//        ∝ H(a) × _2F_1(1/2, 2/3, 5/3, -x)

// 代码使用分段计算:
// - x≈1: 泰勒展开 (避免超几何函数的数值不稳定)
// - x<1: 直接使用超几何函数
// - x>1: 变量变换后的表达式
// - a<0.2: 早期宇宙的级数展开
```

### 性能分析

#### 计算复杂度

每个时间步:

| 操作 | 复杂度 | 备注 |
|------|--------|------|
| Kick | O(N) | 粒子循环,OpenMP并行 |
| Drift | O(N) | 粒子循环,OpenMP并行 |
| 增长因子计算 | O(1) | 预计算或解析 |
| 时间积分 (Sq, Sphi) | O(1) | 预计算或解析近似 |

#### 内存使用

COLA模块本身内存使用很少:

| 变量 | 大小 |
|------|------|
| 全局参数 | ~50 字节 |
| 临时变量 | ~200 字节 |
| **总计** | **< 1 KB** |

主要内存消耗在粒子数据和力场 (由其他模块管理)。

#### 计时统计

从 `timer.c` 模块可以看到,COLA的时间主要用于:

```
timer_start(evolve);   // Kick + Drift
timer_stop(evolve);
```

对于典型模拟 (N=512³, 50步):

- Kick: ~10% 总时间
- Drift: ~10% 总时间
- PM力计算: ~80% 总时间

### 与LPT和PM模块的集成

#### LPT模块接口

```c
// lpt.c 提供:
// - P[i].dx1[3]: 一阶LPT位移 (外推到a=1)
// - P[i].dx2[3]: 二阶LPT位移 (外推到a=1)

// cola.c 使用:
// - LPT力: -(dx1×q1 + dx2×q2)
// - LPT位移: dx1×da1 + dx2×da2
```

#### PM模块接口

```c
// pm.c 提供:
// - particles->force[3]: 残余引力场

// cola.c 使用:
// - PM力: -1.5×Om×f[i]
```

### 待分析项目更新

- [x] PM方法内存分配分析
- [x] **2LPT位移计算方法** ✓
- [x] **PM力计算算法** ✓
- [x] **COLA积分方案** ✓
- [ ] 各模块内部函数实现细节
- [ ] FOF晕寻找算法
- [ ] 光锥输出逻辑细节
