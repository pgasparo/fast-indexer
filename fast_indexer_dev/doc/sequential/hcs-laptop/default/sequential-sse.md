# Runs on HCS Laptop

```
(base) |11:30|1|doc$ lscpu 
CPU(s):                          16
Liste der Online-CPU(s):         0-15
Thread(s) pro Kern:              2
Kern(e) pro Socket:              8
Sockel:                          1
NUMA-Knoten:                     1
Modellname:                      Intel(R) Core(TM) i7-10875H CPU @ 2.30GHz
CPU MHz:                         2300.000
Maximale Taktfrequenz der CPU:   5100.0000
Minimale Taktfrequenz der CPU:   800.0000
L1d Cache:                       256 KiB
L1i Cache:                       256 KiB
L2 Cache:                        2 MiB
L3 Cache:                        16 MiB
```

# Default compilation

Single threaded code, 100% CPU usage (by *top*).

```
(base) |11:30|0|doc$ ldd ../xgandalf-build/libxgandalf.so 
	linux-vdso.so.1 (0x00007ffc79dde000)
	libstdc++.so.6 => /lib/x86_64-linux-gnu/libstdc++.so.6 (0x00007fc8e0d61000)
	libm.so.6 => /lib/x86_64-linux-gnu/libm.so.6 (0x00007fc8e0c12000)
	libgcc_s.so.1 => /lib/x86_64-linux-gnu/libgcc_s.so.1 (0x00007fc8e0bf7000)
	libc.so.6 => /lib/x86_64-linux-gnu/libc.so.6 (0x00007fc8e0a05000)
	/lib64/ld-linux-x86-64.so.2 (0x00007fc8e1065000)
(base) |11:04|0|fast_indexer_dev$ ./test.sh ./test/data/image0_local.txt 100
Reading data ...
Calling xgandalf ...
Found 1 lattices
(-28.8949 20.3254 10.4158) (-28.1623 -60.4518 41.0353) (39.5732 24.9699 62.8447)
1279.29 ms per iteration
```

## perf record
```
  50.95%  xgandalf  libxgandalf.so       [.] xgandalf::function1_periodic
  12.56%  xgandalf  libxgandalf.so       [.] xgandalf::InverseSpaceTransform::performTransform
  11.00%  xgandalf  libxgandalf.so       [.] Eigen::internal::gebp_kernel<float, float, long, Eigen::internal::blas_data_mapper<float, long, 0, 0>, 8, 4, false, false>::
   4.45%  xgandalf  libm-2.31.so         [.] __roundf
   3.74%  xgandalf  libxgandalf.so       [.] xgandalf::InverseSpaceTransform::onePeriodicFunction
   3.54%  xgandalf  libxgandalf.so       [.] xgandalf::LatticeAssembler::computeCandidateLattices
   3.10%  xgandalf  libxgandalf.so       [.] Eigen::internal::gemm_pack_rhs<float, long, Eigen::internal::const_blas_data_mapper<float, long, 0>, 4, 0, false, false>::op
   1.46%  xgandalf  libxgandalf.so       [.] Eigen::internal::general_matrix_vector_product<long, float, Eigen::internal::const_blas_data_mapper<float, long, 1>, 1, fals
   1.21%  xgandalf  libxgandalf.so       [.] xgandalf::function9
   0.73%  xgandalf  libxgandalf.so       [.] Eigen::internal::call_dense_assignment_loop<Eigen::Array<float, -1, -1, 0, -1, -1>, Eigen::CwiseBinaryOp<Eigen::internal::sc
   0.65%  xgandalf  libxgandalf.so       [.] xgandalf::HillClimbingOptimizer::computeStep
   0.55%  xgandalf  libxgandalf.so       [.] xgandalf::Lattice::minimize
   0.47%  xgandalf  libxgandalf.so       [.] 0x000000000001b704
   0.46%  xgandalf  libxgandalf.so       [.] xgandalf::getGradient_detectorAngleMatch
   0.36%  xgandalf  libxgandalf.so       [.] Eigen::internal::call_dense_assignment_loop<Eigen::Matrix<double, 1, -1, 1, 1, -1>, Eigen::PartialReduxExpr<Eigen::CwiseBina
   0.30%  xgandalf  libxgandalf.so       [.] Eigen::internal::call_assignment<Eigen::Matrix<float, 3, -1, 0, 3, -1>, Eigen::Product<Eigen::MatrixWrapper<Eigen::CwiseBina
```

```
     EIGEN_DEVICE_FUNC EIGEN_STRONG_INLINE bool operator()(const LhsScalar& a, const RhsScalar& b) const {return a<b;}
  1.09 │670:   movss     (%rbx,%rax,4),%xmm1
  2.06 │       comiss    %xmm0,%xmm1
       │ _ZNK5Eigen8internal9assign_opIbbE11assignCoeffERbRKb():
  6.43 │       seta      (%rdx,%rax,1)
       │ _ZN5Eigen8internal21dense_assignment_loopINS0_31generic_dense_assignment_kernelINS0_9evaluatorINS_5ArrayIbLin1ELin1ELi0ELin1ELin1EEEEENS3_INS_13CwiseBinary
  1.79 │       add       $0x1,%rax
       │       cmp       %rax,%rbp
  1.11 │     ↑ jne 670
```

```
     _ZN8xgandalfL18function1_periodicERN5Eigen5ArrayIfLin1ELin1ELi0ELin1ELin1EEES3_S3_fRNS1_IbLin1ELin1ELi0ELin1ELin1EEEf():
       │     return (__m128) ((__v4sf)__A + (__v4sf)__B);
  0.05 │       addps     0x50(%rsp),%xmm1
       │     return (__m128) ((__v4sf)__A * (__v4sf)__B);
  0.38 │       mulps     %xmm2,%xmm1
       │     return (__m128) ((__v4sf)__A + (__v4sf)__B);
  1.66 │       addps     0x60(%rsp),%xmm1
       │     return (__m128) ((__v4sf)__A * (__v4sf)__B);
  4.56 │       mulps     %xmm2,%xmm1
  5.51 │       mulps     %xmm15,%xmm1
       │     return (__m128) ((__v4sf)__A + (__v4sf)__B);
  5.91 │       addps     %xmm15,%xmm1
       │     return (__m128) ((__v4sf)__A * (__v4sf)__B);
  0.00 │       movaps    0x70(%rsp),%xmm15
       │     _mm_mul_ps():
  0.01 │       mulps     %xmm2,%xmm15
       │     _ZN8xgandalfL18function1_periodicERN5Eigen5ArrayIfLin1ELin1ELi0ELin1ELin1EEES3_S3_fRNS1_IbLin1ELin1ELi0ELin1ELin1EEEf():
       │     return (__m128) ((__v4sf)__A + (__v4sf)__B);  ◆
       │       addps     0x80(%rsp),%xmm15
       │     return __builtin_ia32_andps (__A, __B);
  1.59 │       andps     %xmm0,%xmm1
       │     return (__m128) ((__v4sf)__A * (__v4sf)__B);
  0.00 │       mulps     %xmm2,%xmm15
       │     return (__m128) ((__v4sf)__A + (__v4sf)__B);
  0.01 │       addps     %xmm6,%xmm15
       │     return (__m128) ((__v4sf)__A * (__v4sf)__B);
       │       mulps     %xmm2,%xmm15
  2.41 │       mulps     %xmm2,%xmm15
  0.00 │       mulps     %xmm5,%xmm2
       │     return (__m128) ((__v4sf)__A - (__v4sf)__B);
  0.01 │       subps     %xmm2,%xmm15
       │     return (__m128) ((__v4sf)__A + (__v4sf)__B);
  5.22 │       addps     %xmm4,%xmm15
       │     return __builtin_ia32_andnps (__A, __B);
  1.83 │       andnps    %xmm15,%xmm0 
       │     return __builtin_ia32_orps (__A, __B);
  0.84 │       orps      %xmm1,%xmm0
       │     return __builtin_ia32_xorps (__A, __B);
  0.39 │       xorps     %xmm3,%xmm0
  1.23 │       xorps     %xmm13,%xmm0
       │     *(__m128 *)__P = __A;
  1.53 │       movaps    %xmm0,0x0(%r13,%rax,4)
```

```
--------------------------------------------------------------------------------
Group 1: MEM_SP
+------------------------------------------+---------+--------------+
|                   Event                  | Counter |  HWThread 7  |
+------------------------------------------+---------+--------------+
|             INSTR_RETIRED_ANY            |  FIXC0  | 769602644697 |
|           CPU_CLK_UNHALTED_CORE          |  FIXC1  | 306938509666 |
|           CPU_CLK_UNHALTED_REF           |  FIXC2  | 283578632160 |
|              PWR_PKG_ENERGY              |   PWR0  |     786.3879 |
|              PWR_DRAM_ENERGY             |   PWR3  |     142.2798 |
| FP_ARITH_INST_RETIRED_128B_PACKED_SINGLE |   PMC0  | 176841417550 |
|    FP_ARITH_INST_RETIRED_SCALAR_SINGLE   |   PMC1  |  15671340907 |
| FP_ARITH_INST_RETIRED_256B_PACKED_SINGLE |   PMC2  |            1 |
|                DRAM_READS                | MBOX0C1 |    170161045 |
|                DRAM_WRITES               | MBOX0C2 |    109200936 |
+------------------------------------------+---------+--------------+

+-----------------------------------+--------------+
|               Metric              |  HWThread 7  |
+-----------------------------------+--------------+
|        Runtime (RDTSC) [s]        |     124.6938 |
|        Runtime unhalted [s]       |     133.2267 |
|            Clock [MHz]            |    2493.6638 |
|                CPI                |       0.3988 |
|             Energy [J]            |     786.3879 |
|             Power [W]             |       6.3066 |
|          Energy DRAM [J]          |     142.2798 |
|           Power DRAM [W]          |       1.1410 |
|            SP [MFLOP/s]           |    5798.4982 |
|          AVX SP [MFLOP/s]         | 6.415714e-08 |
|          Packed [MUOPS/s]         |    1418.2049 |
|          Scalar [MUOPS/s]         |     125.6785 |
|  Memory load bandwidth [MBytes/s] |      87.3364 |
|  Memory load data volume [GBytes] |      10.8903 |
| Memory evict bandwidth [MBytes/s] |      56.0482 |
| Memory evict data volume [GBytes] |       6.9889 |
|    Memory bandwidth [MBytes/s]    |     143.3845 |
|    Memory data volume [GBytes]    |      17.8792 |
|       Operational intensity       |      40.4402 |
+-----------------------------------+--------------+
```

```
--------------------------------------------------------------------------------
Group 1: FLOPS_SP
+------------------------------------------+---------+--------------+
|                   Event                  | Counter |  HWThread 7  |
+------------------------------------------+---------+--------------+
|             INSTR_RETIRED_ANY            |  FIXC0  | 769602644661 |
|           CPU_CLK_UNHALTED_CORE          |  FIXC1  | 306893612244 |
|           CPU_CLK_UNHALTED_REF           |  FIXC2  | 283533062688 |
| FP_ARITH_INST_RETIRED_128B_PACKED_SINGLE |   PMC0  | 176841417550 |
|    FP_ARITH_INST_RETIRED_SCALAR_SINGLE   |   PMC1  |  15671340907 |
| FP_ARITH_INST_RETIRED_256B_PACKED_SINGLE |   PMC2  |            1 |
+------------------------------------------+---------+--------------+

+----------------------+--------------+
|        Metric        |  HWThread 7  |
+----------------------+--------------+
|  Runtime (RDTSC) [s] |     124.6492 |
| Runtime unhalted [s] |     133.2394 |
|      Clock [MHz]     |    2493.0970 |
|          CPI         |       0.3988 |
|     SP [MFLOP/s]     |    5800.5741 |
|   AVX SP [MFLOP/s]   | 6.418011e-08 |
|   Packed [MUOPS/s]   |    1418.7126 |
|   Scalar [MUOPS/s]   |     125.7235 |
|  Vectorization ratio |      91.8596 |
+----------------------+--------------+
```

```
--------------------------------------------------------------------------------
Group 1: CYCLE_STALLS
+-----------------------------------+---------+--------------+
|               Event               | Counter |  HWThread 7  |
+-----------------------------------+---------+--------------+
|         INSTR_RETIRED_ANY         |  FIXC0  | 769602644681 |
|       CPU_CLK_UNHALTED_CORE       |  FIXC1  | 307322114967 |
|        CPU_CLK_UNHALTED_REF       |  FIXC2  | 283890870528 |
|  CYCLE_ACTIVITY_STALLS_L2_PENDING |   PMC0  |   2208142657 |
| CYCLE_ACTIVITY_STALLS_LDM_PENDING |   PMC1  |  12515123887 |
| CYCLE_ACTIVITY_STALLS_L1D_PENDING |   PMC2  |   2842802977 |
|    CYCLE_ACTIVITY_STALLS_TOTAL    |   PMC3  |  19182462474 |
+-----------------------------------+---------+--------------+

+----------------------------------------+-------------+
|                 Metric                 |  HWThread 7 |
+----------------------------------------+-------------+
|           Runtime (RDTSC) [s]          |    124.8074 |
|          Runtime unhalted [s]          |    133.3851 |
|               Clock [MHz]              |   2494.1857 |
|                   CPI                  |      0.3993 |
|         Total execution stalls         | 19182462474 |
|     Stalls caused by L1D misses [%]    |     14.8198 |
|     Stalls caused by L2 misses [%]     |     11.5113 |
|    Stalls caused by memory loads [%]   |     65.2425 |
|        Execution stall rate [%]        |      6.2418 |
|  Stalls caused by L1D misses rate [%]  |      0.9250 |
|   Stalls caused by L2 misses rate [%]  |      0.7185 |
| Stalls caused by memory loads rate [%] |      4.0723 |
+----------------------------------------+-------------+
```
