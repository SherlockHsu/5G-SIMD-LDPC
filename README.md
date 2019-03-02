5G-SIMD-LDPC
==============================

* Vesion:   2.2
* Date:     2019.02.26
* Author:   Xu Yi

File Specification
------------------------------

```
5G-SIMD-LDPC
│   Makefile
│   README.md
│   tb_ldpc.c           // Test main source file
│
└───inc
│   │   simd_ldpc.h     // Head file
│   
└───lib
    │   libsimd_ldpc.a  // Static library for Linux
    │   simd_ldpc.lib   // Static library for Windows
```

Build and Execution Instructions
------------------------------

### Build:
> make

### Execution:
> ./main

### Clean:
> make clean

Decoder Performance
------------------------------------------------------------

### Environment:
* OS: Ubuntu 16.04 xenial
* Kernel: x86_64 Linux 4.4.0-21-generic
* CPU: 4x Intel Xeon Gold 6154 CPU @ 3.001GHz
* ICC: 18.0.2

### Simulation Parameter:

* Base Graph: 5GNR Base Graph 1(i<sub>LS</sub> = 2)
* Code Block Length: 8448
* Code Rate: 948/1024

### Result:
| SIMD          | Throughput    | Code Block Latency    |
| ------------- | ------------- | --------------------- |
| SSE4.1        | 77.54Mbps     | 108.95μs              |
| AVX2          | 145.14Mbps    | 58.21μs               |
| AVX512        | 223.41Mbps    | 37.81μs               |
