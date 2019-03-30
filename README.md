5G-SIMD-LDPC
==============================

* Vesion:   3.1
* Date:     2019.03.27
* Author:   Xu Yi

---
File Specification
------------------------------

```
5G-SIMD-LDPC
│   README.md
│   LICENSE
│
└───inc
│   │   bit.h                       // Data pack and unpack head file
│   │   crc.h                       // CRC head file
│   │   simd_bit.h                  // Data pack and unpack with SIMD head file
│   │   simd_ldpc.h                 // SIMD LDPC head file
│   │   thread_pool.h               // Thread pool head file
│   
└───src
│   │   thread_pool.c               // Thread pool source file
│
└───lib
│   │   libsimd_5gfec.a             // Static library for x86_64 Linux
│   
└───example
    └───test_multicore              // example for multi-core test
    │   │   test_multicore.c
    │   │   Makefile
    │    
    └───test_simd                   // example for different SIMD instruction sets
    │   │   test_simd.c
    │   │   Makefile
    │    
    └───test_stop                   // example for early stopping strategy
        │   test_stip.c
        │   Makefile
```

Build and Execution Instructions
------------------------------

### Build:
> cd example/[example name]/

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
* Code Rate: 5/6

### Result:
| SIMD          | Throughput    | Code Block Latency    |
| ------------- | ------------- | --------------------- |
| SSE4.1        | 62.34Mbps     | 135.52μs              |
| AVX2          | 130.23Mbps    | 64.87μs               |
| AVX512        | 223.29Mbps    | 37.83μs               |
