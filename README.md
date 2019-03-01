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
