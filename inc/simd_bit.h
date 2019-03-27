#ifndef SIMD_BIT
#define SIMD_BIT

#include <stdint.h>
#include <immintrin.h>

int fast_decide_byte_avx512(const int8_t *input, int len, uint8_t *output);

int fast_extend_avx512(const uint8_t *input, int len, uint8_t *output);

int	fast_decide_bit_avx512(const int8_t *input, int len, uint8_t *output);

int	fast_decide_bit_avx2(const int8_t *input, int len, uint8_t *output);

int	fast_decide_bit_sse(const int8_t *input, int len, uint8_t *output);

#endif