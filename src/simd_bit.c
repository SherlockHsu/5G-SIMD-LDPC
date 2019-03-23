#include "simd_bit.h"

#include <stdio.h>
#include <stdint.h>
#include <string.h>
#include <immintrin.h>

int fast_decide_byte_avx512(const int8_t *input, int len, uint8_t *output)
{
    int i;
    int units = len / 64;
    int remain = len % 64;
    __m512i zero_zmm = _mm512_setzero_si512();
    __m512i input_zmm;
    __mmask64 output_m64;
    const int8_t *pi = input;
    __mmask64 *po = (__mmask64 *)output;

    for (i = 0; i < units; ++i)
    {
        input_zmm = _mm512_loadu_si512(pi);
        output_m64 = _mm512_cmp_epi8_mask(input_zmm, zero_zmm, _MM_CMPINT_LT);
        _store_mask64(po, output_m64);
        pi += 64;
        ++po;
    }
    if (remain)
    {
        units = (remain - 1) / 8 + 1;
        input_zmm = _mm512_loadu_si512(pi);
        output_m64 = _mm512_cmp_epi8_mask(input_zmm, zero_zmm, _MM_CMPINT_LT);
        memcpy(po, &output_m64, units);
    }

    return 0;
}

int fast_extend_avx512(const uint8_t *input, int len, uint8_t *output)
{
    int i;
    int units = len / 8;
    int remain = len % 8;
    __m512i output_zmm;
    __mmask64 input_m64;
    __mmask64 mask_m64;
    const __mmask64 *pi = (const __mmask64 *)input;
    uint8_t *po = output;

    for (i = 0; i < units; ++i)
    {
        input_m64 = _load_mask64(pi);
        output_zmm = _mm512_maskz_set1_epi8(input_m64, 1);
        _mm512_storeu_si512(po, output_zmm);
        ++pi;
        po += 64;
    }
    if (remain)
    {
        mask_m64 = (uint64_t)1 << remain - 1;
        input_m64 = _load_mask64(pi);
        output_zmm = _mm512_maskz_set1_epi8(input_m64, 1);
        _mm512_mask_storeu_epi8(po, mask_m64, output_zmm);
    }

    return 0;
}

int	fast_decide_bit_avx512(const int8_t *input, int len, uint8_t *output)
{
    int i;
    int units = len / 64;
    int remain = len % 64;
    __m512i zero_zmm = _mm512_setzero_si512();
    __m512i input_zmm;
    __mmask64 output_m64, mask_m64;
    __m512i output_zmm;
    const int8_t *pi = input;
    uint8_t *po = output;

    for (i = 0; i < units; ++i)
    {
        input_zmm = _mm512_loadu_si512(pi);
        output_m64 = _mm512_cmp_epi8_mask(input_zmm, zero_zmm, _MM_CMPINT_LT);
        output_zmm = _mm512_maskz_set1_epi8(output_m64, 1);
        _mm512_storeu_si512(po, output_zmm);
        pi += 64;
        po += 64;
    }
    if (remain)
    {
        mask_m64 = (uint64_t)1 << remain - 1;
        input_zmm = _mm512_loadu_si512(pi);
        output_m64 = _mm512_cmp_epi8_mask(input_zmm, zero_zmm, _MM_CMPINT_LT);
        output_zmm = _mm512_maskz_set1_epi8(output_m64, 1);
        _mm512_mask_storeu_epi8(po, mask_m64, output_zmm);
    }

    return 0;
}

int fast_crc_checksum_avx512(const int8_t *input, int len, uint32_t crc_poly, int crc_order, uint8_t *output)
{
}