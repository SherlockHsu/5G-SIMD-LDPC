#ifndef NR5G_BIT_H
#define NR5G_BIT_H

#include <stdint.h>

void nr5g_bit_unpack_vector(uint8_t *packed,
                              uint8_t *unpacked,
                              int nof_bits);

void nr5g_bit_pack_vector(uint8_t *unpacked,
                            uint8_t *packed,
                            int nof_bits);

void nr5g_bit_unpack(uint32_t value,
                       uint8_t **bits,
                       int nof_bits);

uint32_t nr5g_bit_pack(uint8_t **bits,
                         int nof_bits);

#endif