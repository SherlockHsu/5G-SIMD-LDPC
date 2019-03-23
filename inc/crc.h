#ifndef NR5G_CRC_H
#define NR5G_CRC_H

#include <stdint.h>

typedef struct
{
    uint64_t table[256];
    int polynom;
    int order;
    uint64_t crcinit;
    uint64_t crcmask;
    uint64_t crchighbit;
    uint32_t nr5g_crc_out;
} nr5g_crc_t;

int nr5g_crc_init(nr5g_crc_t *h,
                  uint32_t nr5g_crc_poly,
                  int nr5g_crc_order);

int nr5g_crc_set_init(nr5g_crc_t *h,
                      uint64_t init_value);

uint32_t nr5g_crc_attach_byte(nr5g_crc_t *h,
                              uint8_t *data,
                              int len);

static inline void nr5g_crc_checksum_put_byte(nr5g_crc_t *h, uint8_t byte)
{

    // Polynom order 8, 16, 24 or 32 only.
    int ord = h->order - 8;
    uint64_t crc = h->crcinit;

    crc = (crc << 8) ^ h->table[((crc >> (ord)) & 0xff) ^ byte];
    h->crcinit = crc;
}

static inline uint64_t nr5g_crc_checksum_get(nr5g_crc_t *h)
{
    return (h->crcinit & h->crcmask);
}

uint32_t nr5g_crc_checksum_byte(nr5g_crc_t *h,
                                uint8_t *data,
                                int len);


uint32_t nr5g_crc_checksum_16(nr5g_crc_t *h, uint8_t *data, int len);

static inline void nr5g_crc_checksum_put_16(nr5g_crc_t *h, uint16_t byte)
{

    // Polynom order 8, 16, 24 or 32 only.
    int ord = h->order - 16;
    uint64_t crc = h->crcinit;

    crc = (crc << 16) ^ h->table[((crc >> (ord)) & 0xffff) ^ byte];
    h->crcinit = crc;
}

#endif