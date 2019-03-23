#include "crc.h"
#include <stdio.h>
#include <string.h>
#include <stdlib.h>

void nr5g_gen_crc_table(nr5g_crc_t *h)
{

    int i, j, ord = (h->order - 8);
    uint64_t bit, crc;

    for (i = 0; i < 256; i++)
    {
        crc = ((uint64_t)i) << ord;
        for (j = 0; j < 8; j++)
        {
            bit = crc & h->crchighbit;
            crc <<= 1;
            if (bit)
                crc ^= h->polynom;
        }
        h->table[i] = crc & h->crcmask;
    }
}

int nr5g_crc_set_init(nr5g_crc_t *crc_par, uint64_t crc_init_value)
{

    crc_par->crcinit = crc_init_value;
    if (crc_par->crcinit != (crc_par->crcinit & crc_par->crcmask))
    {
        printf("ERROR, invalid crcinit in crc_set_init().\n");
        return -1;
    }
    return 0;
}

int nr5g_crc_init(nr5g_crc_t *h, uint32_t crc_poly, int crc_order)
{

    // Set crc working default parameters
    h->polynom = crc_poly;
    h->order = crc_order;
    h->crcinit = 0x00000000;

    // Compute bit masks for whole CRC and CRC high bit
    h->crcmask = ((((uint64_t)1 << (h->order - 1)) - 1) << 1) | 1;
    h->crchighbit = (uint64_t)1 << (h->order - 1);

    // check parameters
    if (h->order % 8 != 0)
    {
        fprintf(stderr, "ERROR, invalid order=%d, it must be 8, 16, 24 or 32.\n",
                h->order);
        return -1;
    }

    if (nr5g_crc_set_init(h, h->crcinit))
    {
        fprintf(stderr, "Error setting CRC init word\n");
        return -1;
    }

    // generate lookup table
    nr5g_gen_crc_table(h);

    return 0;
}

// len is multiple of 8
uint32_t nr5g_crc_checksum_byte(nr5g_crc_t *h, uint8_t *data, int len)
{
    int i;
    uint32_t crc = 0;

    nr5g_crc_set_init(h, 0);

    // Calculate CRC
    for (i = 0; i < len / 8; i++)
        nr5g_crc_checksum_put_byte(h, data[i]);

    crc = (uint32_t)nr5g_crc_checksum_get(h);

    return crc;
}

uint32_t nr5g_crc_attach_byte(nr5g_crc_t *h, uint8_t *data, int len)
{
    uint32_t checksum = nr5g_crc_checksum_byte(h, data, len);

    // Add CRC
    for (int i = 0; i < h->order / 8; i++)
    {
        data[len / 8 + (h->order / 8 - i - 1)] = (checksum & (0xff << (8 * i))) >> (8 * i);
    }
    return checksum;
}

uint32_t nr5g_crc_checksum_16(nr5g_crc_t *h, uint8_t *data, int len)
{
    int i;
    uint32_t crc = 0;

    nr5g_crc_set_init(h, 0);

    uint16_t *data16 = (uint16_t *)data;

    // Calculate CRC
    for (i = 0; i < len / 16; i++)
        nr5g_crc_checksum_put_16(h, data16[i]);
        
    crc = (uint32_t)nr5g_crc_checksum_get(h);

    return crc;
}