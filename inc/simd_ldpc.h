// version 3.1
#ifndef SIMD_LDPC_H
#define SIMD_LDPC_H

#include <stdint.h>
#include <stdlib.h>
#include <immintrin.h>
#include "crc.h"

#define SIMD_MODE_AUTO 0
#define SIMD_MODE_SSE 1
#define SIMD_MODE_AVX2 2
#define SIMD_MODE_AVX512 3

#define DECODER_MODE_OMS 1
#define DECODER_MODE_NMS 2

#define EARLY_STOP_OFF 0	// early stop function off
#define EARLY_STOP_ON 1		// early stop function on

typedef struct nr5g_ldpc_simd_t
{
	/* General Parameters */
	int32_t B;		 // the length of TBS, 'B' in TS38.212
	int32_t R;		 // code rate*1024
	int32_t K_b;	 // block number of message bits, 'K_b' in TS38.212
	int32_t K_p;	 // 'K_+' in TS38.212
	int32_t K_n;	 // 'K_-' in TS38.212
	int32_t C;		 // number of code block after CBS, 'C' in TS38.212
	int32_t L;		 // CRC lenth, 'L' in TS38.212
	int32_t iLS;	 // LDPC lifting size
	int32_t BG_sel;  // number of base graph
	int32_t Z_c;	 // block size, 'Z_c' in TS38.212
	int32_t K;		 // message bits length, 'K' in TS38.212
	int32_t N;		 // matrix length, 'N' in TS38.212
	int32_t col_hbg; // column number of base graph
	int32_t row_hbg; // row number of base graph
	int32_t K_cb;	// maximum code block size, 'K_cb' in TS38.212
	int32_t E;		 // transmit block length
	int32_t k0;		 // start point of rate matching
	int32_t N_cb;	// decode bits block length
	int32_t Nd;		 // decode H matrix col num
	int32_t Md;		 // decode H matrix row num
	int32_t G;		 // total number of coded bits for transmission
	int32_t E_p;	 // 'E_+'
	int32_t E_n;	 // 'E_-'
	int16_t *H_BG;   // base graph
	int32_t simd_mode;

	/* Encoder Parameters */
	int8_t *p0, *p1, *p2, *p3;
	int8_t *cbs_bits;
	int8_t *coded_bits;
	int8_t *crc_serial;

	/* Decoder Parameters */
	int32_t M_whole;
	int32_t col_hbg_d;
	int32_t row_hbg_d;
	float *rdmed_llr;
	int8_t *decoded_bits;
	int8_t *degree;				 // number of connective check nodes(length:M/hbg_row_d)
	__m128i *cn_msg_sse;		 // sse message from cn to vn(length:M_whole)
	__m128i *vn_msg_sse;		 // temp sse message from vn to cn(length:19)
	__m256i *cn_msg_avx2;		 // avx2 message from cn to vn(length:M_whole)
	__m256i *vn_msg_avx2;		 // temp avx2 message from vn to cn(length:19)
	__m512i *cn_msg_avx512;		 // avx512 message from cn to vn(length:M_whole)
	__m512i *vn_msg_avx512;		 // temp avx512 message from vn to cn(length:19)
	int8_t *llr_fixed;			 // fixed llr info(length:2*REG_SIZE+Nd)
	int32_t units;				 // floor(Z_c/REG_SIZE)
	int32_t whole_degree;		 // sum of degree of every hbg_row_d
	int8_t **llr_addr_start;	 // llr address(length:whole_degree*units)
	int8_t *llr_addr_flag;		 // flag for access type(0:no mask;1:1 mask;2:2 mask)
	int8_t **llr_addr_pre;		 // extra llr address for flag=2(length:whole_degree)
	__m128i *mask_sse;			 // mask1 for flag=2(length:whole_degree)
	__m128i *mask_pre_sse;		 // mask2 for flag=2(length:whole_degree)
	__m128i endmask_sse;		 // mask for flag=1
	__m256i *mask_avx2;			 // mask1 for flag=2(length:whole_degree)
	__m256i *mask_pre_avx2;		 // mask2 for flag=2(length:whole_degree)
	__m256i endmask_avx2;		 // mask for flag=1
#ifdef PAST_METHOD	
	__m512i *mask_avx512;		 // mask1 for flag=2(length:whole_degree)
	__m512i *mask_pre_avx512;	 // mask2 for flag=2(length:whole_degree)
	__m512i endmask_avx512;		 // mask for flag=1
#else
	__mmask64 *mmask_avx512;	 // mmask1 for flag=2(length:whole_degree)
	__mmask64 *mmask_pre_avx512; // mmask2 for flag=2(length:whole_degree)
	__mmask64 endmmask_avx512;   // mmask for flag=1
#endif

	nr5g_crc_t *crc_t;
	uint8_t *byte_list;

} nr5g_ldpc_simd_t;

/*************************************************************************************/
/*                            Declare LDPC initial functions                         */
/*************************************************************************************/

void nr5g_ldpc_simd_init(nr5g_ldpc_simd_t *h, int32_t B, int32_t R, int32_t simd_mode);

/* initial LDPC parameter */
void nr5g_ldpc_simd_mode_init(nr5g_ldpc_simd_t *h);
void nr5g_ldpc_simd_param_init(nr5g_ldpc_simd_t *h, int32_t B, int32_t R);
void nr5g_ldpc_simd_rvid_param_init(nr5g_ldpc_simd_t *h, int32_t rvid);

/* initial base graph matrix */
void nr5g_ldpc_simd_matrix_init(nr5g_ldpc_simd_t *h);

/* initial encoder parameters */
void nr5g_ldpc_simd_encoder_mem_init(nr5g_ldpc_simd_t *h);

/* initial decoder parameters */
void nr5g_ldpc_simd_decoder_mem_init(nr5g_ldpc_simd_t *h);
void nr5g_ldpc_simd_decoder_param_init(nr5g_ldpc_simd_t *h);

/*************************************************************************************/
/*                              Declare LDPC free functions                          */
/*************************************************************************************/

void free_nr5g_ldpc_simd_t(nr5g_ldpc_simd_t *h);
void free_nr5g_ldpc_encoder(nr5g_ldpc_simd_t *h);
void free_nr5g_ldpc_decoder(nr5g_ldpc_simd_t *h);

/*************************************************************************************/
/*                     Declare LDPC code block segmentation functions                */
/*************************************************************************************/

void nr5g_ldpc_simd_cbs(const int8_t *input_bits, nr5g_ldpc_simd_t *h, int8_t *output_bits);

void nr5g_ldpc_simd_cbs_scb(const int8_t *input_bits, nr5g_ldpc_simd_t *h, int8_t *output_bits, int32_t r);

void nr5g_fec_crc_encode(const int8_t *input_bits, int32_t len, int32_t L, int8_t *output_bits);

void nr5g_ldpc_simd_decbs(const int8_t *input_bits, nr5g_ldpc_simd_t *h, int8_t *output_bits);

void nr5g_ldpc_simd_decbs_scb(const int8_t *input_bits, nr5g_ldpc_simd_t *h, int8_t *output_bits, int32_t r);

/*************************************************************************************/
/*                            Declare LDPC encode functions                          */
/*************************************************************************************/

void nr5g_ldpc_simd_encoder(const int8_t *info_bits, nr5g_ldpc_simd_t *h, int8_t *coded_bits);

void nr5g_ldpc_simd_encoder_scb(const int8_t *info_bits, nr5g_ldpc_simd_t *h, int8_t *coded_bits);

/*************************************************************************************/
/*                        Declare LDPC rate matching functions                       */
/*************************************************************************************/

void nr5g_ldpc_simd_rate_matching(const int8_t *coded_bits, nr5g_ldpc_simd_t *h, int8_t *rmed_bits);

void nr5g_ldpc_simd_rate_matching_scb(const int8_t *coded_bits, nr5g_ldpc_simd_t *h, int8_t *rmed_bits, int32_t r);

void nr5g_ldpc_simd_rate_dematching(const float *llr, nr5g_ldpc_simd_t *h, float *rdmed_llr);

void nr5g_ldpc_simd_rate_dematching_scb(const float *llr, nr5g_ldpc_simd_t *h, float *rdmed_llr, int32_t r);

/*************************************************************************************/
/*                            Declare LDPC decode functions                          */
/*************************************************************************************/

void nr5g_ldpc_simd_decoder(const float *llr, nr5g_ldpc_simd_t *h, int32_t I_max, float coef, int32_t decoder_mode, int32_t early_stop, int8_t *decoded_bits, float *decoded_llr);

void nr5g_ldpc_fast_load_llr_simd_scb(const float *llr, nr5g_ldpc_simd_t *h, int32_t r);

void nr5g_ldpc_simd_decoder_scb(nr5g_ldpc_simd_t *h, int32_t I_max, float coef, int32_t decoder_mode, int32_t early_stop, int8_t *decoded_bits, float *decoded_llr, int32_t r);

/*************************************************************************************/
/*                             Declare LDPC combo functions                          */
/*************************************************************************************/

void nr5g_ldpc_simd_cbs_enc_rm(const int8_t *info_bits, nr5g_ldpc_simd_t *h, int8_t *rmed_bits);

void nr5g_ldpc_simd_rdm_dec_decbs(const float *llr, nr5g_ldpc_simd_t *h, int32_t I_max, int32_t decoder_mode, int32_t early_stop, float coef, int8_t *decbs_bits, float *decoded_llr);

int is_ldpc_code(nr5g_ldpc_simd_t *h, uint8_t *code);

#endif // !SIMD_LDPC_H
