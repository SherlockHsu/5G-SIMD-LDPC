// version 3.1
#include "simd_ldpc.h"
#include "simd_bit.h"
#include "crc.h"
#include "bit.h"
#include "thread_pool.h"

#ifndef _GNU_SOURCE
#define _GNU_SOURCE
#include <sched.h>
#endif

#include <stdlib.h>
#include <stdint.h>
#include <string.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <mmintrin.h>
#include <mkl.h>

#if defined(_MSC_VER)
#include <windows.h>
#else
#include <sys/time.h>
#endif

#define CORE_NUM 18
#define TEST_CORE_NUM 18

#define CRC_24A 0x1864CFB
#define BLOCK_SIZE 10000

pthread_mutex_t mutex;
pthread_mutex_t demutex;
volatile int cnt;
double decode_run_time;

typedef struct ldpc_decoder_thrd_t
{
	const float *llr;
	nr5g_ldpc_simd_t *h;
	int32_t I_max;
	float coef;
	int32_t decoder_mode;
	int8_t *decoded_bits;
	float *decoded_llr;
} ldpc_decoder_thrd_t;

void ldpc_decoder_thrd(void *arg)
{
#if defined(_MSC_VER)
	LARGE_INTEGER num;
	long long start, end, freq;
#else
	struct timeval start, end;
	long timeuse;
#endif
	ldpc_decoder_thrd_t *h = (ldpc_decoder_thrd_t *)arg;

#if defined(_MSC_VER)
	QueryPerformanceFrequency(&num);
	freq = num.QuadPart;
	QueryPerformanceCounter(&num);
	start = num.QuadPart;
#else
	gettimeofday(&start, NULL);
#endif
	for (int i = 0; i < BLOCK_SIZE; ++i)
		nr5g_ldpc_simd_decoder(h->llr, h->h, h->I_max, h->coef, h->decoder_mode, EARLY_STOP_OFF, h->decoded_bits, h->decoded_llr);
#if defined(_MSC_VER)
	QueryPerformanceCounter(&num);
	end = num.QuadPart;
	pthread_mutex_lock(&demutex);
	decode_run_time += (double)(end - start) / freq;
	pthread_mutex_unlock(&demutex);
#else
	gettimeofday(&end, NULL);
	timeuse = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
	pthread_mutex_lock(&demutex);
	decode_run_time += (double)timeuse / 1000000.0;
	pthread_mutex_unlock(&demutex);
#endif
	pthread_mutex_lock(&mutex);
	cnt++;
	pthread_mutex_unlock(&mutex);
}

#define B_NUM 2
#define R_NUM 3
#define EBN0_SIZE 6

int main()
{
	int B_list[B_NUM] = {8448, 3840};
	int R_list[R_NUM] = {853, 768, 512};
	int core_list[TEST_CORE_NUM] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18};
	int j, k;
	FILE *fp, *fq;
	fp = fopen("mutlicore_tp.txt", "a");
	fq = fopen("mutlicore_latency.txt", "a");
	for (int indx_core = 0; indx_core < TEST_CORE_NUM; indx_core++)
	{
		printf("==================================================\n");
		printf("CORE NUM: %d\n", core_list[indx_core]);
		printf("--------------------------------------------------\n");
		for (int j = 0; j < B_NUM; j++)
			for (int k = 0; k < R_NUM; k++)
			{
				double encode_run_time;
#if defined(_MSC_VER)
				LARGE_INTEGER num;
				long long start, end, freq;
#else
				struct timeval start, end;
				long timeuse;
#endif
				double avg_tp, avg_latency;

				int32_t i, indx_block, indx_ebn0, sum_err_bits, test_size;
				int32_t B, R, I_max, decoder_mode;
				float coef;
				nr5g_ldpc_simd_t *ldpc_arg[CORE_NUM];
				VSLStreamStatePtr stream_g;
				VSLStreamStatePtr stream_b;
				int8_t *info_bits[CORE_NUM];
				int8_t *info_byte[CORE_NUM];
				int8_t *rmed_bits[CORE_NUM];
				float *mapped_sig[CORE_NUM];
				float *noise[CORE_NUM];
				float *llr[CORE_NUM];
				float *decoded_llr[CORE_NUM];
				int8_t *decbs_bits[CORE_NUM];
				float EbN0, sigma2, sigma;
				int32_t err_bits[BLOCK_SIZE];

				ldpc_decoder_thrd_t *ldpct[CORE_NUM];
				pool_init(0, core_list[indx_core], 0);

				nr5g_crc_t crc_t;
				nr5g_crc_init(&crc_t, CRC_24A, 24);

				/* set parameters */
				B = B_list[j];
				R = R_list[k];
				I_max = 10;
				decoder_mode = DECODER_MODE_OMS;

				switch (decoder_mode)
				{
				case DECODER_MODE_OMS:
					coef = (float)2 / 4; // beta for oms
					break;
				case DECODER_MODE_NMS:
					coef = (float)25 / 32; // alpha for nms
					break;
				default:
					printf("ERROR(main): SIMD MODE %d IS NOT EXISTED.\n", decoder_mode);
					exit(0);
					break;
				}

				float EbN0_list[EBN0_SIZE] = {0.0, 1.0, 2.0, 3.0, 4.0, 5.0};
				test_size = EBN0_SIZE;
				pthread_mutex_init(&demutex, NULL);
				pthread_mutex_init(&mutex, NULL);
				cnt = 0;

				/* initialize */
				for (int c = 0; c < core_list[indx_core]; ++c)
				{
					ldpc_arg[c] = (nr5g_ldpc_simd_t *)malloc(sizeof(nr5g_ldpc_simd_t));
					nr5g_ldpc_simd_init(ldpc_arg[c], B, R, SIMD_MODE_AUTO);
					ldpc_arg[c]->crc_t = &crc_t;

					info_bits[c] = (int8_t *)malloc(sizeof(int8_t) * ldpc_arg[c]->B);
					info_byte[c] = (int8_t *)malloc(sizeof(int8_t) * ldpc_arg[c]->B / 8);
					rmed_bits[c] = (int8_t *)malloc(sizeof(int8_t) * ldpc_arg[c]->G);
					mapped_sig[c] = (float *)malloc(sizeof(float) * ldpc_arg[c]->G);
					noise[c] = (float *)malloc(sizeof(float) * ldpc_arg[c]->G);
					llr[c] = (float *)malloc(sizeof(float) * ldpc_arg[c]->G);
					decoded_llr[c] = (float *)malloc(sizeof(float) * ldpc_arg[c]->G);
					decbs_bits[c] = (int8_t *)malloc(sizeof(int8_t) * ldpc_arg[c]->B);

					ldpct[c] = (ldpc_decoder_thrd_t *)malloc(sizeof(ldpc_decoder_thrd_t));
					ldpct[c]->llr = ldpc_arg[c]->rdmed_llr;
					ldpct[c]->h = ldpc_arg[c];
					ldpct[c]->I_max = I_max;
					ldpct[c]->coef = coef;
					ldpct[c]->decoder_mode = decoder_mode;
					ldpct[c]->decoded_bits = ldpc_arg[c]->decoded_bits;
					ldpct[c]->decoded_llr = decoded_llr[c];
				}

				vslNewStream(&stream_g, VSL_BRNG_MCG31, 0);
				vslNewStream(&stream_b, VSL_BRNG_MCG31, 1);

				avg_tp = 0.0;
				avg_latency = 0.0;

				/* test loop start */
				for (indx_ebn0 = 0; indx_ebn0 < test_size; indx_ebn0++)
				{
					EbN0 = EbN0_list[indx_ebn0];
					// sigma2 = (float)(1 / (pow(10, (double)EbN0 / 10) * 2 * R / 1024));
					sigma2 = (float)(1 / pow(10, (double)EbN0 / 10));
					sigma = (float)sqrt(sigma2);
					sum_err_bits = 0;
					encode_run_time = 0.0;
					decode_run_time = 0.0;

					for (indx_block = 0; indx_block < 1; indx_block++)
					{
						err_bits[indx_block] = 0;

						/* generate random tbs */
						for (int c = 0; c < core_list[indx_core]; ++c)
						{
							/* generate random tbs */
							viRngUniformBits(VSL_RNG_METHOD_UNIFORMBITS_STD, stream_b, (B - 1) / 32 + 1, (uint32_t *)info_byte[c]);
							nr5g_crc_attach_byte(&crc_t, info_byte[c], B - 24);
							nr5g_bit_unpack_vector(info_byte[c], info_bits[c], B);

							/* cbs */
							nr5g_ldpc_simd_cbs(info_bits[c], ldpc_arg[c], ldpc_arg[c]->cbs_bits);

							/* encode */
#if defined(_MSC_VER)
							QueryPerformanceFrequency(&num);
							freq = num.QuadPart;
							QueryPerformanceCounter(&num);
							start = num.QuadPart;
#else
							gettimeofday(&start, NULL);
#endif
							nr5g_ldpc_simd_encoder(ldpc_arg[c]->cbs_bits, ldpc_arg[c], ldpc_arg[c]->coded_bits);
#if defined(_MSC_VER)
							QueryPerformanceCounter(&num);
							end = num.QuadPart;
							encode_run_time += (double)(end - start) / freq;
#else
							gettimeofday(&end, NULL);
							timeuse = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
							encode_run_time += (double)timeuse / 1000000.0;
#endif

							/* rate matching */
							nr5g_ldpc_simd_rate_matching(ldpc_arg[c]->coded_bits, ldpc_arg[c], rmed_bits[c]);

							/* BPSK map */
							for (i = 0; i < ldpc_arg[c]->G; i++)
								mapped_sig[c][i] = (float)2 * rmed_bits[c][i] - 1;

							/* pass AWGN channel */
							vsRngGaussian(VSL_RNG_METHOD_GAUSSIANMV_BOXMULLER, stream_g, ldpc_arg[c]->G, noise[c], 0.0, sigma);
							for (i = 0; i < ldpc_arg[c]->G; i++)
								llr[c][i] = 2 * (mapped_sig[c][i] + noise[c][i]) / sigma2;

							/* rate dematching */
							nr5g_ldpc_simd_rate_dematching(llr[c], ldpc_arg[c], ldpc_arg[c]->rdmed_llr);
						}

						/* decode */
						for (int c = 0; c < core_list[indx_core]; ++c)
							pool_add_task(ldpc_decoder_thrd, (void *)ldpct[c], 0);
						while (cnt < core_list[indx_core])
							;
						cnt = 0;

						for (int c = 0; c < core_list[indx_core]; ++c)
						{
							/* decbs */
							nr5g_ldpc_simd_decbs(ldpc_arg[c]->decoded_bits, ldpc_arg[c], decbs_bits[c]);

							/* statistics */
							for (i = 0; i < ldpc_arg[c]->B; i++)
								err_bits[indx_block] += (info_bits[c][i] == decbs_bits[c][i] ? 0 : 1);

							sum_err_bits += err_bits[indx_block];
						}
					}

					/* print results */
					avg_tp += (double)B * BLOCK_SIZE * core_list[indx_core] / (decode_run_time / core_list[indx_core]) / 1e6;
					avg_latency += decode_run_time / core_list[indx_core] * 1e6 / BLOCK_SIZE;
				}

				avg_tp /= test_size;
				avg_latency /= test_size;
				printf("B = %d, R = %d\n", B, R);
				printf("Average Throughput:\t%.2lfMbps\n", avg_tp);
				printf("Average Latency:\t%.2lfus\n", avg_latency);

				fprintf(fp, "%.2lf\t", avg_tp);
				fprintf(fq, "%.2lf\t", avg_latency);
				printf("--------------------------------------------------\n");

				for (int c = 0; c < core_list[indx_core]; ++c)
				{
					free(info_byte[c]);
					free(info_bits[c]);
					free(rmed_bits[c]);
					free(mapped_sig[c]);
					free(noise[c]);
					free(llr[c]);
					free(decoded_llr[c]);
					free(decbs_bits[c]);
					free_nr5g_ldpc_simd_t(ldpc_arg[c]);
				}
				pool_destroy(0);
			}
		fprintf(fp, "\n");
		fprintf(fq, "\n");
		printf("==================================================\n\n");
	}

	fclose(fp);
	fclose(fq);

	return 0;
}