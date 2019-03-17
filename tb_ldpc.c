// version 2.2
#include "inc/simd_ldpc.h"

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

#define BLOCK_SIZE 10000
#define EBN0_SIZE 6

int main()
{
	int B_list[2] = {8448, 3840};
	int R_list[3] = {853, 768, 512};
	int j, k;
	for (int j = 0; j < 1; j++)
		for (int k = 0; k < 2; k++)
		{

			double encode_run_time;
			double decode_run_time;
#if defined(_MSC_VER)
			LARGE_INTEGER num;
			long long start, end, freq;
#else
			struct timeval start, end;
			long timeuse;
#endif
			double avg_tp, avg_latency;

			int32_t i, indx_block, indx_ebn0, sum_err_bits, test_size;
			int32_t B, R, I_max, decoder_mode, G;
			float coef;
			nr5g_ldpc_simd_t *ldpc_arg;
			VSLStreamStatePtr stream_g;
			VSLStreamStatePtr stream_b;
			int8_t *info_bits;
			int32_t *info_bits_32;
			int8_t *rmed_bits;
			float *mapped_sig;
			float *noise;
			float *llr;
			float *decoded_llr;
			int8_t *decbs_bits;
			float EbN0, sigma2, sigma;
			int32_t err_bits[BLOCK_SIZE];
			FILE *fp;

			/* set parameters */
			// B = 8448;
			// R = 853;//512,768,853
			B = B_list[j];
			R = R_list[k]; //512,768,853
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

			/* initialize */
			ldpc_arg = (nr5g_ldpc_simd_t *)malloc(sizeof(nr5g_ldpc_simd_t));
			nr5g_ldpc_simd_init(ldpc_arg, B, R);

			G = ldpc_arg->G;
			info_bits_32 = (int32_t *)malloc(sizeof(int32_t) * ldpc_arg->B);
			info_bits = (int8_t *)malloc(sizeof(int8_t) * ldpc_arg->B);
			rmed_bits = (int8_t *)malloc(sizeof(int8_t) * ldpc_arg->G);
			mapped_sig = (float *)malloc(sizeof(float) * ldpc_arg->G);
			noise = (float *)malloc(sizeof(float) * ldpc_arg->G);
			llr = (float *)malloc(sizeof(float) * ldpc_arg->G);
			decoded_llr = (float *)malloc(sizeof(float) * ldpc_arg->G);
			decbs_bits = (int8_t *)malloc(sizeof(int8_t) * ldpc_arg->B);

			vslNewStream(&stream_g, VSL_BRNG_MCG31, 0);
			vslNewStream(&stream_b, VSL_BRNG_MCG31, 1);

			fp = fopen("BER.txt", "a");

			avg_tp = 0.0;
			avg_latency = 0.0;

			/* test loop start */
			for (indx_ebn0 = 0; indx_ebn0 < test_size; indx_ebn0++)
			{
				EbN0 = EbN0_list[indx_ebn0];
				sigma2 = (float)(1 / (pow(10, (double)EbN0 / 10) * 2 * R / 1024));
				sigma = (float)sqrt(sigma2);
				sum_err_bits = 0;
				encode_run_time = 0.0;
				decode_run_time = 0.0;

				for (indx_block = 0; indx_block < BLOCK_SIZE; indx_block++)
				{
					err_bits[indx_block] = 0;

					/* generate random tbs */
					viRngBernoulli(VSL_RNG_METHOD_BERNOULLI_ICDF, stream_b, B, info_bits_32, 0.5);
					for (i = 0; i < B; i++)
						info_bits[i] = (int8_t)info_bits_32[i];

					/* cbs */
					nr5g_ldpc_simd_cbs(info_bits, ldpc_arg, ldpc_arg->cbs_bits);

					/* encode */
#if defined(_MSC_VER)
					QueryPerformanceFrequency(&num);
					freq = num.QuadPart;
					QueryPerformanceCounter(&num);
					start = num.QuadPart;
#else
					gettimeofday(&start, NULL);
#endif
					nr5g_ldpc_simd_encoder(ldpc_arg->cbs_bits, ldpc_arg, ldpc_arg->coded_bits);
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
					nr5g_ldpc_simd_rate_matching(ldpc_arg->coded_bits, ldpc_arg, rmed_bits);

					/* BPSK map */
					for (i = 0; i < G; i++)
						mapped_sig[i] = (float)2 * rmed_bits[i] - 1;

					/* pass AWGN channel */
					vsRngGaussian(VSL_RNG_METHOD_GAUSSIANMV_BOXMULLER, stream_g, G, noise, 0.0, sigma);
					for (i = 0; i < G; i++)
						llr[i] = 2 * (mapped_sig[i] + noise[i]) / sigma2;

					/* rate dematching */
					nr5g_ldpc_simd_rate_dematching(llr, ldpc_arg, ldpc_arg->rdmed_llr);

					/* decode */
#if defined(_MSC_VER)
					QueryPerformanceFrequency(&num);
					freq = num.QuadPart;
					QueryPerformanceCounter(&num);
					start = num.QuadPart;
#else
					gettimeofday(&start, NULL);
#endif
					nr5g_ldpc_simd_decoder_param_init(ldpc_arg);
					nr5g_ldpc_simd_decoder(ldpc_arg->rdmed_llr, ldpc_arg, I_max, coef, decoder_mode, ldpc_arg->decoded_bits, decoded_llr);
#if defined(_MSC_VER)
					QueryPerformanceCounter(&num);
					end = num.QuadPart;
					decode_run_time += (double)(end - start) / freq;
#else
					gettimeofday(&end, NULL);
					timeuse = 1000000 * (end.tv_sec - start.tv_sec) + end.tv_usec - start.tv_usec;
					decode_run_time += (double)timeuse / 1000000.0;
#endif

					/* decbs */
					nr5g_ldpc_simd_decbs(ldpc_arg->decoded_bits, ldpc_arg, decbs_bits);

					/* statistics */
					for (i = 0; i < ldpc_arg->B; i++)
						err_bits[indx_block] += (info_bits[i] == decbs_bits[i] ? 0 : 1);

					sum_err_bits += err_bits[indx_block];
				}

				/* print results */
				// printf("Eb/N0:%.2f:\tBER:\t%.2e(%d/%d)\n", EbN0_list[indx_ebn0], (float)sum_err_bits / B / BLOCK_SIZE, sum_err_bits, B * BLOCK_SIZE);
				// printf("encode_Latency:%lfus\n", encode_run_time * 1e6 / BLOCK_SIZE);
				// printf("encode_Throughput:%.2lfMbps\n", (double)B * BLOCK_SIZE / encode_run_time / 1e6);
				// printf("decode_Latency:%lfus\n", decode_run_time * 1e6 / BLOCK_SIZE);
				// printf("decode_Throughput:%.2lfMbps\n", (double)B * BLOCK_SIZE / decode_run_time / 1e6);
				// fprintf(fp, "%.2e\t", (float)sum_err_bits / B / BLOCK_SIZE);

				avg_tp += (double)B * BLOCK_SIZE / decode_run_time / 1e6;
				avg_latency += decode_run_time * 1e6 / BLOCK_SIZE;
			}
			fprintf(fp, "\n");
			fclose(fp);

			avg_tp /= test_size;
			avg_latency /= test_size;
			printf("\n");
			printf("B = %d, R = %d\n", B, R);
			printf("Average Throughput:\t%.2lfMbps\n", avg_tp);
			printf("Average Latency:\t%.2lfus\n", avg_latency);

			free(info_bits_32);
			free(info_bits);
			free(rmed_bits);
			free(mapped_sig);
			free(noise);
			free(llr);
			free(decoded_llr);
			free(decbs_bits);
			free_nr5g_ldpc_simd_t(ldpc_arg);
		}

	return 0;
}