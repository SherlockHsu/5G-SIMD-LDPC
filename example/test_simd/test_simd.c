// version 3.1
#include "simd_ldpc.h"
#include "simd_bit.h"
#include "bit.h"
#include "crc.h"

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

#define CRC_24A 0x1864CFB
#define BLOCK_SIZE 10000

#define B_NUM 2
#define R_NUM 3
#define SIMD_MODE_NUM 3
#define EBN0_SIZE 21

int main()
{
	/* test parameters */
	int simd_list[SIMD_MODE_NUM] = {SIMD_MODE_SSE, SIMD_MODE_AVX2, SIMD_MODE_AVX512};
	int B_list[B_NUM] = {8448, 3840};
	int R_list[R_NUM] = {853, 768, 512};
	float EbN0_list[EBN0_SIZE] = {0.0, 0.5, 1.0, 1.5, 2.0, 2.5, 3.0, 3.5, 4.0, 4.5, 5.0, 5.5, 6.0, 6.5, 7.0, 7.5, 8.0, 8.5, 9.0, 9.5, 10.0};
	int indx_B, indx_R, indx_block, indx_ebn0, indx_simd;

	/* file parameters */
	FILE *fber[SIMD_MODE_NUM], *fbler[SIMD_MODE_NUM], *ftp[SIMD_MODE_NUM], *fla[SIMD_MODE_NUM], *fatp[SIMD_MODE_NUM], *fala[SIMD_MODE_NUM];
	
	fber[0] = fopen("sse_BER.txt", "a");
	fbler[0] = fopen("sse_BLER.txt", "a");
	ftp[0] = fopen("sse_tp.txt", "a");
	fla[0] = fopen("sse_latency.txt", "a");
	fatp[0] = fopen("sse_avg_tp.txt", "a");
	fala[0] = fopen("sse_avg_latency.txt", "a");
	fber[1] = fopen("avx2_BER.txt", "a");
	fbler[1] = fopen("avx2_BLER.txt", "a");
	ftp[1] = fopen("avx2_tp.txt", "a");
	fla[1] = fopen("avx2_latency.txt", "a");
	fatp[1] = fopen("avx2_avg_tp.txt", "a");
	fala[1] = fopen("avx2_avg_latency.txt", "a");
	fber[2] = fopen("avx512_BER.txt", "a");
	fbler[2] = fopen("avx512_BLER.txt", "a");
	ftp[2] = fopen("avx512_tp.txt", "a");
	fla[2] = fopen("avx512_latency.txt", "a");
	fatp[2] = fopen("avx512_avg_tp.txt", "a");
	fala[2] = fopen("avx512_avg_latency.txt", "a");
	

	for (indx_B = 0; indx_B < B_NUM; indx_B++)
		for (indx_R = 0; indx_R < R_NUM; indx_R++)
			for (indx_simd = 0; indx_simd < SIMD_MODE_NUM; indx_simd++)
			{
				printf("==================================================\n");
				printf("B:%d\tR:%d\t%s\n", B_list[indx_B], R_list[indx_R],
					   (simd_list[indx_simd] == SIMD_MODE_SSE ? "SIMD_MODE_SSE" : (simd_list[indx_simd] == SIMD_MODE_AVX2 ? "SIMD_MODE_AVX2" : "SIMD_MODE_AVX512")));
				printf("--------------------------------------------------\n");
				/* time parameters */
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

				/* mkl stream parameters */
				VSLStreamStatePtr stream_g;
				VSLStreamStatePtr stream_b;

				int32_t i, sum_err_bits, test_size;
				int32_t B, R, I_max, decoder_mode;
				float coef;
				float EbN0, sigma2, sigma;

				nr5g_ldpc_simd_t *ldpc_arg;
				int8_t *info_byte;
				int8_t *info_bits;
				int8_t *rmed_bits;
				float *mapped_sig;
				float *noise;
				float *llr;
				float *decoded_llr;
				int8_t *decbs_bits;

				int32_t err_bits[BLOCK_SIZE];
				int32_t err_bl;

				nr5g_crc_t crc_t;
				nr5g_crc_init(&crc_t, CRC_24A, 24);

				/* set parameters */
				B = B_list[indx_B];
				R = R_list[indx_R];
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

				test_size = EBN0_SIZE;

				/* initialize */
				ldpc_arg = (nr5g_ldpc_simd_t *)malloc(sizeof(nr5g_ldpc_simd_t));
				nr5g_ldpc_simd_init(ldpc_arg, B, R, simd_list[indx_simd]);
				ldpc_arg->crc_t = &crc_t;

				info_bits = (int8_t *)malloc(sizeof(int8_t) * ldpc_arg->B);
				info_byte = (int8_t *)malloc(sizeof(int8_t) * ldpc_arg->B / 8);
				rmed_bits = (int8_t *)malloc(sizeof(int8_t) * ldpc_arg->G);
				mapped_sig = (float *)malloc(sizeof(float) * ldpc_arg->G);
				noise = (float *)malloc(sizeof(float) * ldpc_arg->G);
				llr = (float *)malloc(sizeof(float) * ldpc_arg->G);
				decoded_llr = (float *)malloc(sizeof(float) * ldpc_arg->G);
				decbs_bits = (int8_t *)malloc(sizeof(int8_t) * ldpc_arg->B);

				vslNewStream(&stream_g, VSL_BRNG_MCG31, 0);
				vslNewStream(&stream_b, VSL_BRNG_MCG31, 1);

				avg_tp = 0.0;
				avg_latency = 0.0;

				/* test loop start */
				for (indx_ebn0 = 0; indx_ebn0 < test_size; indx_ebn0++)
				{
					EbN0 = EbN0_list[indx_ebn0];
					sigma2 = (float)(1 / (pow(10, (double)EbN0 / 10) * 2 * R / 1024));
					// sigma2 = (float)(1 / pow(10, (double)EbN0 / 10));
					sigma = (float)sqrt(sigma2);
					sum_err_bits = 0;
					err_bl = 0;
					encode_run_time = 0.0;
					decode_run_time = 0.0;

					for (indx_block = 0; indx_block < BLOCK_SIZE; indx_block++)
					{
						err_bits[indx_block] = 0;

						/* generate random tbs */
						viRngUniformBits(VSL_RNG_METHOD_UNIFORMBITS_STD, stream_b, (B - 1) / 32 + 1, (uint32_t *)info_byte);
						nr5g_crc_attach_byte(&crc_t, info_byte, B - 24);
						// fast_extend_avx512(info_byte, B / 8, info_bits);
						nr5g_bit_unpack_vector(info_byte, info_bits, B);

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
						for (i = 0; i < ldpc_arg->G; i++)
							mapped_sig[i] = (float)2 * rmed_bits[i] - 1;

						/* pass AWGN channel */
						vsRngGaussian(VSL_RNG_METHOD_GAUSSIANMV_BOXMULLER, stream_g, ldpc_arg->G, noise, 0.0, sigma);
						for (i = 0; i < ldpc_arg->G; i++)
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

						nr5g_ldpc_simd_decoder(ldpc_arg->rdmed_llr, ldpc_arg, I_max, coef, decoder_mode, EARLY_STOP_OFF, ldpc_arg->decoded_bits, decoded_llr);

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
						if (err_bits[indx_block])
							err_bl++;

						sum_err_bits += err_bits[indx_block];
					}

					/* print results */
					float ber = (float)sum_err_bits / (B * BLOCK_SIZE);
					float bler = (float)err_bl / BLOCK_SIZE;
					printf("Eb/N0:%.2f:\tBER:%.2e(%d/%d)\tBLER:%.2e(%d/%d)\n", EbN0_list[indx_ebn0], ber, sum_err_bits, B * BLOCK_SIZE, bler, err_bl, BLOCK_SIZE);
					printf("encode_Latency:%lfus\n", encode_run_time / BLOCK_SIZE * 1e6);
					printf("encode_Throughput:%.2lfMbps\n", (double)B * BLOCK_SIZE / encode_run_time / 1e6);
					printf("decode_Latency:%lfus\n", decode_run_time / BLOCK_SIZE * 1e6);
					printf("decode_Throughput:%.2lfMbps\n", (double)B * BLOCK_SIZE / decode_run_time / 1e6);
					printf("--------------------------------------------------\n");
					fprintf(fber[indx_simd], "%.2e\t", ber);
					fprintf(fbler[indx_simd], "%.2e\t", bler);
					fprintf(ftp[indx_simd], "%.2lf\t", (double)B * BLOCK_SIZE / decode_run_time / 1e6);
					fprintf(fla[indx_simd], "%.2lf\t", (double)decode_run_time * 1e6 / BLOCK_SIZE);

					avg_tp += (double)B * BLOCK_SIZE / decode_run_time / 1e6;
					avg_latency += decode_run_time * 1e6 / BLOCK_SIZE;
				}

				avg_tp /= test_size;
				avg_latency /= test_size;
				fprintf(fatp[indx_simd], "%.2lf\t", avg_tp);
				fprintf(fala[indx_simd], "%.2lf\t", avg_latency);

				printf("B:%d\tR:%d\t%s\n", B_list[indx_B], R_list[indx_R], (indx_simd == 1 ? "SIMD_MODE_SSE" : (indx_simd == 2 ? "SIMD_MODE_AVX2" : "SIMD_MODE_AVX512")));
				printf("Average Throughput:\t%.2lfMbps\n", avg_tp);
				printf("Average Latency:\t%.2lfus\n", avg_latency);
				printf("==================================================\n\n");

				fprintf(fber[indx_simd], "\n");
				fprintf(fbler[indx_simd], "\n");
				fprintf(ftp[indx_simd], "\n");
				fprintf(fla[indx_simd], "\n");

				free(info_byte);
				free(info_bits);
				free(rmed_bits);
				free(mapped_sig);
				free(noise);
				free(llr);
				free(decoded_llr);
				free(decbs_bits);
				free_nr5g_ldpc_simd_t(ldpc_arg);
			}
	for (indx_simd = 0; indx_simd < SIMD_MODE_NUM; ++indx_simd)
	{
		fprintf(fber[indx_simd], "\n");
		fprintf(fbler[indx_simd], "\n");
		fprintf(ftp[indx_simd], "\n");
		fprintf(fla[indx_simd], "\n");
		fprintf(fatp[indx_simd], "\n");
		fprintf(fala[indx_simd], "\n");
		fclose(fber[indx_simd]);
		fclose(fbler[indx_simd]);
		fclose(ftp[indx_simd]);
		fclose(fla[indx_simd]);
	}
	return 0;
}