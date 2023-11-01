/**
 * @file BatchMandelCalculator.cc
 * @author FULL NAME <xlogin00@stud.fit.vutbr.cz>
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over small batches
 * @date DATE
 */

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>

#include <stdlib.h>
#include <stdexcept>

#include "BatchMandelCalculator.h"

BatchMandelCalculator::BatchMandelCalculator (unsigned matrixBaseSize, unsigned limit) :
	BaseMandelCalculator(matrixBaseSize, limit, "BatchMandelCalculator") {
	data = (int *)(aligned_alloc(64, height * width * sizeof(int)));
	xvals = (float *)(aligned_alloc(64, tile_size * tile_size * sizeof(float)));
	yvals = (float *)(aligned_alloc(64, tile_size * tile_size * sizeof(float)));
	xcalc = (float *)(aligned_alloc(64, tile_size * tile_size * sizeof(float)));
	ycalc = (float *)(aligned_alloc(64, tile_size * tile_size * sizeof(float)));
}

BatchMandelCalculator::~BatchMandelCalculator() {
	free(data);
	free(xvals);
	free(yvals);
	free(xcalc);
	free(ycalc);
	data=NULL;
	xvals=NULL; 
	yvals=NULL;
	xcalc=NULL;
	ycalc=NULL;
}

int * BatchMandelCalculator::calculateMandelbrot () {
	int *pdata = data, start, end, global_start, row_count, v, j_w;

	float *xs = xvals, *ys = yvals, *xc = xcalc, *yc = ycalc, xb, yb, r2, i2;
	const int n_tiles = width / tile_size, tile_size_s = tile_size*tile_size;
	
	
    for (int i = 0; i < (height / (2 * tile_size)); i++) {
		const int r = i * width;
		const int i_t_s = i * tile_size;
		const int i_t_s_w = i_t_s * width;

		for (int t = 0; t < n_tiles; t++) {
			start = t * tile_size;
			end = start + tile_size;
			const int i_t_s_w_t = i_t_s_w + start;

			#pragma omp simd aligned(xs, ys, xc, yc, pdata: 64) simdlen(16)
			for (int j = 0; j < tile_size; j++) {
				yb = y_start + ((i_t_s + j) * dy);
				row_count = j * tile_size;
				j_w = (j * width) + i_t_s_w_t;
				for (int k = 0; k < tile_size; k++) {
					pdata[j_w + k] = 0;
					xb = x_start + ((start + k) * dx);
					v = row_count + k; 
					xs[v] = xb;
					ys[v] = yb;
					xc[v] = xb;
					yc[v] = yb;
				}
			}

			for (int q = 0, escaped=0; (q < limit) &&
			 (escaped < (tile_size_s)); ++q) {
				escaped = 0;
				for (int j = 0; j < tile_size; j++) {
					row_count = j * tile_size;
					j_w = (j * width) + i_t_s_w_t;
					#pragma omp simd aligned(xs, ys, xc, yc, pdata:64) reduction(+:escaped) simdlen(16)
					for (int k = 0; k < tile_size; k++) {
						v = row_count + k;
						r2 = xs[v] * xs[v];
						i2 = ys[v] * ys[v];

						(r2 + i2 < 4.0f) ? 
							pdata[j_w + k]++ : 
							escaped++;
						
						ys[v] = 2.0f * xs[v] * ys[v] + yc[v];
						xs[v] = r2 - i2 + xc[v];
					}
				}
        	}
		}
		for (int p = 0; p < tile_size; p++) {
			std::memcpy(&pdata[(((height/tile_size) - i) * tile_size * width) - ((p+1)*width)],
						   &pdata[i*tile_size*width + p * width],
						    width*sizeof(int));
		}
    }
	
    return data;
}