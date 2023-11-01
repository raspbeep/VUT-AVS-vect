/**
 * @file LineMandelCalculator.cc
 * @author FULL NAME <xlogin00@stud.fit.vutbr.cz>
 * @brief Implementation of Mandelbrot calculator that uses SIMD paralelization over lines
 * @date DATE
 */
#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <cstring>

#include <stdlib.h>

#include "LineMandelCalculator.h"


LineMandelCalculator::LineMandelCalculator (unsigned matrixBaseSize, unsigned limit) :
	BaseMandelCalculator(matrixBaseSize, limit, "LineMandelCalculator") {
	data = (int *)(aligned_alloc(64, height * width * sizeof(int)));
	xvals = (float *)(aligned_alloc(64, width * sizeof(float)));
	yvals = (float *)(aligned_alloc(64, width * sizeof(float)));
	xcalc = (float *)(aligned_alloc(64, width * sizeof(float)));
}

LineMandelCalculator::~LineMandelCalculator() {
	free(data);
	free(xvals);
	free(yvals);
	free(xcalc);
	data=NULL;
	xvals=NULL;
	yvals=NULL;
	xcalc=NULL;
}

int *LineMandelCalculator::calculateMandelbrot() {
    int *pdata = data;
	float *xs = xvals, *ys = yvals, *xc = xcalc, xb, r2, i2;

    for (int i = 0; i < height/2; i++) {
		const float y = y_start + i * dy;
		const int r = i*width;

		#pragma omp simd aligned(xs, ys, xc, pdata: 64) simdlen(16)
		for (int j = 0; j < width; j++) {
			pdata[r + j] = 0;
			xb = x_start + j * dx;
			xs[j] = xb;
			ys[j] = y;
			xc[j] = xb;
		}
		
		for (int k = 0, escaped=0; k < limit && escaped < width; ++k) {
			escaped = 0;
			#pragma omp simd aligned(xs, ys, xc, pdata:64) reduction(+:escaped) simdlen(16)
			for (int j = 0; j < width; j++) {
				r2 = xs[j] * xs[j];
				i2 = ys[j] * ys[j];

				(r2 + i2 < 4.0f) ? 
					pdata[r + j]++ : 
					escaped++;
				
				ys[j] = 2.0f * xs[j] * ys[j] + y;
				xs[j] = r2 - i2 + xc[j];
			}
        }
		std::memcpy(&pdata[(height-i-1) * width], &pdata[r], width*sizeof(int));
    }
    return data;
}
