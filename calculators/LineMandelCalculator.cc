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
    int *pdata = data, r;
	float *xs = xvals, *ys = yvals, *xc = xcalc, xb, r2, i2, y;

	r = 0;
    for (int i = 0; i < height/2; i++) {
		y = y_start + i * dy;
		#pragma omp simd aligned(xs, ys, xc, pdata: 64)
		for (int j = 0; j < width; j++) {
			pdata[r + j] = 0;
			xb = x_start + (j * dx);
			xs[j] = xb;
			ys[j] = y;
			xc[j] = xb;
		}
<<<<<<< HEAD
		#pragma omp simd aligned(xs, ys, xc, pdata:64) simdlen(16)
		for (int k = 0, escaped=0; k < limit; ++k) {
=======
		int first = 0;
		bool f;
		for (int k = 0, escaped=0; k < limit && escaped < width; ++k) {
>>>>>>> e7e224b0078864a19ccc6ca585b2b0686e020431
			escaped = 0;
			f = true;
			#pragma omp simd aligned(xs, ys, xc, pdata:64) simdlen(16)
			for (int j = 0; j < width; j++) {

				r2 = xs[j] * xs[j];
				i2 = ys[j] * ys[j];

				if ((r2 + i2 < 4.0f)) {
					pdata[r + j]++;
					if (f && first <= j) {
						f = false;
						first = j;
					}
				} else {
					escaped++;
				}
				
				ys[j] = 2.0f * xs[j] * ys[j] + y;
				xs[j] = r2 - i2 + xc[j];
			}
			if (escaped >= width) break;
        }
<<<<<<< HEAD
		std::memcpy(&pdata[(height-1) * width - r], &pdata[r], width*sizeof(int));

=======
		std::memcpy(&pdata[(height-i-1) * width], &pdata[r], width*sizeof(int));
		r += width;
>>>>>>> e7e224b0078864a19ccc6ca585b2b0686e020431
    }
    return data;
}
