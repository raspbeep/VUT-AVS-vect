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
	line_b = (int *)(aligned_alloc(64, width * sizeof(int)));
	xvals = (float *)(aligned_alloc(64, width * sizeof(float)));
	yvals = (float *)(aligned_alloc(64, width * sizeof(float)));
	xcalc = (float *)(aligned_alloc(64, width * sizeof(float)));
}

LineMandelCalculator::~LineMandelCalculator() {
	free(data);
	free(xvals);
	free(yvals);
	free(xcalc);
	free(line_b);
	data=NULL;
	xvals=NULL;
	yvals=NULL;
	xcalc=NULL;
	line_b=NULL;
}

int *LineMandelCalculator::calculateMandelbrot() {
    int *pdata = data, r, w = sizeof(float) * width;
	float *xs = xvals, *ys = yvals, *xc = xcalc, xb, r2, i2, y;

	r = 0;
    for (int i = 0; i < height / 2; i++) {
		y = y_start + i * dy;
		#pragma omp simd aligned(xs: 64)
		for (int j = 0; j < width; j++) {
			xs[j] = x_start + (j * dx);
		}
		std::memcpy(xc, xs, w);
		std::uninitialized_fill(ys, ys + width, y);
		std::uninitialized_fill(line_b, line_b + width, 0);

		int first = 0;
		bool f;
		for (int k = 0, escaped=0; k < limit && escaped < width; ++k) {
			escaped = 0;
			f = true;
			#pragma omp simd aligned(xs, ys, xc, pdata:64) simdlen(16)
			for (int j = 0; j < width; j++) {

				r2 = xs[j] * xs[j];
				i2 = ys[j] * ys[j];

				if ((r2 + i2 < 4.0f)) {
					line_b[j]++;
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
		std:memcpy(&pdata[i * width], line_b, w);
		std::memcpy(&pdata[(height-i-1) * width], &pdata[r], w);
		r += width;
    }
    return data;
}
