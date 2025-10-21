#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>
#include <stdlib.h>

double fifth(double x)
{
	return asin(1.0 / sqrt(1.1 + x));
}

double
left_rect_calc(double begin, double end, int num_steps, double (*func)(double))
{
	double sum = 0.0;
	double step = (end - begin) / (double)num_steps;

#pragma omp parallel for reduction(+ : sum)
	for (int i = 0; i <= num_steps; ++i) {
		double x = begin + i * step;
		double h = func(x);
		sum += h * step;
	}

	return sum;
}

void left_rect(double epsilon, double begin, double end, double (*func)(double))
{
	double prev = INFINITY;

	double start_time = omp_get_wtime();
	for (int num_steps = 1;; num_steps++) {
		double res = left_rect_calc(begin, end, num_steps, func);

		if (fabs(res - prev) < epsilon) {
			double dt = omp_get_wtime() - start_time;
			printf("Методом левых прямоугольников точность %.1le "
			       "достигнута за %lf с. Число шагов: %d\n",
			       epsilon, dt, num_steps);
			printf("Результат: %.10lf\n", res);
			return;
		}
		prev = res;
	}
}

int main(void)
{
	double const BEGIN = 1.0, END = 3.0;

	left_rect(0.00001, BEGIN, END, fifth);
	left_rect(0.000001, BEGIN, END, fifth);
	left_rect(0.0000001, BEGIN, END, fifth);
	return 0;
}
