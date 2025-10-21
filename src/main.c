#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>

double fifth(double x)
{
	return asin(1.0 / sqrt(1.1 + x));
}

double left_rect(double begin, double end, int num_steps, double(*func)(double))
{
	double sum = 0.0;
	double step = (end - begin)/ (double)num_steps;

	#pragma omp parallel for reduction(+:sum)
	for(int i = 0; i <= num_steps; ++i){
		double x = begin + i * step;
		double h = func(x);
		sum += h * step;
	}

	return sum;
}

#define BEGIN 1.0
#define END 3.0

int main(void)
{
	printf("%lf\n", left_rect(BEGIN, END, 1488, fifth));
	return 0;
}
