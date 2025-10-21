#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>

#define BEGIN 1.0
#define END 3.0

double func(double x)
{
	return asin(1.0 / sqrt(1.1 + x));
}

int main(void)
{
	int num_steps = 10;
	double step = (END - BEGIN)/ (double)num_steps;

	// #pragma omp parallel for
	for(int i = 0; i <= num_steps; ++i){
		double x = BEGIN + i * step;
		printf("f(%lf) = %lf\n", x, func(x));
	}
	return 0;
}
