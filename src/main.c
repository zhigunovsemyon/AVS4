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
	for(int i = 0; (double)i < END; ++i){
		printf("f(%d) = %lf\n", i, func(i));
	}
	return 0;
}
