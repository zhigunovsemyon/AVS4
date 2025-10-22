#include <math.h>
#include <omp.h>
#include <stdint.h>
#include <stdio.h>

static double fifth(double x)
{
	return asin(1.0 / sqrt(1.1 + x));
}

static double
simpson_calc(double begin, double end, int num_steps, double (*func)(double))
{
	double sum = 0.0;
	double step = (end - begin) / (double)num_steps;

#pragma omp parallel for reduction(+ : sum)
	for (int i = 1; i <= num_steps; ++i) {
		double x2 = begin + i * step;
		double x1 = x2 - step;
		double x12 = x1 + step * 0.5;

		sum += func(x2) + func(x1) + func(x12);
	}

	return sum * step / 3.0;
}

static double
trapezoid_calc(double begin, double end, int num_steps, double (*func)(double))
{
	double sum = 0.0;
	double step = (end - begin) / (double)num_steps;

#pragma omp parallel for reduction(+ : sum)
	for (int i = 1; i <= num_steps; ++i) {
		double x1 = begin + i * step;
		double x2 = x1 - step;

		sum += (func(x2) + func(x1));
	}

	return sum * 0.5 * step;
}

static double
right_rect_calc(double begin, double end, int num_steps, double (*func)(double))
{
	double sum = 0.0;
	double step = (end - begin) / (double)num_steps;

#pragma omp parallel for reduction(+ : sum)
	for (int i = 1; i <= num_steps; ++i) {
		double x = begin + i * step;
		double h = func(x);
		sum += h * step;
	}

	return sum;
}

static double centre_rect_calc(double begin,
			       double end,
			       int num_steps,
			       double (*func)(double))
{
	double sum = 0.0;
	double step = (end - begin) / (double)num_steps;

#pragma omp parallel for reduction(+ : sum)
	for (int i = 0; i < num_steps; ++i) {
		double x = begin + i * step;
		double h = func(x + 0.5 * step);
		sum += h * step;
	}

	return sum;
}

static double
left_rect_calc(double begin, double end, int num_steps, double (*func)(double))
{
	double sum = 0.0;
	double step = (end - begin) / (double)num_steps;

#pragma omp parallel for reduction(+ : sum)
	for (int i = 0; i < num_steps; ++i) {
		double x = begin + i * step;
		double h = func(x);
		sum += h * step;
	}

	return sum;
}

static void calculate(double epsilon,
		      double begin,
		      double end,
		      double (*method)(double begin,
				       double end,
				       int num_steps,
				       double (*)(double)),
		      double (*math_func)(double))
{
	double res, prev = INFINITY;

	double start_time = omp_get_wtime();
	for (int num_steps = 1;; num_steps++, prev = res) {
		res = method(begin, end, num_steps, math_func);

		if (fabs(res - prev) < epsilon) {
			double dt = omp_get_wtime() - start_time;
			printf("Точность %.1le достигнута за %lf с. Число "
			       "шагов: %d\n",
			       epsilon, dt, num_steps);
			printf("Результат: %.10lf\n", res);
			return;
		}
	}
}

int main(void)
{
	double const BEGIN = 1.0, END = 3.0;

	puts("Метод левых прямоугольников:");
	calculate(0.00001, BEGIN, END, left_rect_calc, fifth);
	calculate(0.000001, BEGIN, END, left_rect_calc, fifth);
	calculate(0.0000001, BEGIN, END, left_rect_calc, fifth);
	putchar('\n');

	puts("Метод правых прямоугольников:");
	calculate(0.00001, BEGIN, END, right_rect_calc, fifth);
	calculate(0.000001, BEGIN, END, right_rect_calc, fifth);
	calculate(0.0000001, BEGIN, END, right_rect_calc, fifth);
	putchar('\n');

	for (double e = 0.00001; e > 0.000000000000001; e /= 10) {

		puts("Метод средних прямоугольников:");
		calculate(e, BEGIN, END, centre_rect_calc, fifth);
		putchar('\n');

		puts("Метод трапеций:");
		calculate(e, BEGIN, END, trapezoid_calc, fifth);
		putchar('\n');

		puts("Метод Симпсона:");
		calculate(e, BEGIN, END, simpson_calc, fifth);
		putchar('\n');
	}

	return 0;
}
