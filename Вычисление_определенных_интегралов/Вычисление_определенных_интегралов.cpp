#include <iostream>
#include <math.h>
#include <iomanip>
using namespace std;
const double eps1 = 0.0001;
const double eps2 = 0.00001;

double myFunc(double x)
{
	return (3 * x / pow(1 + x * x * x, 0.5));
}

double function(double x, double y)
{
	return (4-x*x-y*y);
}

double MethodOfTrap(double a, double b)
{
	double sum = 0;
	double sum_prev = 0;
	double h = (b - a);
	int kol_vo = 0;

	do {
		sum_prev = sum;              // для нахождения кол-ва итераций
		sum = 0;					// зануление суммы
		for (int i = 1; i <= ((b - a) / h) - 1; i++)        //(b-a)/h = n 
		{
			sum += 2 * myFunc(a + h * i);
		}
		sum = (myFunc(a) + sum + myFunc(b)) * h / 2;		//Формула трапеций - для подсчета определенного интеграла
		h /= 2;
		kol_vo++;
	} while (abs(sum - sum_prev) > 3*eps1);							//критерий завершения вычисительного процесса
	cout << "myTrampezia: " << endl << kol_vo << " iter" << endl;
	return sum;
}

double MethodOfSimpson(double a, double b)
{
	double sum = 0; double sum_prev = 0;
	double h = (b - a) / 2;                         //количество интервалов четное 
	int kol_vo = 0;

	do {
		sum_prev = sum;
		sum = 0;
		for (int i = 1; i <= ((b - a) / h); i += 2)
		{
			sum += 4 * myFunc(a + h * i);
		}
		for (int i = 2; i <= ((b - a) / h) - 1; i += 2)
		{
			sum += 2 * myFunc(a + h * i);
		}
		sum = (myFunc(a) + sum + myFunc(b)) * h / 3;
		h /= 2;                                  // чётность 
		kol_vo++;
	} while (abs(sum - sum_prev) > 15*eps2);
	cout << "mySimpson: " << endl << kol_vo << " iter" << endl;
	return sum;
}

double MethodOfCubeSimpson(double a, double b, double c, double d)
{
	int m = 2; int n = m;
	int kol_vo = 0;
	double sum = 0; double sum_prev = 0;

	do {
		sum_prev = sum;
		sum = 0;

		double hx = (b - a) / (2*n);
		double hy = (d - c) / (2*m);

		double xi = a;
		double yi = c;

		double* Xi = new double[2 * n + 1];
		Xi[0] = xi;

		for (int i = 1; i <= 2 * n; i++)
			Xi[i] = Xi[i - 1] + hx;

		double* Yi = new double[2 * m + 1];
		Yi[0] = yi;

		for (int j = 1; j <= 2 * m; j++)
			Yi[j] = Yi[j - 1] + hy;

		for (int i = 0; i < n; i++)
		{
			for (int j = 0; j < m; j++)
			{
				sum += function(Xi[2 * i], Yi[2 * j]);
				sum += 4 * function(Xi[2 * i + 1], Yi[2 * j]);
				sum += function(Xi[2 * i + 2], Yi[2 * j]);
				sum += 4 * function(Xi[2 * i], Yi[2 * j + 1]);
				sum += 16 * function(Xi[2 * i + 1], Yi[2 * j + 1]);
				sum += 4 * function(Xi[2 * i + 2], Yi[2 * j + 1]);
				sum += function(Xi[2 * i], Yi[2 * j + 2]);
				sum += 4 * function(Xi[2 * i + 1], Yi[2 * j + 2]);
				sum += function(Xi[2 * i + 2], Yi[2 * j + 2]);
			}
		}
		sum *= (hx * hy / 9);
		n *= 2;
		m *= 2;
		kol_vo++;
	} while (abs(sum - sum_prev) > eps2);
	cout << "Cube Simpson: " << endl << kol_vo << " iter" << endl;
	return sum;
}

void main()
{
	double start = 0;
	double end = 1.075;

	double trap;
	double simp;
	double cube_simp;

	trap = MethodOfTrap(start, end);
	cout << trap << endl << endl;

	simp = MethodOfSimpson(start, end);
	cout << simp << endl << endl;

	double a = -1;
	double b = 1;
	double c = -1;
	double d = 1;
	cube_simp = MethodOfCubeSimpson(a, b, c, d);
	cout << cube_simp << endl;
}