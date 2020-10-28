#include <iostream>
#include <fstream>
#include <stdio.h>
#include <iomanip>
#include <cmath>
using namespace std;

const int var = 24;
//24 f(x) = 2e^x-2x+3
//25 f(x) = 3^x + 2-x
const int N = 6;

const int givenN = 5;
const double a = 1.0;
const double b = 2.0;
const double h = (b - a) / givenN;

const char* fileOutput = "testout.txt"; //"output2.txt";//
ofstream output;

double getFunctionResult(double x)
{
	/*if (var == 24)
		return 2 * exp(x) - 2 * x + 3;
	else
		return pow(3, x) + 2 - x;*/
	return pow(3, x - 1) + 4 - x;
}

void printHead(int method)
{
	switch (method)
	{
	case 1:
		output << setw(15) << "x" << setw(15) << "f(x)" << setw(15) << "Pn(x)" << setw(20) << "Delta" << setw(20) << "Оценка \n";
		break;
	default:
		output << "Неизвестный метод";
		break;
	}
}

void printSerpDiffTable(double** table)
{
	output << " Таблица разделенных разностей \n";
	for (int i = 0; i < N; i++)
	{
		for (int j = 0; j <= N; j++)
		{
			//if (table[i][j] > 0)
			if (table[i][j] == 0)
				output << setw(15) << " ";
			//output << setw(15) << setprecision(7) << fixed << table[i][j];
			else
				output << setw(15) << setprecision(7) << fixed << table[i][j];
		}
		output << "\n";
	}
}

//Таблица разделенных разностей
double** createSerpDiffTable(bool reverse)
{
	double** table = new double* [N];
	for (int i = 0; i < N; i++)
		table[i] = new double[N + 1];

	for (int i = 0; i < N; i++)
		for (int j = 0; j <= N; j++)
			table[i][j] = 0;

	for (int i = 0; i < N; i++)
	{
		double x = 1 + i * h;
		if (!reverse)
		{
			table[i][0] = x;
			table[i][1] = getFunctionResult(x);
		}
		else
		{
			table[i][0] = getFunctionResult(x); //- 4.2320508;
			table[i][1] = x;
		}
	}

	for (int j = 2; j <= N; j++)
		for (int k = 0; k <= N - j; k++)
			table[k][j] = (table[k + 1][j - 1] - table[k][j - 1]) / (table[k + j - 1][0] - table[k][0]);

	return table;
}

//интерполяция с помощью формулы Ньютона;
void newtonInterpolation()
{
	output << " Интерполяционная формула Ньютона \n";
	double** SepDiffTable = createSerpDiffTable(false);
	printSerpDiffTable(SepDiffTable);

	//pow(3, x - 1) + 4 - x;
	double M6 = 3 * pow(log(3), 6);
	double x = 0;
	double fx = 0;
	double pnx = 0;
	double delta = 0;
	double error = 0;
	output << "\n M6 = " << M6 << "\n\n";

	printHead(1);
	int factorial = 1;
	for (int j = 2; j < N + 1; j++)
		factorial *= j;

	for (int i = 0; i < N - 1; i++)
	{
		x = a + (i + 0.5) * h;
		//x = a + i * h;
		double omega = 1;
		pnx = SepDiffTable[0][1];
		for (int k = 1; k < N; k++)
		{
			omega *= x - SepDiffTable[k - 1][0];
			pnx += omega * SepDiffTable[0][k + 1];
		}
		fx = getFunctionResult(x);
		error = 0;
		omega *= (x - SepDiffTable[5][0]);

		error = M6 * abs(omega) / factorial;
		output << setw(15) << setprecision(7) << fixed << x;
		output << setw(15) << setprecision(7) << fixed << fx;
		output << setw(15) << setprecision(7) << fixed << pnx;
		output << setw(20) << setprecision(10) << fixed << scientific << abs(pnx - fx);
		output << setw(20) << setprecision(10) << error << "\n";
	}

	for (int i = 0; i < N; i++)
		delete[] SepDiffTable[i];
	delete[] SepDiffTable;
}

void splainInterpolation()
{
	output << " Интерполяционная формула Ньютона \n";
	double M5 = 3 * pow(log(3), 5);
	double M4 = 3 * pow(log(3), 4);


}
double** createDiffTable()
{
	double** table = new double* [N];
	for (int i = 0; i < N; i++)
		table[i] = new double[2];

	for (int i = 0; i < N; i++)
	{
		double x = 1 + i * h;
		table[i][0] = x;
		table[i][1] = getFunctionResult(x);
	}

	return table;
}
void reverseNewton()
{
	output << "\n Решение уравнения методом обратной интерполяции\n";
	double c = 4.2320508;

	double** SepDiffTable = createSerpDiffTable(true);

	//printSerpDiffTable(SepDiffTable);
	double solution;

	double omega = 1;
	solution = SepDiffTable[0][1];
	for (int k = 1; k < N; k++)
	{
		omega *= c - SepDiffTable[k - 1][0]; //эт первый столбец с уже вычетом из C
		solution += omega * SepDiffTable[0][k + 1];
	}

	int j = 0;
	for (j = 0; SepDiffTable[j][0] < 0; j++);

	//полукостыль для вывода
	for (int i = 0; i < N; i++)
	{
		SepDiffTable[i][0] -= c;
	}
	printSerpDiffTable(SepDiffTable);

	output << "\n c = " << c;
	output << "\n Корень = " << solution << " j = " << j;
	output << "\n Неявязка = Abs(f(x)-c) = " << abs(getFunctionResult(solution) - c);
	for (int i = 0; i < N; i++)
		delete[] SepDiffTable[i];
	delete[] SepDiffTable;
}

int main()
{
	output.open(fileOutput);

	newtonInterpolation();
	//splainInterpolation();
	reverseNewton();
	output.close();
}

