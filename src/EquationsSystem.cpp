#ifndef _EQUATIONS_SYSTEM_CPP_
#define _EQUATIONS_SYSTEM_CPP_
#include "EquationsSystem.hpp"
#include <iostream>
#include <math.h>

EquationsSystem::EquationsSystem(const Matrix &A, const Matrix &b) {
    _A = A;
    if (b.getM() != _A.getM() || b.getN() != 1)
        throw "Wrong b size";
    _b = b;
}

EquationsSystem::~EquationsSystem() {}

void EquationsSystem::print() {
    for (int i = 0; i < _A.getM(); i++) {
        for (int j = 0; j < _A.getN(); j++) {
            std::cout.width(ALIGN_SIZE);
            std::cout << std::right << _A[i][j];
        }
        std::cout.width(ALIGN_SIZE);
        std::cout << '|';
        std::cout.width(ALIGN_SIZE);
        std::cout << _b[i][0] << '\n';
    }
}

const Matrix &EquationsSystem::getA() const { return _A; }

const Matrix &EquationsSystem::getB() const { return _b; }

Matrix &thomas(const EquationsSystem &eq) {
    int n = eq.getA().getN();
    double *c = new double[n];
    double *a = new double[n - 1];
    double *b = new double[n - 1];
    double *f = new double[n];
    double *x = new double[n];

    for (int i = 0; i < n; i++) {
        c[i] = eq.getA().at(i, i);
    }

    for (int i = 0; i < n - 1; i++) {
        a[i] = eq.getA().at(i + 1, i);
    }

    for (int i = 0; i < n - 1; i++) {
        b[i] = eq.getA().at(i, i + 1);
    }

    for (int i = 0; i < n; i++) {
        f[i] = eq.getB().at(i, 0);
    }

    double m;
    for (int i = 1; i < n; i++) {
        m = a[i - 1] / c[i - 1];
        c[i] = c[i] - m * b[i - 1];
        f[i] = f[i] - m * f[i - 1];
    }

    x[n - 1] = f[n - 1] / c[n - 1];

    for (int i = n - 2; i >= 0; i--) {
        x[i] = (f[i] - b[i] * x[i + 1]) / c[i];
    }

    Matrix *result = new Matrix(1, n);
    for (int i = 0; i < n; i++) {
        (*result)[i][0] = x[i];
    }

    return *result;
}

Matrix &gauss(const EquationsSystem &eq) {
    double *x, max;
    Matrix a = eq.getA();
    Matrix y = eq.getB();
    int n = eq.getA().getN();
    int k, index;
    x = new double[n];
    k = 0;
    while (k < n) {
        max = abs(a[k][k]);
        index = k;
        for (int i = k + 1; i < n; i++) {
            if (abs(a[i][k]) > max) {
                max = abs(a[i][k]);
                index = i;
            }
        }
        if (max < eps) {
            throw "Solution impossible due to null column";
        }
        for (int j = 0; j < n; j++) {
            double temp = a[k][j];
            a[k][j] = a[index][j];
            a[index][j] = temp;
        }
        double temp = y[k][0];
        y[k][0] = y[index][0];
        y[index][0] = temp;
        for (int i = k; i < n; i++) {
            double temp = a[i][k];
            if (abs(temp) < eps)
                continue;
            for (int j = 0; j < n; j++)
                a[i][j] = a[i][j] / temp;
            y[i][0] = y[i][0] / temp;
            if (i == k)
                continue;
            for (int j = 0; j < n; j++)
                a[i][j] = a[i][j] - a[k][j];
            y[i][0] = y[i][0] - y[k][0];
        }
        k++;
    }
    for (k = n - 1; k >= 0; k--) {
        x[k] = y[k][0];
        for (int i = 0; i < k; i++)
            y[i][0] = y[i][0] - a[i][k] * x[k];
    }
    Matrix *result = new Matrix(x, 1, n);
    delete [] x;
    return *result;
}

Matrix& jacobi (const EquationsSystem& eq, bool display)
{
    int N = eq.getA().getN();
    Matrix A = eq.getA();
    Matrix F = eq.getB();
    double* X = new double[N];
	double* TempX = new double[N];
	double norm;
    double j = 1;
    if(display) std::cout << "\n_____________________________\n";
	do {
		for (int i = 0; i < N; i++) {
			TempX[i] = F[i][0];
			for (int g = 0; g < N; g++) {
				if (i != g)
					TempX[i] -= A[i][g] * X[g];
			}
			TempX[i] /= A[i][i];
		}
        norm = fabs(X[0] - TempX[0]);
		for (int h = 0; h < N; h++) {
			if (fabs(X[h] - TempX[h]) > norm)
				norm = fabs(X[h] - TempX[h]);
			X[h] = TempX[h];
		}
        if(display) {
            std::cout.width(4);
            std::cout << std::right << j++ << "    ";
            for(int i = 0; i < N; i++) {
                std::cout.width(10);
                std::cout << std::right << TempX[i] << '|';
            }
            std::cout << '\n';

        }
	} while (norm > eps);
    Matrix* result = new Matrix(X, 1, N);
	delete[] TempX;
    delete[] X;
    return *result;
}

#endif