#ifndef _EQUATIONS_SYSTEM_CPP_
#define _EQUATIONS_SYSTEM_CPP_
#include "EquationsSystem.hpp"
#include <iostream>

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

    for(int i = 0; i < n; i++) {
        c[i] = eq.getA().at(i, i);
    }

    for(int i = 0; i < n - 1; i++) {
        a[i] = eq.getA().at(i + 1, i);
    }

    for(int i = 0; i < n - 1; i++) {
        b[i] = eq.getA().at(i, i + 1);
    }

    for(int i = 0; i < n; i++) {
        f[i] = eq.getB().at(i, 0);
    }

    double m;
	for (int i = 1; i < n; i++)
	{
		m = a[i - 1]/c[i-1];
		c[i] = c[i] - m*b[i-1];
		f[i] = f[i] - m*f[i-1];
	}

	x[n-1] = f[n-1]/c[n-1];

	for (int i = n - 2; i >= 0; i--)
    {
		x[i]=(f[i]-b[i]*x[i+1])/c[i];
    }
    

    Matrix *result = new Matrix(1, n);
    for(int i = 0; i < n; i++) {
        (*result)[i][0] = x[i];
    }

    return *result;
}

#endif