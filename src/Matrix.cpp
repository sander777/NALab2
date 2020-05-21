#ifndef _MATRIX_CPP_
#define _MATRIX_CPP_

#include "Matrix.hpp"
#include <cmath>
#include <ctime>
#include <iostream>
#include <random>

double _abs(double a) {
    return (a < 0) ? a * -1: a;
}

Matrix::Matrix() { _matrix = new double *[1]; }

Matrix::Matrix(int n, int m) {
    this->_matrix = new double *[m];
    for (int i = 0; i < m; i++)
        this->_matrix[i] = new double[n]{};
    _n = n;
    _m = m;
}

Matrix::Matrix(const Matrix &copy) {
    _n = copy._n;
    _m = copy._m;
    _matrix = new double *[_m];
    for (int i = 0; i < _m; i++)
        this->_matrix[i] = new double[_n]{};

    for (int i = 0; i < _m; i++) {
        for (int j = 0; j < _n; j++) {
            _matrix[i][j] = copy._matrix[i][j];
        }
    }
}

Matrix::Matrix(double *arr, int n, int m) {
    this->_matrix = new double *[m];
    for (int i = 0; i < m; i++)
        this->_matrix[i] = new double[n]{};
    _n = n;
    _m = m;

    for (int i = 0; i < _m; i++) {
        for (int j = 0; j < _n; j++) {
            _matrix[i][j] = arr[i * _n + j];
        }
    }
}

Matrix::~Matrix() {
    for (int i = 0; i < _m; i++)
        delete this->_matrix[i];
    delete _matrix;
}

double Matrix::at(int m, int n) const { return _matrix[m][n]; }

const Matrix &Matrix::operator=(const Matrix &copy) {
    for (int i = 0; i < _m; i++)
        delete[] _matrix[i];
    delete[] _matrix;

    _n = copy._n;
    _m = copy._m;
    _matrix = new double *[_m];
    for (int i = 0; i < _m; i++)
        this->_matrix[i] = new double[_n]{};

    for (int i = 0; i < _m; i++) {
        for (int j = 0; j < _n; j++) {
            _matrix[i][j] = copy._matrix[i][j];
        }
    }
}

double *Matrix::operator[](int m) { return _matrix[m]; }

const Matrix &Matrix::operator*(const Matrix &mx) {
    if (_n != mx.getM())
        throw "Wrong size of multiplicator";

    Matrix *result = new Matrix(mx.getN(), _m);
    for (int i = 0; i < _m; i++) {
        for (int j = 0; j < mx.getN(); j++) {
            double sum = 0;
            for (int p = 0; p < _n; p++) {
                sum += _matrix[i][p] * mx.at(p, j);
            }
            (*result)[i][j] = sum;
        }
    }

    return *result;
}

const Matrix &Matrix::operator*=(const Matrix &mx) {
    if (_n != mx.getM() || _m != mx.getN())
        throw "Wrong size of multiplicator";

    Matrix *result = new Matrix(_n, _m);
    for (int i = 0; i < _m; i++) {
        for (int j = 0; j < _n; j++) {
            double sum = 0;
            for (int p = 0; p < _n; p++) {
                sum += _matrix[i][p] * mx.at(p, j);
            }
            (*result)[i][j] = sum;
        }
    }
    *this = *result;
    return *this;
}

void Matrix::print() const {
    for (int i = 0; i < _m; i++) {
        for (int j = 0; j < _n; j++) {
            std::cout.width(ALIGN_SIZE);
            std::cout << std::right << _matrix[i][j];
        }
        std::cout << '\n';
    }
}

int Matrix::getN() const { return _n; }
int Matrix::getM() const { return _m; }

Matrix &randomMatrix(int n, int m) {
    Matrix *result = new Matrix(n, m);
    srand(time(0));
    for (int i = 0; i < m; i++) {
        for (int j = 0; j < n; j++) {
            (*result)[i][j] = rand() % 21 - 10;
        }
    }
    return *result;
}

Matrix &transpose(const Matrix &mx) {
    Matrix *result = new Matrix(mx.getM(), mx.getN());
    for (int i = 0; i < mx.getM(); i++) {
        for (int j = 0; j < mx.getN(); j++) {
            (*result)[j][i] = mx.at(i, j);
        }
    }
    return *result;
}

Matrix &cholesky(const Matrix &mx) {
    Matrix *result = new Matrix(mx.getN(), mx.getM());
    (*result)[0][0] = sqrt(mx.at(0, 0));
    for (int j = 1; j < mx.getN(); j++) {
        (*result)[j][0] = mx.at(j, 0) / (*result)[0][0];
    }
    for (int i = 0; i < mx.getN(); i++) {
        double sum = 0;
        for (int p = 0; p < i; p++) {
            sum += pow(result->at(i, p), 2);
        }
        (*result)[i][i] = sqrt(mx.at(i, i) - sum);
        for (int j = i + 1; j < mx.getN() && i < mx.getN() - 1; j++) {
            double sum = 0;
            for (int p = 0; p < i; p++) {
                sum += result->at(i, p) * result->at(j, p);
            }
            (*result)[j][i] = (mx.at(j, i) - sum) / result->at(i, i);
        }
    }

    return *result;
}

double norm(const Matrix &mx) {
    double max = 0;
    for (int i = 0; i < mx.getM(); i++) {
        double temp = 0;
        for (int j = 0; j < mx.getM(); j++) {
            temp += _abs(mx.at(i, j));
        }
        if (temp > max) {
            max = temp;
        }
    }
    return max;
}

Matrix &inverse(const Matrix &mx) {
    int n = mx.getN();
    Matrix a(mx);
    Matrix &ia = *(new Matrix(n, n));
    for (int i = 0; i < n; i++) {
        ia[i][i] = 1;
    }

    for (int i = 0; i < n; i++) {
        double d = a[i][i];
        for (int j = 0; j < n; j++) {
            a[i][j] /= d;
            ia[i][j] /= d;
        }
        for (int j = 0; j < n; j++) {
            if (j == i)
                continue;
            double mult = a[j][i];
            for (int k = 0; k < n; k++) {
                a[j][k] -= a[i][k] * mult;
                ia[j][k] -= ia[i][k] * mult;
            }
        }
    }
    return ia;
}

double conditionNumber(const Matrix &mx) { return norm(mx) * norm(inverse(mx)); }

#endif