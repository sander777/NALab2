#ifndef _MATRIX_HPP_
#define _MATRIX_HPP_

#define ALIGN_SIZE 10

class Matrix {
  private:
    int _n, _m;
    double **_matrix;

  public:
    Matrix();
    Matrix(int, int);
    Matrix(const Matrix &);
    Matrix(const double[], int, int);
    double *operator[](int);
    const Matrix &operator*(const Matrix &);
    const Matrix &operator*=(const Matrix &);
    const Matrix &operator=(const Matrix &);
    void print() const;
    double at(int, int) const;

    int getN() const;
    int getM() const;
    ~Matrix();
};

Matrix &randomMatrix(int, int);
Matrix &transpose(const Matrix &);
Matrix &cholesky(const Matrix &);

#endif