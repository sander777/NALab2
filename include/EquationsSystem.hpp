#ifndef _EQUATIONS_SYSTEM_HPP_
#define _EQUATIONS_SYSTEM_HPP_

#include "Matrix.hpp"
#include "EquationsSystem.hpp"

#define eps 0.00001

class EquationsSystem {
  private:
    Matrix _A;
    Matrix _b;
  public:
    EquationsSystem() = delete;
    EquationsSystem(const Matrix&, const Matrix&);
    ~EquationsSystem();

    void print();

    const Matrix& getA() const;
    const Matrix& getB() const;
};

Matrix& thomas(const EquationsSystem&);
Matrix& gauss(const EquationsSystem&);
Matrix& jacobi(const EquationsSystem&, bool display);
#endif