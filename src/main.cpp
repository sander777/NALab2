#include <iostream>
#include "Matrix.hpp"
#include "EquationsSystem.hpp"

int main(int argc, char *argv[]) {
	double choleskyTestArray[] = {
		4, 3, 2, 1,
		3, 6, 4, 2,
		2, 4, 6, 3,
		1, 2, 3, 4,
	};
	Matrix choleskyTestMatrix(choleskyTestArray, 4, 4);
	Matrix a = cholesky(choleskyTestMatrix);
	a.print();
	std::cout << '\n';
	double thomasTestArrayA[] = {
		3, 2, 0, 0,
		1, 3, 2, 0,
		0, 1, 3, 2,
		0, 0, 1, 3,
	};

	double thomasTestArrayB[] = { 1, 2, 3, 4 };
	Matrix thomasTestMatrixA = Matrix(thomasTestArrayA, 4, 4);
	Matrix thomasTestMatrixB = Matrix(thomasTestArrayB, 1, 4);
	EquationsSystem thomasTestEquation(thomasTestMatrixA, thomasTestMatrixB);

	thomasTestEquation.print();
	thomas(thomasTestEquation).print();
}