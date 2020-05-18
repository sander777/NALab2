#include <iostream>
#include <cstdlib>
#include "Matrix.hpp"
#include "EquationsSystem.hpp"

int main(int argc, char *argv[]) {
	system("cls");
	double choleskyTestArray[] = {
		4, 3, 2, 1,
		3, 6, 4, 2,
		2, 4, 6, 3,
		1, 2, 3, 4,
	};
	Matrix choleskyTestMatrix(choleskyTestArray, 4, 4);
	std::cout << "A:\n";
	choleskyTestMatrix.print();
	std::cout << "A = L * transpose(L)\nL:\n";
	Matrix a = cholesky(choleskyTestMatrix);
	a.print();
	std::cout << '\n';
	double det = 1;
	for (int i = 0; i < a.getN(); i++) {
		det *= a[i][i];
	}
	std::cout << "det(A) = det(L)^2 = " << det * det << '\n';
	system("pause");


	system("cls");
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
	std::cout << "A|b:\n";
	thomasTestEquation.print();
	std::cout << "x:\n";
	thomas(thomasTestEquation).print();
}