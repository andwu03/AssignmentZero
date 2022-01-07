/*
 * LeastSquares.cpp
 *
 *  Created on: Sep. 10, 2020
 *      Author: Andrew Wu
 */

#include "LeastSquares.h"


#include <iostream>
#include <cmath>

using namespace std;

void printMatrix(int r, int c, double matrix[]) {//static method to print a given matrix
	//uses row, column, and flattened matrix. Matrices are flattened to effectively pass a 2D array into a function

	for (int k = 0; k < r; k++) {//loops through the matrix
		for (int i = 0; i < c; i++) {
			cout << matrix[c * k + i] << " ";//c*k + i accesses an element from the flattened 1D array that
			//represents an element in the original 2D array
		}
		cout << "\n";
	}
}

double getDeterminant(double matrix[], int cSize) {//static method to find the determinant of a given matrix
	//uses flattened matrix and its current dimension
	double unflatten[cSize][cSize];//create matrix to unflatten the flattened matrix

	for (int k = 0; k < cSize; k++) {//for loop to unflatten matrix
		for (int i = 0; i < cSize; i++) {
			unflatten[k][i] = matrix[cSize * k + i];//stores element from given matrix into unflatten
		}
	}

	double det = 0;//variable to store determinant

	if (cSize == 1) {//if current dimension of matrix is 1, return det = element
		return matrix[0];
	} else {

		double temp[cSize - 1][cSize - 1];//temp array to store co-factors

		double multi = 1;//multiplication co-efficient

		for (int k = 0; k < cSize; k++) {//for loop to access each row of the given matrix

			int i = 0, j = 0, cRow = 0, cCol = k;//i and j store the row and col of the cofactor matrix.
			//cRow and cCol store the current row and col of the given matrix

			for (int r = 0; r < cSize; r++) {
				for (int c = 0; c < cSize; c++) {

					if (r != cRow && c != cCol) {//add the element from the given matrix only if it
						//does not belong to the given row and col of the given matrix
						temp[i][j++] = unflatten[r][c];//stores into temp

						if (j == cSize - 1) {//once the row is filled, move onto next one and reset col
							j = 0;
							i++;
						}
					}
				}
			}

			double flatten[(cSize - 1) * (cSize - 1)];//create array to store flattened cofactor matrix

			for (int j = 0; j < cSize - 1; j++) {
				for (int i = 0; i < cSize - 1; i++) {
					flatten[(cSize - 1) * j + i] = temp[j][i];
				}
			}

			det += multi * unflatten[0][k] * getDeterminant(flatten, cSize - 1);//calculates determinant recursively
			//using the first row of the given matrix

			multi = -multi;//changes the co-efficient as the sign needs to alternate
		}

		return det;//returns det value
	}

}

int main() {
	int row, col, input;//declare variables

	//Data Collection + Vandermonde Matrix

	cout << "Please input the degree of the polynomial you want to fit: ";
	cin >> input;
	col = input + 1;//number of colums of Vandermonde matrix is one more than the degree of interest

	cout << "Please input the number of ordered pairs: ";
	cin >> row;//number of ordered pairs determines the row

//	double xValues[9] = { 1.01, 2.2, 2.9, 4.03, 5.32, 6.22, 8.56, 9.09 , 10};
//	double yValues[9][1] = { 18.5, 76.2, 150.5, 365, 780, 1265, 3250, 7099 , 16000};

	double xValues[row];//1D matrices to store values
	double yValues[row][1];
	double Vmatrix[row][col];

	cout << "Please input the X values: \n";//receiving input from user

	for (int k = 0; k < row; k++) {
		cout << "Value " << k + 1 << " is: ";
		cin >> xValues[k];//storing value into matrix
	}

	cout << "Please input the Y values: \n";

	for (int k = 0; k < row; k++) {
		cout << "Value " << k + 1 << " is: ";
		cin >> yValues[k][0];
	}

	for (int k = 0; k < row; k++) {//for loop that loops through rows and cols of V matrix
		for (int i = 0; i < col; i++) {
			Vmatrix[k][i] = std::pow(xValues[k], i);//raises col to relevant power, in this case it is current col number
		}
	}

	double VmOneD[row * col];//flattens V matrix for use elsewhere in program

	for (int j = 0; j < row; j++) {
		for (int i = 0; i < col; i++) {
			VmOneD[col * j + i] = Vmatrix[j][i];
		}
	}

	//Transpose Matrix

	double Tmatrix[col][row];

	for (int k = 0; k < col; k++) {//for loop that switches, so this one loops through cols and rows to transpose the V matrix
		for (int i = 0; i < row; i++) {
			Tmatrix[k][i] = std::pow(xValues[i], k);//uses k, which represents row number to raise the power
		}
	}

	double TpOneD[row * col];//flattens T matrix for use elsewhere in program

	for (int j = 0; j < col; j++) {
		for (int i = 0; i < row; i++) {
			TpOneD[row * j + i] = Tmatrix[j][i];
		}
	}

	//Multiplied Matrix

	double VtVmatrix[col][col];//creates a new 2D array of required dimensions

	double value = 0.0;//value to store the product of the elements

	for (int j = 0; j < col; j++) {//for loop to put necessary elements into the VtV matrix
		for (int k = 0; k < col; k++) {
			for (int i = 0; i < row; i++) {
				value = value + Tmatrix[k][i] * Vmatrix[i][j];//to get the product, we move along the
				//given row of T matrix while we move down the column of V matrix. The column of V matrix does not change
				//until we have multiplied all of T matrix to it
			}

			VtVmatrix[k][j] = value;//stores the multiplied value in the relevant position in VtV matrix

			value = 0.0;//resets value
		}
	}

	double VtVOneD[col * col];//flattens VtV matrix for use elsewhere in the program

	for (int k = 0; k < col; k++) {
		for (int i = 0; i < col; i++) {
			VtVOneD[col * k + i] = VtVmatrix[k][i];
		}
	}

	if (getDeterminant(VtVOneD, col) == 0) {//if the det of VtV is 0, we cannot find the inverse
		cout << "Inverse does not exist. Determinant of VtV is 0.";
	} else {//if det is non-zero, we move on with the computation

		//Adjoint Matrix

		double ADJmatrix[col][col];//create a 2D array to store the adjoint matrix

		if (col == 1) {//if the col of VtV happens to be one
			ADJmatrix[0][0] = 1;//adjdoint is then 1
		} else {

			double sign = 1;//store the co-efficient to multiply by

			double temp[col - 1][col - 1];//temp array to store cofactor

			for (int i = 0; i < col; i++) {
				for (int j = 0; j < col; j++) {//find the cofactor

					int q = 0, k = 0, cRow = i, cCol = j;

					for (int r = 0; r < col; r++) {
						for (int c = 0; c < col; c++) {

							if (r != cRow && c != cCol) {

								temp[q][k++] = VtVmatrix[r][c];

								if (k == col - 1) {
									k = 0;
									q++;
								}
							}
						}
					}

					if ((i + j) % 2 == 1) {
						sign = -1;//sign of the adjoin matrix is -1 if row and col index sum is odd
					} else {
						sign = 1;
					}

					double flatten[(col - 1) * (col - 1)];//flatten the cofactor matrix

					for (int a = 0; a < col - 1; a++) {
						for (int b = 0; b < col - 1; b++) {
							flatten[(col - 1) * a + b] = temp[a][b];
						}
					}

					ADJmatrix[j][i] = (sign)
							* (getDeterminant(flatten, col - 1));//finds the necessary det value
					// and also transposes it, note that j and i are reversed to transpose
				}
			}

		}

		double ADJOneD[col * col];//flattens adjoint matrix for future use

		for (int k = 0; k < col; k++) {
			for (int i = 0; i < col; i++) {
				ADJOneD[col * k + i] = ADJmatrix[k][i];
			}
		}

		//Inverse Matrix

		double INVmatrix[col][col];//create new matrix to store the inverse matrix
		double detVtV = getDeterminant(VtVOneD, col);//calculates the determinant of VtV

		for (int k = 0; k < col; k++) {
			for (int i = 0; i < col; i++) {//inverse if found by mapping ADJ to INV but dividing each
				//element of VtV
				INVmatrix[k][i] = static_cast<double>(ADJmatrix[k][i] / detVtV);
			}
		}

		double INVOneD[col * col];//flattens inverse matrix for future use

		for (int k = 0; k < col; k++) {
			for (int i = 0; i < col; i++) {
				INVOneD[col * k + i] = INVmatrix[k][i];
			}
		}

		//Co-efficient Calculation

		double VtYmatrix[col][1];//create matrix to aid in calculation of fit co-efficients

		value = 0.0;//stores value used during multiplication of matrices

		for (int k = 0; k < col; k++) {//multiplies T matrix with Y values
			for (int i = 0; i < row; i++) {
				value = value + Tmatrix[k][i] * yValues[i][0];
			}

			VtYmatrix[k][0] = value;

			value = 0.0;
		}

		double INVVtY[col][1];//create matrix to store co-efficients

		value = 0.0;

		for (int k = 0; k < col; k++) {//multiplies INV matrix with VtY matrix
			for (int i = 0; i < col; i++) {
				value = value + INVmatrix[k][i] * VtYmatrix[i][0];
			}

			INVVtY[k][0] = value;

			value = 0.0;
		}

		//Output
		cout << "\nX Values are: \n";
		for (int k = 0; k < row; k++) {
			cout << xValues[k] << "\n";
		}

		cout << "\nY Values are: \n";
		for (int k = 0; k < row; k++) {
			cout << yValues[k][0] << "\n";
		}

		cout << "\nV Matrix is \n";
		printMatrix(row, col, VmOneD);

		cout << "\nVt matrix is \n";
		printMatrix(col, row, TpOneD);

		cout << "\nVtV matrix is \n";
		printMatrix(col, col, VtVOneD);

		cout << "\nADJ matrix is \n";
		printMatrix(col, col, ADJOneD);

		cout << "\nINVmatrix is \n";
				printMatrix(col, col, INVOneD);

		cout << "\nFinal Co-eff:\n";
		for (int k = 0; k < col; k++) {
			cout << INVVtY[k][0] << "\n";
		}
	}
}
