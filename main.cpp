#include "decomp.hpp"

using namespace std;

#define INDX(i, j)  ((i) * mdim + (j))

void matrixInit(double *matrix, double p)
{
	matrix[0] = p + 27;		matrix[1] = -6;		matrix[2] = -1;		 matrix[3] = -6;	matrix[4] = -3;		 matrix[5] = -4;	matrix[6] = -3;		matrix[7] = -4;
	matrix[8] = -6;			matrix[9] = 35;		matrix[10] = -1;	 matrix[11] = -6;	matrix[12] = -5;	 matrix[13] = -6;	matrix[14] = -3;	matrix[15] = -8;
	matrix[16] = -1;		matrix[17] = -1;	matrix[18] = 19;	 matrix[19] = -6;	matrix[20] = -8;	 matrix[21] = -2;	matrix[22] = 0;		matrix[23] = -1;
	matrix[24] = -6;		matrix[25] = -6;	matrix[26] = -6;	 matrix[27] = 36;	matrix[28] = -4;	 matrix[29] = -3;	matrix[30] = -4;	matrix[31] = -7;
	matrix[32] = -3;		matrix[33] = -5;	matrix[34] = -8;	 matrix[35] = -4;	matrix[36] = 25;	 matrix[37] = 0;	matrix[38] = -1;	matrix[39] = -4;
	matrix[40] = -4;		matrix[41] = -6;	matrix[42] = -2;	 matrix[43] = -3;	matrix[44] = 0;		 matrix[45] = 28;	matrix[46] = -8;	matrix[47] = -5;
	matrix[48] = -3;		matrix[49] = -3;	matrix[50] = 0;		 matrix[51] = -4;	matrix[52] = -1;	 matrix[53] = -8;	matrix[54] = 21;	matrix[55] = -2;
	matrix[56] = -4;		matrix[57] = -8;	matrix[58] = -1;	 matrix[59] = -7;	matrix[60] = -4;	 matrix[61] = -5;	matrix[62] = -2;	matrix[63] = 31;
}

void matrixBInit(double * matrix, double p)
{
	matrix[0] = 140 + 8 * p;
	matrix[1] = -91;
	matrix[2] = -7;
	matrix[3] = 142;
	matrix[4] = 7;
	matrix[5] = -99;
	matrix[6] = 25;
	matrix[7] = -117;
}

void inverse_matrix(double * matrix, double * matrix1, double * cond) // получение обратной матрицы
{
	int pivot[mdim], flag;
	double matrix_copy[mdim * mdim];
	double b[mdim] = { 0 };

	for (int i = 0; i < mdim; i++)
	{
		for (int j = 0; j < mdim; j++)
		{
			matrix_copy[INDX(i, j)] = matrix[INDX(i, j)];
		}
	}

	decomp(mdim, mdim, matrix_copy, cond, pivot, &flag);

	for (int k = 0; k < mdim; k++)
	{
		b[k] = 1;
		solve(mdim, mdim, matrix_copy, b, pivot);
		for (int i = 0; i < mdim; i++) {
			matrix1[INDX(i, k)] = b[i];
			b[i] = 0;
		}
	}
}

void mult_matrix(double *matrix_a, double *matrix_b, double * matrix_x)
{
	double temp;
	for (int i = 0; i < mdim; ++i)
	{
		temp = 0;
		for (int j = 0; j < mdim; ++j)
		{
			temp += matrix_a[INDX(i, j)] * matrix_b[j];
		}
		matrix_x[i] = temp;
	}
}

void print_matrix_b(double * matrix, ofstream & stream)
{
	stream << endl;
	for (int i = 0; i < mdim; ++i)
	{
		stream << matrix[i] << endl;
	}
}

void print_matrix_x(double * matrix1,double *matrix2,ofstream & stream)
{
	stream << endl;
	for (int i = 0; i < mdim; ++i)
	{
		stream << matrix1[i]<< endl;
	}
	stream << endl;
	for (int i = 0; i < mdim; ++i)
	{
		stream << matrix2[i]  << endl;
	}
	stream << endl;
	for (int i = 0; i < mdim; ++i)
	{
		stream << abs(matrix1[i] - matrix2[i]) << endl;
	}
	stream << endl;
}

void print_matrix(double *matrix, ofstream & stream) // печать матрицы
{
	
	stream << endl;
	for (int i = 0; i < mdim; ++i)
	{
		for (int j = 0; j < mdim; ++j)
		{
			stream << matrix[INDX(i, j)] << "  ";
		}
		stream << endl;
	}
}

int main()
{
	double matrix_a[mdim * mdim];
	double matrix_a_inverted[mdim * mdim];
	double b_matrix[mdim];
	double x[mdim];
	double p[5] = { 1.0,0.1,0.01,0.0001,0.000001 };
	int pivot[mdim];
	double cond1,cond2;
	int flag;
	double maxx, maxerr;

	ofstream fin("output.txt");
	
	fin.setf(ios::fixed);


	for (int i = 0; i < 5; ++i)
	{
		matrixInit(matrix_a, p[i]);
		matrixBInit(b_matrix, p[i]);
		fin.precision(6);

		fin << "P = " << p[i] << endl<<endl;
		
		//fin  << "   Matrix A" << endl;
		//print_matrix(matrix_a,fin);
		//fin << endl  << "   Matrix B" << endl;
		//print_matrix_b(b_matrix,fin);

		inverse_matrix(matrix_a, matrix_a_inverted,&cond1);

		//fin << endl << "   Matrix A^(-1)" << endl;
		//print_matrix(matrix_a_inverted,fin);

		mult_matrix(matrix_a_inverted, b_matrix, x);
		decomp(mdim, mdim, matrix_a, &cond2, pivot, &flag);
		solve(mdim, mdim, matrix_a, b_matrix, pivot);
		fin.precision(20);
		fin  << endl << "First x        Second x      First x - Second x" << endl;
		print_matrix_x(x,b_matrix,fin);
		fin.precision(20);
		fin << endl << "Cond 1    Cond 2" << endl;
		fin << cond1 << "    " << cond2 << endl<<endl;

		maxx = x[0];
		maxerr = abs(x[0] - b_matrix[0]);
		for (int i = 0; i < mdim; ++i)
		{
			if (maxx < x[i])
			{
				maxx = x[i];
			}
			if (maxerr < abs(x[i] - b_matrix[i]))
			{
				maxerr = abs(x[i] - b_matrix[i]);
			}
		}
		fin.precision(20);
		fin << "d=" << maxerr / maxx << endl << endl;
	}
	//system("pause");
	return 0;

}
