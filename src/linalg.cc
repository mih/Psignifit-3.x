#include "linalg.h"

Matrix::Matrix ( const std::vector< std::vector<double> >& A )
	: nrows(A.size()), ncols(A[0].size())
{
	data = new double [nrows*ncols];

	int i,j;
	for (i=0; i<nrows; i++) {
		for (j=0; j<ncols; j++) {
			this->operator() (i,j) = A[i][j];
		}
	}
}

Matrix::Matrix ( unsigned int n, unsigned int m )
	: nrows(n), ncols(m)
{
	data = new double [nrows*ncols];
	int i;
	for (i=0; i<nrows*ncols; i++)
		data[i] = 0;
}

Matrix::Matrix ( const Matrix& A )
	: nrows(A.getnrows()), ncols(A.getncols())
{
	data = new double [nrows*ncols];

	int i,j;
	for (i=0; i<nrows; i++) {
		for (j=0; j<ncols; j++) {
			this->operator() (i,j) = A(i,j);
		}
	}
}

double& Matrix::operator() ( unsigned row, unsigned col ) const
{
	//row --; col --;
	if (row>=nrows || col>=ncols)
		throw MatrixError();
	return data[row+nrows*col];
}

void Matrix::print ( void ) {
	int i,j;
	std::cout << "[ ";
	for (i=0; i<nrows; i++) {
		std::cout << "[";
		for (j=0; j<ncols; j++) {
			std::cout << " " << std::setprecision(3) << std::setw(5) << data[i+nrows*j] << (j!=ncols-1 ? " , " : (i==nrows-1 ? "] ]\n" : "],\n  ") );
		}
	}
}

Matrix* Matrix::cholesky_dec ( void ) {
	if (nrows!=ncols)
		throw MatrixError();

	Matrix *L = new Matrix (nrows,nrows);
	int k,j,K;

	for (K=0; K<nrows; K++) {
		// l_{KK} = (a_{KK} - \sum_{k=1}^{K-1} l_{KK}^2)^{1/2}
		(*L)(K,K) = this->operator() (K,K);
		for (k=0; k<K; k++)
			(*L)(K,K) -= (*L)(K,k) * (*L)(K,k);
		(*L)(K,K) = sqrt((*L)(K,K));

		// l_{jK} = (a_{jK} - \sum_{k=1}^{K-1} l_{jk} l_{Kk})/l_{KK}, j >= K+1
		for (j=K+1; j<nrows; j++) {
			(*L)(j,K) = this->operator() (j,K);
			for (k=0; k<K; k++) {
				(*L)(j,K) -= (*L)(j,k)* (*L)(K,k);
			}
			(*L)(j,K) /= (*L)(K,K);
		}
	}

	return L;
}

Matrix* Matrix::lu_dec ( void ) {
	if (nrows!=ncols)
		throw MatrixError ();

	Matrix *LU = new Matrix (*this);

	int i,j,k;
	double c;
	int pivotindex;
	double pivot;

	for (i=0; i<nrows-1; i++) {
		for (k=i+1; k<nrows; k++) {
			c = (*LU)(k,i) / (*LU)(i,i);
			(*LU)(k,i) = c;
			for (j=i+1; j<nrows; j++)
				(*LU)(k,j) = (*LU)(k,j) - c * (*LU)(i,j);
		}
	}

	return LU;
}

std::vector<double> Matrix::solve ( const std::vector<double>& b ) {
	Matrix *LU = lu_dec();

	std::vector<double> x ( nrows );
	std::vector<double> y ( nrows );
	int i,k;
	double s;

	y = forward ( LU, b);
	x = backward ( LU, y );

	delete LU;

	return x;
}

std::vector<double> Matrix::forward ( const Matrix *LU, const std::vector<double>& b ) {
	int i,k;
	double s;
	std::vector<double> y (nrows);

	for (i=0; i<nrows; i++) {
		s = b[i];
		for (k=0; k<i; k++) {
			s -= (*LU)(i,k) * y[k];
		}
		y[i] = s;
	}
	return y;
}

std::vector<double> Matrix::backward ( const Matrix *LU, const std::vector<double>& y ) {
	int i,k;
	double s;
	std::vector<double> x (nrows);

	for (i=nrows-1; i>=0; i--) {
		s = y[i];
		for (k=i+1; k<nrows; k++) {
			s -= (*LU)(i,k) * x[k];
		}
		x[i] = s/(*LU)(i,i);
	}
	return x;
}

Matrix* Matrix::inverse ( void ) {
	if (nrows!=ncols)
		throw MatrixError();

	Matrix *LU = lu_dec ();
	Matrix *Inv = new Matrix ( nrows, nrows );
	std::vector<double> x ( nrows, 0);
	std::vector<double> y ( nrows, 0);
	int i,j;

	for ( i=0; i<ncols; i++ ) {
		for (j=0; j<nrows; j++)
			x[j] = 0;
		x[i] = 1;

		y = forward ( LU, x );
		x = backward ( LU, y );

		for (j=0; j<nrows; j++)
			(*Inv)(j,i) = x[j];
	}

	delete LU;

	return Inv;
}
