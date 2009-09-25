#ifndef LINALG_H
#define LINALG_H

#include <vector>
#include <iostream>
#include <iomanip>
#include <cmath>

/***************************************************************************************
 * This file defines vary basic linear algebra facilities
 *
 * These are not particularly fast or fancy implementations. They are just intended
 * to keep the number of dependencies low.
 */

class MatrixError
{
};

/** A very simple matrix class */
class Matrix
{
	private:
		double *data;
		unsigned int nrows;
		unsigned int ncols;
		// given a decomposition A=LU these two methods help solving the equation Ax=b:
		std::vector<double> forward ( const Matrix* LU, const std::vector<double>& b );    // forward solution of Ly=b
		std::vector<double> backward ( const Matrix*LU, const std::vector<double>& y );    // backward solution of Ux=y
	public:
		Matrix ( const std::vector< std::vector<double> >& A );      ///< Construct a matrix from a vector of vectors
		Matrix ( unsigned int n, unsigned int m );                   ///< Construct a matrix with given dimensions initialized to 0
		Matrix ( const Matrix& A );                                  ///< copy a matrix
		~Matrix ( void ) { delete [] data; }                         ///< delete a matrix
		double& operator() ( unsigned i, unsigned j ) const;         ///< data access to the element in row i and column j (indices starting with 0!)
		void print ( void );                                         ///< print the matrix to stdout
		unsigned int getnrows ( void ) const { return nrows; }       ///< get the number of rows
		unsigned int getncols ( void ) const { return ncols; }       ///< get the number of columns
		Matrix* cholesky_dec ( void ) const;                         ///< return a pointer to a newly allocated Matrix L, such that L*L^T = A
		Matrix* lu_dec ( void ) const;                               ///< return a pointer to a newly allocated Matrix LU. If you construct a matrix L from the elements of LU below the diagnonal and L_ii = 1, and a matrix U from the elements above the diagnonal, then A=LU
		std::vector<double> solve ( const std::vector<double>& b );  ///< solve a linear equation Ax=b for x
		Matrix* inverse ( void );                                    ///< return a pointer to the (newly allocated) inverse Matrix
		std::vector<double> operator* ( std::vector<double>& x );    ///< right multiplication with a vector
		void scale ( double a );                                     ///< scale the whole matrix by a constant factor
		bool symmetric ( void );                                     ///< check whether the matrix is symmetric
};

#endif
