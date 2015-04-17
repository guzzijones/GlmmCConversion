#include "TriangleMatrix.h"
#ifndef TRIANGLEMATRIX_CPP
#define TRIANGLEMATRIX_CPP
template <typename T>
void TriangleMatrix<T>::solveLForLDLT(const std::vector<std::vector<T>> & L, std::vector<T> & b , int n ){
		T sum = 0;
    	for( int i = 0; i < n; i++ )
    	{
    	  sum = b[i];
    	  for( int k = 0; k < i; k++)
    	  	  sum -= L[i][k] * b[k];
    	  b[i] = sum;
    	}

}
template <typename T>
void TriangleMatrix<T>::solveUForLDLT( const std::vector<std::vector<double>> & U, std::vector<double> & b , int n )
    {
        for( int i = n - 1; i >= 0; i-- )
        {
            double sum = b[i];
            for( int j = i+1; j <n; j++ )
                 sum -= U[i][j] * b[j];
            b[i] = sum;
        }
    }
#endif
