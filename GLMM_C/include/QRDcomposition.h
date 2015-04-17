#ifndef QRD_H
#define QRD_H
#include <vector>
#include <stdexcept>
#include <cmath>
#include "QSMatrix.h"


class QRDcomposition{
	private:
	static const long serialVersionUID = 1L;
   std::vector<std::vector<double>> QR;

   /** Row and column dimensions.
   @serial column dimension.
   @serial row dimension.
   */
   int m, n;

   /** Array for internal storage of diagonal of R.
   @serial diagonal of R.
   */
   std::vector<double> Rdiag;


public:

	QSMatrix<double> solve (const QSMatrix<double> & B) ;
	QRDcomposition (const QSMatrix<double>& A); 
	bool isFullRank () ;
};

#endif
