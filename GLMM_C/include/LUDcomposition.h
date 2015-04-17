#ifndef LUD_H
#define LUD_H
#include <vector>
#include <stdexcept>
#include <cmath>
#include "QSMatrix.h"
#include "FileIOHelper.h"


class LUDcomposition{

	private:

	static const long serialVersionUID = 1L;
	std::vector<std::vector<double>> LU;
    int m, n, pivsign;

    std::vector<int> piv;

	protected:

	public:

	LUDcomposition(const QSMatrix<double> & A);
	QSMatrix<double> solve (const QSMatrix<double>& B) ;
	bool isNonsingular ();


};
#endif

