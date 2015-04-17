#ifndef LDLT_H
#define LDLT_H
#include "TriangleMatrix.h"
#include <vector>
#include <iostream>
#include <cmath>
class LDLT{
	private:
		std::vector<std::vector<double> > & A;
		std::vector<double> D;
		int n;


	protected:


	public:
	LDLT(std::vector<std::vector<double>> &aM);
	void decompose();
	std::vector<std::vector<double>> inverse();
	std::vector<double> getD();
	double getLodDet();
};

#endif

