#ifndef TRIANGLEMATRIX_H
#define TRIANGLEMATRIX_H
#include <vector>
template <typename T>
class TriangleMatrix
{
	private:

	protected:

	public:
	static void solveLForLDLT(const std::vector<std::vector<T>> & L, std::vector<T> & b , int n );
	static void solveUForLDLT( const std::vector<std::vector<double>> & U, std::vector<double>  & b , int n );



};

#include "../src/TriangleMatrix.cpp"
#endif
