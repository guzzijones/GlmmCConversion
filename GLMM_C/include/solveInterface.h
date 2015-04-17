#ifndef SOLVEINTERFACE_H
#define SOLVEINTERFACE_H
#include "QSMatrix.h"
#include "LUDcomposition.h"
#include "QRDcomposition.h"
#include <iostream>

class solveInterface{
	private:
		QSMatrix<double> A;
		int rows;
		int cols;

	protected:

	public:
	solveInterface(const QSMatrix<double> & in){
	A=in;
	rows=A.get_rows();
	cols=A.get_cols();
	};

	QSMatrix<double> solve(const QSMatrix<double> & rhs) {

		return (rows == cols ? (LUDcomposition(A)).solve(rhs) :
                       (QRDcomposition(A)).solve(rhs));

      }
	QSMatrix<double> inverse () {
    //  std::cout << "testbefore:" <<A.print(1,20) << std::endl;

      //  QSMatrix<double> test=A.identity(rows,rows);
     //   std::unique_ptr<std::ofstream> hf= FileIOHelper::openFileForWrite("testinverse.txt",false);
       //     FileIOHelper::writeLine(hf,test.print(1,20)) ;
      //std::cout << "test:" <<test.print(1,20) << std::endl;
      return solve(A.identity(rows,rows));
   }





};
#endif
