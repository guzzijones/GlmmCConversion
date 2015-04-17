#include "LDLT.h"
#include "FileIOHelper.h"
#include "QSMatrix.h"
LDLT::LDLT(std::vector<std::vector<double>> & aM):A(aM){

	n = A.size();//The order of matrix A.

}

void LDLT::decompose(){
        //    std::unique_ptr<std::ofstream> outvm3= FileIOHelper::openFileForWrite("vMindecompose.txt",false);
          //  QSMatrix<double> QA(A);
            //FileIOHelper::writeLine(outvm3,QA.print(1,20)) ;
        D.resize(n);
		std::vector<double> v;
		v.resize(n);
		for(int j = 0; j < n; j++)
		{
			if(j > 0)
			{
				for(int k = 0; k <= j - 1; k++)
					v[k] = A[j][k] * D[k];

				double sum = 0;
				for(int k = 0; k <= j - 1; k++)
					sum += A[j][k] * v[k];
				v[j] = A[j][j] - sum;
				D[j] = v[j];

				if(j < n - 1)
					for(int k = j + 1; k < n; k++)
					{
						sum = 0;
						for(int w = 0; w <= j - 1; w++)
							sum += A[k][w] * v[w];
						A[k][j] = (A[j][k] - sum) / (double)v[j];
					}
			}
			else
			{
				//Initialize...checked.
				v[0] = A[0][0];
				D[0] = v[0];
				for(int k = 1; k < n; k++)
					A[k][0] /= (double)v[0];
			}
		}//End of this loop

 //std::unique_ptr<std::ofstream> outvm5= FileIOHelper::openFileForWrite("vMbeforesetIndompose.txt",false);
   //         QSMatrix<double> QA5(A);
     //       FileIOHelper::writeLine(outvm5,QA5.print(1,20)) ;
		for(int i = 0; i < n; i++){
			A[i][i] = 1.0;
			}

// std::unique_ptr<std::ofstream> outvm4= FileIOHelper::openFileForWrite("vMinsidefuncendompose.txt",false);
  //          QSMatrix<double> QA2(A);
 //           FileIOHelper::writeLine(outvm4,QA2.print(1,20)) ;

}

std::vector<std::vector<double>> LDLT::inverse(){
	//Copy lower triangle data to upper triangle of matrix A.
		//So, the lower triangle part is L
		//The upper triangle part is L' which will be used later.
		//The diagonal of matrix A is 1.
		for (int j = 0; j < n; j++)
			for (int i = j + 1; i < n; i++)
			    A[j][i] = A[i][j];

		//Solve AZ = B with B is an n x n identify matrix.
		//This is equivalent to invert matrix A, but
		//by using solving we save computational time.
		std::vector<double> b;
		b.resize(n);
		for(int i = 0; i < n; i++)
		{
			for(int j = 0; j < n; j++)
				b[j]=(j == i)?1.0:0.0;
				TriangleMatrix<double>::solveLForLDLT(A, b, n);//Update b.
			//Note here we only update the lower triangle of A.
			for(int j = i; j < n; j++)
				A[j][i] = b[j];
		}

		//Solve DY = Z;
		//This is the fastest way of solving this equation given
		//D is diagonal matrix.
	   for (int i = 0; i < n; i++)
	      for(int j = 0; j <= i; j++)
	         A[i][j] /= (double)D[i];//Only update the lower triangle of matrix A.

	   //Solve L'X = Y with L' saved as upper triangle of A matrix.
	   	for(int i = 0; i < n; i++)
		{
		    for(int j = 0; j < n; j++)
			    b[j] = (j < i)? 0.0: A[j][i];
		   	TriangleMatrix<double>::solveUForLDLT(A, b, n);//Update b.
			for(int j = i; j < n; j++)
				A[j][i] = b[j];//Update upper triangle of matrix A.
		}

		//Make A a symmetric matrix.
		for (int j = 0; j < n; j++)
		   for (int i = j + 1; i < n; i++)
		     A[j][i] = A[i][j];
	   return A;

}
std::vector<double> LDLT::getD(){
return D;
}
double LDLT::getLodDet(){
    double logD = 0; //Math.log(d);
	    for (int j = 0; j < n; j++)
	      if(D[j] > 0)
	    	  logD += std::log(D[j]);
	    return logD;
}

