#include "LUDcomposition.h"
#include <iostream>
LUDcomposition::LUDcomposition(const QSMatrix<double> & A){
// Use a "left-looking", dot-product, Crout/Doolittle algorithm.
// Use a "left-looking", dot-product, Crout/Doolittle algorithm.

      LU = A.getArrayRef();

      m = A.get_rows();
      n = A.get_cols();
      piv.resize(m);
      for (int i = 0; i < m; i++) {
         piv[i] = i;
      }
      pivsign = 1;

      std::vector<double> LUcolj;
      LUcolj.resize(m);

      // Outer loop.

      for (int j = 0; j < n; j++) {

         // Make a copy of the j-th column to localize references.

         for (int i = 0; i < m; i++) {
            LUcolj[i] = LU[i][j];
         }

         // Apply previous transformations.

         for (int i = 0; i < m; i++) {
            std::vector<double> & LUrowi = LU[i];

            // Most of the time is spent in the following dot product.

            int kmax = std::min(i,j);
            double s = 0.0;
            for (int k = 0; k < kmax; k++) {
               s += LUrowi[k]*LUcolj[k];
            }

            LUrowi[j] = LUcolj[i] -= s;
         }

         // Find pivot and exchange if necessary.

         int p = j;
         for (int i = j+1; i < m; i++) {
            if (std::abs(LUcolj[i]) > std::abs(LUcolj[p])) {
               p = i;
            }
         }
         if (p != j) {
            for (int k = 0; k < n; k++) {
               double t = LU[p][k]; LU[p][k] = LU[j][k]; LU[j][k] = t;
            }
            int k = piv[p]; piv[p] = piv[j]; piv[j] = k;
            pivsign = -pivsign;
         }

         // Compute multipliers.

         if (j < m & LU[j][j] != 0.0) {
            for (int i = j+1; i < m; i++) {
               LU[i][j] /= (double)LU[j][j];
            }
         }
      }
       QSMatrix<double> tmp(LU);
//  std::unique_ptr<std::ofstream> out= FileIOHelper::openFileForWrite("LU.txt",false);
 //  FileIOHelper::writeLine(out,tmp.print(1,10)) ;
  //    std::cout << tmp.print(1,5);

   }

bool LUDcomposition::isNonsingular () {
      for (int j = 0; j < n; j++) {
		//	std::cout << LU[j][j] << std::endl;
         if (LU[j][j] == 0){
			std::cout << std::endl<< "Singular found: " << j << " value: "<< LU[j][j] << std::endl;

            return false;
			}
      }

      return true;
   }

QSMatrix<double> LUDcomposition::solve (const QSMatrix<double> & B) {

      if (B.get_rows() != m) {
         throw new std::runtime_error("Matrix row dimensions must agree.");
      }
      if (!isNonsingular()) {
         throw new std::runtime_error("Matrix is singular.");
      }

      // Copy right hand side with pivoting
      int nx = B.get_cols();
      QSMatrix<double> Xmat = B.getMatrix(piv,0,nx-1);
//        std::cout << "xmat:" << std::endl;
//        std::cout << Xmat.print(1,10);

      std::vector<std::vector<double>> & X = Xmat.getArrayRefEdit();
  //    std::unique_ptr<std::ofstream> out= FileIOHelper::openFileForWrite("X.txt",false);
//      std::cout << "X" << std::endl;
//     QSMatrix<double> Xtmp(X);
 //     std::cout << Xtmp.print(1,10);
//      FileIOHelper::writeLine(out,Xmat.print(1,10)) ;
      // Solve L*Y = B(piv,:)
      for (int k = 0; k < n; k++) {
         for (int i = k+1; i < n; i++) {
            for (int j = 0; j < nx; j++) {
               X[i][j] -= X[k][j]*LU[i][k];
            }
         }
      }
      // Solve U*X = Y;
      for (int k = n-1; k >= 0; k--) {
         for (int j = 0; j < nx; j++) {
            X[k][j] /= (double)LU[k][k];
         }
         for (int i = 0; i < k; i++) {
            for (int j = 0; j < nx; j++) {
               X[i][j] -= X[k][j]*LU[i][k];
            }
         }
      }

      return Xmat;

   }

