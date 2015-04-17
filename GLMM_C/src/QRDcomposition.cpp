#include "QRDcomposition.h"
QRDcomposition::QRDcomposition (const QSMatrix<double>& A) {
      // Initialize.
      QR = A.getArrayRef();
      m = A.get_rows();
      n = A.get_cols();
      Rdiag ;
	  Rdiag.resize(n);

      // Main loop.
      for (int k = 0; k < n; k++) {
         // Compute 2-norm of k-th column without under/overflow.
         double nrm = 0;
         for (int i = k; i < m; i++) {
            nrm = std::hypot(nrm,QR[i][k]);
         }

         if (nrm != 0.0) {
            // Form k-th Householder vector.
            if (QR[k][k] < 0) {
               nrm = -nrm;
            }
            for (int i = k; i < m; i++) {
               QR[i][k] /= (double)nrm;
            }
            QR[k][k] += 1.0;

            // Apply transformation to remaining columns.
            for (int j = k+1; j < n; j++) {
               double s = 0.0;
               for (int i = k; i < m; i++) {
                  s += QR[i][k]*QR[i][j];
               }
               s = -s/QR[k][k];
               for (int i = k; i < m; i++) {
                  QR[i][j] += s*QR[i][k];
               }
            }
         }
         Rdiag[k] = -nrm;
      }
   }
QSMatrix<double> QRDcomposition::solve (const QSMatrix<double>& B) {
      if (B.get_rows() != m) {
         throw std::runtime_error("Matrix row dimensions must agree.");
      }
      if (!isFullRank()) {
         throw std::runtime_error("Matrix is rank deficient.");
      }

      // Copy right hand side
      int nx = B.get_cols();
      std::vector<std::vector<double>> X = B.getArrayRef();

      // Compute Y = transpose(Q)*B
      for (int k = 0; k < n; k++) {
         for (int j = 0; j < nx; j++) {
            double s = 0.0;
            for (int i = k; i < m; i++) {
               s += QR[i][k]*X[i][j];
            }
            s = -s/QR[k][k];
            for (int i = k; i < m; i++) {
               X[i][j] += s*QR[i][k];
            }
         }
      }
      // Solve R*X = Y;
      for (int k = n-1; k >= 0; k--) {
         for (int j = 0; j < nx; j++) {
            X[k][j] /= Rdiag[k];
         }
         for (int i = 0; i < k; i++) {
            for (int j = 0; j < nx; j++) {
               X[i][j] -= X[k][j]*QR[i][k];
            }
         }
      }
      return (QSMatrix<double>(X,n,nx).getMatrix(0,n-1,0,nx-1));
   }

bool QRDcomposition::isFullRank () {
      for (int j = 0; j < n; j++) {
         if (Rdiag[j] == 0)
            return false;
      }
      return true;
   }


