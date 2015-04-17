#include "QSMatrix.h"
#ifndef QSMATRIX_CPP
#define QSMATRIX_CPP
// Parameter Constructor
template<typename T>
QSMatrix<T>::QSMatrix(unsigned _rows, unsigned _cols, const T& _initial) {
/*
  mat.resize(_rows);
  for (unsigned i=0; i<mat.size(); i++) {
    mat[i].resize(_cols, _initial);
  }
*/
    mat.reserve(_rows);
 for(unsigned i=0;i<_rows;++i){
        std::vector<T> cols;
        cols.reserve(_cols);
        for(unsigned j=0;j<_cols;++j){

            cols.push_back(_initial);
        }
        mat.push_back(cols);
    }
  rows = _rows;
  cols = _cols;
}
//vector constructor
template<typename T>
 QSMatrix<T>::QSMatrix(const std::vector<std::vector<T> > & inVector) {
	rows=inVector.size();
	mat.reserve(inVector.size());
  for(unsigned int i=0; i<inVector.size(); ++i) {
		mat.push_back(inVector[i]);
		if(i==0)
			cols=inVector[i].size();
		}

}


// Copy Constructor
template<typename T>
 QSMatrix<T>::QSMatrix(const QSMatrix<T>& rhs) {
  mat = rhs.mat;
  rows = rhs.get_rows();
  cols = rhs.get_cols();
}

// (Virtual) Destructor
template<typename T>
 QSMatrix<T>::~QSMatrix() {}
// Assignment Operator
template<typename T>
 QSMatrix<T>& QSMatrix<T>::operator=(const QSMatrix<T>&& rhs) {

    if (this!=&rhs)
    {
        mat=std::move(rhs.mat);
        rows=std::move(rhs.rows);
        cols=std::move(rhs.cols);
        return *this;

    }
    else
        return *this;

 }
// Assignment Operator
template<typename T>
 QSMatrix<T>& QSMatrix<T>::operator=(const QSMatrix<T>& rhs) {
  if (&rhs == this)
    return *this;

  unsigned new_rows = rhs.get_rows();
  unsigned new_cols = rhs.get_cols();
  mat=rhs.getArrayRef();
/*
 mat.resize(new_rows);
 // for (unsigned i=0; i<mat.size(); i++) {
   // mat[i].resize(new_cols);
  //}

  for (unsigned i=0; i<new_rows; i++) {
    for (unsigned j=0; j<new_cols; j++) {
        mat[i].resize(new_cols);
      mat[i][j] = rhs(i, j);
    }
  }

*/
  rows = new_rows;
  cols = new_cols;

  return *this;
}

// Addition of two matrices
template<typename T>
QSMatrix<T> QSMatrix<T>::operator+(const QSMatrix<T>& rhs) {
  QSMatrix<T> result(rows, cols, 0.0);

  for (unsigned i=0; i<rows; i++) {
    for (unsigned j=0; j<cols; j++) {
      result(i,j) = this->mat[i][j] + rhs(i,j);
    }
  }

  return result;
}

// Cumulative addition of this matrix and another
template<typename T>
QSMatrix<T>& QSMatrix<T>::operator+=(const QSMatrix<T>& rhs) {
  unsigned rows = rhs.get_rows();
  unsigned cols = rhs.get_cols();

  for (unsigned i=0; i<rows; i++) {
    for (unsigned j=0; j<cols; j++) {
      this->mat[i][j] += rhs(i,j);
    }
  }

  return *this;
}

// Subtraction of this matrix and another
template<typename T>
QSMatrix<T> QSMatrix<T>::operator-(const QSMatrix<T>& rhs) {
  unsigned rows = rhs.get_rows();
  unsigned cols = rhs.get_cols();
  QSMatrix<T> result(rows, cols, 0.0);

  for (unsigned i=0; i<rows; i++) {
    for (unsigned j=0; j<cols; j++) {
      result(i,j) = this->mat[i][j] - rhs(i,j);
    }
  }

  return result;
}

// Cumulative subtraction of this matrix and another
template<typename T>
QSMatrix<T>& QSMatrix<T>::operator-=(const QSMatrix<T>& rhs) {
  unsigned rows = rhs.get_rows();
  unsigned cols = rhs.get_cols();

  for (unsigned i=0; i<rows; i++) {
    for (unsigned j=0; j<cols; j++) {
      this->mat[i][j] -= rhs(i,j);
    }
  }

  return *this;
}

// Left multiplication of this matrix and another
template<typename T>
inline QSMatrix<T> QSMatrix<T>::operator*(const QSMatrix<T> & B) {
  // std::cout << "this rows: " << rows << " cols: " << cols << std::endl;
  // std::cout << "rhs rows: " << B.get_rows() << " cols: " << B.get_cols()<< std::endl;
  // std::cout << std::endl;
  			 // std::unique_ptr<std::ofstream> outvm= FileIOHelper::openFileForWrite("vMinside.txt",false);
           // FileIOHelper::writeLine(outvm,print(1,5)) ;
            //			  std::unique_ptr<std::ofstream> outvm2= FileIOHelper::openFileForWrite("Xinside.txt",false);
           // FileIOHelper::writeLine(outvm2,B.print(1,5)) ;
 if (B.get_rows() != cols) {
         throw new std::runtime_error("Matrix inner dimensions must agree.");
      }
      QSMatrix<T> X(rows,B.get_cols(),0);
      std::vector<std::vector<double>> & C = X.getArrayRefEdit();
     // std::vector<double> Bcolj;
     // Bcolj.resize(cols);

      for (int j = 0; j < B.get_cols(); j++) {
         std::vector<double> Bcolj;
         Bcolj.reserve(cols);
         for (int k = 0; k < cols; k++) {
          //  Bcolj[k] = B.mat[k][j];
            T value=B.mat[k][j];
            Bcolj.push_back(value);
         //   std::cout << "bcolj1" << Bcolj[k] << std::endl;
         }
         for (int i = 0; i < rows; i++) {
            const std::vector<double> &Arowi = mat[i];

            double s = 0;
            for (int k = 0; k < cols; k++) {
            //  double temp= Arowi[k];
              //  double tmp2= Bcolj[k] ;
               s += Arowi[k]*Bcolj[k];

            }
            C[i][j] = s;
         //   std::cout << C[i][j] << std::endl;
         }
      }
      return X;
}
/*
QSMatrix<T> QSMatrix<T>::operator*(const QSMatrix<T>& rhs) {
  unsigned rows = rhs.get_rows();
  unsigned cols = rhs.get_cols();
  QSMatrix<T> result(rows, cols, 0.0);
   std::cout << "in operator *" << std::endl;
  for (unsigned i=0; i<rows; i++) {
    for (unsigned j=0; j<cols; j++) {
      for (unsigned k=0; k<rows; k++) {
        result(i,j) += this->mat[i][k] * rhs(k,j);
      }
    }
  }
   std::cout << "done w mult" << std::endl;
  return result;
}
*/

// Cumulative left multiplication of this matrix and another
template<typename T>
QSMatrix<T>& QSMatrix<T>::operator*=(const QSMatrix<T>& rhs) {
  QSMatrix<T> result = (*this) * rhs;
  (*this) = result;
  return *this;
}

// Calculate a transpose of this matrix
template<typename T>
 QSMatrix<T> QSMatrix<T>::transpose() {
  QSMatrix<T> result(cols, rows, 0.0);
 // std::cout << "this: rows: " << rows << "cols: " << cols << std::endl;
 // std::cout << "result: rows: " << result.get_rows() << "cols: " << result.get_cols() << std::endl;

  for (unsigned i=0; i<rows; i++) {
    for (unsigned j=0; j<cols; j++) {
      result(j,i) = this->mat[i][j];
    }
  }
 //  std::cout << "done w transpose" << std::endl;
  return result;
}

// Matrix/scalar addition
template<typename T>

QSMatrix<T> QSMatrix<T>::operator+(const T& rhs) {
  QSMatrix<T> result(rows, cols, 0.0);

  for (unsigned i=0; i<rows; i++) {
    for (unsigned j=0; j<cols; j++) {
      result(i,j) = (T)this->mat[i][j] + rhs;
    }
  }

  return result;
}

// Matrix/scalar subtraction
template<typename T>
QSMatrix<T> QSMatrix<T>::operator-(const T& rhs) {
  QSMatrix<T> result(rows, cols, 0.0);

  for (unsigned i=0; i<rows; i++) {
    for (unsigned j=0; j<cols; j++) {
      result(i,j) = this->mat[i][j] - rhs;
    }
  }

  return result;
}

// Matrix/scalar multiplication
template<typename T>
 QSMatrix<T> QSMatrix<T>::operator*(const T& rhs) {
  QSMatrix<T> result(rows, cols, 0.0);

  for (unsigned i=0; i<rows; i++) {
    for (unsigned j=0; j<cols; j++) {
      result(i,j) = (T)this->mat[i][j] * rhs;
    }
  }

  return result;
}

// Matrix/scalar division
template<typename T>
 QSMatrix<T> QSMatrix<T>::operator/(const T& rhs) {
  QSMatrix<T> result(rows, cols, 0.0);

  for (unsigned i=0; i<rows; i++) {
    for (unsigned j=0; j<cols; j++) {
      result(i,j) = (T)this->mat[i][j] / rhs;
    }
  }

  return result;
}

// Multiply a matrix with a vector
template<typename T>
 std::vector<T> QSMatrix<T>::operator*(const std::vector<T>& rhs) {
  std::vector<T> result(rhs.size(), 0.0);

  for (unsigned i=0; i<rows; i++) {
    for (unsigned j=0; j<cols; j++) {
      result[i] = (T)this->mat[i][j] * rhs[j];
    }
  }

  return result;
}

// Obtain a vector of the diagonal elements
template<typename T>
std::vector<T> QSMatrix<T>::diag_vec() {
  std::vector<T> result(rows, 0.0);

  for (unsigned i=0; i<rows; i++) {
    result[i] = this->mat[i][i];
  }

  return result;
}

// Access the individual elements
template<typename T>
 T& QSMatrix<T>::operator()(const unsigned& row, const unsigned& col) {
  return this->mat[row][col];
}

// Access the individual elements (const)
template<typename T>
const T& QSMatrix<T>::operator()(const unsigned& row, const unsigned& col) const {
  return this->mat[row][col];
}
// Access the individual elements (const)
template<typename T>
const std::vector<T> & QSMatrix<T>::operator()(const unsigned& row) const{
  return this->mat[row];
}
// Get the number of rows of the matrix
template<typename T>
 unsigned QSMatrix<T>::get_rows() const {
  return this->rows;
}

// Get the number of columns of the matrix
template<typename T>
 unsigned QSMatrix<T>::get_cols() const {
  return this->cols;
}
// return x = A^-1 b, assuming A is square and has full rank
template<typename T>
 void QSMatrix<T>::set(int i, int j, T s){
	mat[i][j]=s;
}

template<typename T>
 T QSMatrix<T>::get(int i, int j){
	return mat[i][j];
}

template<typename T>
T QSMatrix<T>::norm1() {
      T f = 0;
      for (int j = 0; j < cols; j++) {
         T s = 0;
         for (int i = 0; i < rows; i++) {
            s += std::abs(mat[i][j]);
         }
         f = std::max(f,s);
      }
      return f;
}

template<typename T>
QSMatrix<T>  QSMatrix<T>::minusInPlace (QSMatrix<T> & B)
   {
	   for (int i = 0; i < rows; i++)
	   {
	   	 for (int j = 0; j < cols; j++)
	   	 {
	   		 mat[i][j] = mat[i][j] - B.mat[i][j];
	   	 }
	   }
	   return mat;
   }

template<typename T>
QSMatrix<T> QSMatrix<T>::identity (int m, int n) {
      QSMatrix<T> A(m,n,0);
      std::vector<std::vector<double>> & X = A.getArrayRefEdit();
      for (int i = 0; i < m; i++) {
         for (int j = 0; j < n; j++) {
            X[i][j] = (i == j ? 1.0 : 0.0);
         }
      }
      return A;
   }

template<typename T>
std::vector<double> QSMatrix<T>::getColumnPackedCopy () {
      std::vector<double> vals;
	  vals.resize(rows*cols);
      for (int i = 0; i < rows; i++) {
         for (int j = 0; j < cols; j++) {
            vals[(i+j*rows)] = mat[i][j];
         }
      }
      return vals;
}
template<typename T>
QSMatrix<double> QSMatrix<T>::getMatrix (const std::vector<int> & r, int j0, int j1) const{
      QSMatrix<double> X (r.size(),j1-j0+1,0);
      std::vector<std::vector<double>> & B = X.getArrayRefEdit();
      try {
         for (int i = 0; i < r.size(); i++) {
            for (int j = j0; j <= j1; j++) {
               B[i][j-j0] = mat[r[i]][j];
            }
         }
      } catch(std::exception & e) {
         throw std::runtime_error("Submatrix indices");
      }
      return X;
   }
template<typename T>
QSMatrix<T> QSMatrix<T>::getMatrix (const std::vector<T> & r, int j0, int j1) const{
      QSMatrix<T> X (r.size(),j1-j0+1,0);
      std::vector<std::vector<T>> B = X.getArrayRefEdit();
      try {
         for (int i = 0; i < r.size(); i++) {
            for (int j = j0; j <= j1; j++) {
               B[i][j-j0] = mat[r[i]][j];
            }
         }
      } catch(std::exception & e) {
         throw std::runtime_error("Submatrix indices");
      }
      return X;
   }

template<typename T>
QSMatrix<T> QSMatrix<T>::getMatrix (int i0, int i1, int j0, int j1) {
      QSMatrix X (i1-i0+1,j1-j0+1,0);
      std::vector<std::vector<double>> B = X.getArrayRef();
      try {
         for (int i = i0; i <= i1; i++) {
            for (int j = j0; j <= j1; j++) {
               B[i-i0][j-j0] = mat[i][j];
            }
         }
      } catch(std::exception & e) {
         throw std::runtime_error("Submatrix indices");
      }
      return X;
   }

template<typename T>
QSMatrix<T>::QSMatrix (const std::vector<std::vector<T>> & A, int m, int n) {
      mat = A;
      rows = m;
      cols = n;
   }

template<typename T>
QSMatrix<T>::QSMatrix (const std::vector<double> &vals, int m) {
      rows = m;
      cols = (m != 0 ? vals.size()/m : 0);
      if (m*cols != vals.size()) {
         throw std::runtime_error("Array length must be a multiple of m.");
      }
      std::vector<std::vector<double>>A ;
	  A.resize(m);
      for (int i = 0; i < m; i++) {
         for (int j = 0; j < cols; j++) {
            A[i].push_back(vals[i+j*m]);
         }
      }
   }

template<typename T>
std::string QSMatrix<T>::print (int padding, int width) const {
      std::stringstream output;
      for (int i = 0; i < get_rows(); i++) {
         for (int j = 0; j < get_cols(); j++) {
            std::stringstream s;
               s << mat[i][j]; // format the number
            std::string tmp =s.str();
            int padding = width-tmp.length(); // At _least_ 1 space
            for (int k = 0; k < padding; k++)
               output <<" ";
            output<<tmp;
         }
         output << "\n" ;
      }
         output << "\n";
      return output.str();
}
#endif
