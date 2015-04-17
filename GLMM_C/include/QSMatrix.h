#ifndef __QS_MATRIX_H
#define __QS_MATRIX_H
#include <vector>
#include <stdexcept>
#include <cmath>
#include <iostream>
#include <sstream>
#include <string>
#include "FileIOHelper.h"
template <typename T> class QSMatrix {
 private:
  std::vector<std::vector<T> > mat;
  unsigned rows;
  unsigned cols;

 public:

	static QSMatrix<T> identity (int m, int n) ;
    QSMatrix<double> getMatrix (const std::vector<int> & r, int j0, int j1)const;
	QSMatrix<T> getMatrix (const std::vector<T> & r, int j0, int j1)const;
	QSMatrix<T> getMatrix (int i0, int i1, int j0, int j1) ;
	QSMatrix (const std::vector<double> &vals, int m) ;
    QSMatrix(QSMatrix<T> && in)  { mat=std::move(in.mat);rows=std::move(in.rows);cols=std::move(in.cols);}
 std::vector<double> getColumnPackedCopy () ;
  QSMatrix<T> inverse () ;
	QSMatrix(const std::vector<std::vector<T>> & A, int m, int n) ;
  QSMatrix(unsigned _rows, unsigned _cols, const T& _initial);
  QSMatrix(const std::vector<std::vector<T> > & _vector);
  QSMatrix(const QSMatrix<T>& rhs);
  virtual ~QSMatrix();
  QSMatrix(){};

  // Operator overloading, for "standard" mathematical matrix operations
  QSMatrix<T>& operator=(const QSMatrix<T>& rhs);

  // Matrix mathematical operations
  QSMatrix<T> operator+(const QSMatrix<T>& rhs);
  QSMatrix<T>& operator+=(const QSMatrix<T>& rhs);
  QSMatrix<T> operator-(const QSMatrix<T>& rhs);
  QSMatrix<T>& operator-=(const QSMatrix<T>& rhs);
  QSMatrix<T> operator*(const QSMatrix<T>& rhs);
  QSMatrix<T>& operator*=(const QSMatrix<T>& rhs);
  QSMatrix<T>& operator=(const QSMatrix<T> && rhs);
  QSMatrix<T> transpose();

  // Matrix/scalar operations
  QSMatrix<T> operator+(const T& rhs);
  QSMatrix<T> operator-(const T& rhs);
  QSMatrix<T> operator*(const T& rhs);
  QSMatrix<T> operator/(const T& rhs);

  // Matrix/vector operations
  std::vector<T> operator*(const std::vector<T>& rhs);
  std::vector<T> diag_vec();

  // Access the individual elements
  T& operator()(const unsigned& row, const unsigned& col);
  const std::vector<T> & operator()(const unsigned& row) const;
  const T& operator()(const unsigned& row, const unsigned& col) const;


  // Access the row and column sizes
  unsigned get_rows() const;
  unsigned get_cols() const;
  const std::vector<std::vector<T> > & getArrayRef()const{const std::vector<std::vector<T>> & matref=mat;return matref;};
  std::vector<std::vector<T> > & getArrayRefEdit(){ std::vector<std::vector<T>> & matref=mat;return matref;};
  //mutators
  void set(int i, int j, T s);
  T get(int i, int j);
	T norm1();

std::string print (int padding=1, int width=1) const;
	QSMatrix<T>  minusInPlace(QSMatrix<T> & B);
};

#include "../src/QSMatrix.cpp"
#endif
