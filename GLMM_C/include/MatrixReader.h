#ifndef MATRIXREADER_H
#define MATRIXREADER_H
#include "FileIOHelper.h"
#include <string>
#include <vector>
#include <memory>
#include "nstring.h"
#include <sstream>

class MatrixReader {
	private:
		std::string _fileName;
		std::vector< std::vector<double> > _X;
	protected:

	public:

	MatrixReader(std::string in);
	const std::vector < std::vector<double> > & getData()const;
	void read();


	//accessors


};
#endif
