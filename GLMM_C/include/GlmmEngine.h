#ifndef GLMMENGINE_H
#define GLMMENGINE_H
#include <string>
#include <vector>
#include "MatrixReader.h"
#include "QSMatrix.h"
#include "Glmm.h"
#include "time.h"
class GlmmEngine{
	private:

	protected:

	public:

	
	static void doAnalysis(std::string YFielName, 
			std::string XFielName, 
			std::string RFielName, 
			std::vector<std::string> zFileNames, 
			std::vector<std::string> gFileNames, 
			std::string outputFileName) ;
	
	

};

#endif

