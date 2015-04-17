#include "GlmmEngine.h"
#include <fstream>

void GlmmEngine::doAnalysis(std::string YFielName,
			std::string XFielName,
			std::string RFielName,
			std::vector<std::string> zFileNames,
			std::vector<std::string> gFileNames,
			std::string outputFileName){

		std::cout << "Running REML to solve GLMM..." << std::endl;
      std::cout << "y file: " << YFielName << std::endl;
		MatrixReader dataReader(YFielName);



		dataReader.read();
		QSMatrix<double> Y(dataReader.getData());


		dataReader=MatrixReader(XFielName);
		dataReader.read();
		QSMatrix<double> X = dataReader.getData();
     // std::unique_ptr<std::ofstream> outX= FileIOHelper::openFileForWrite("Testx_inengine.txt",false);
     // FileIOHelper::writeLine(outX,X.print()) ;

		dataReader=MatrixReader(RFielName);
		dataReader.read();
		std::vector< std::vector<double>> R = dataReader.getData();
		std::vector<double> eM;
		for(unsigned int i =0; i < R.size();++i)
			eM.push_back(R[i][0]);

		std::vector<QSMatrix<double> > Z;
		for(unsigned int i=0;i< zFileNames.size();i++)
		{
			dataReader=MatrixReader(zFileNames.at(i));
			dataReader.read();
			QSMatrix<double> zM(dataReader.getData());
			Z.push_back(zM);
		}

		std::vector<QSMatrix<double> > G;
		for(unsigned int i = 0; i < gFileNames.size(); i++){
			dataReader=MatrixReader(gFileNames.at(i));
			dataReader.read();
			QSMatrix<double> gM(dataReader.getData());
			G.push_back(gM);
		}
		std::cout << "\nComputing status with REML:\n";

		std::cout << "Total number of entries in Y: " << Y.get_rows() << std::endl;
		std::cout << "Number of fixed effects in X: " << X.get_cols() << std::endl;
		std::cout << "Number of random groups: " << gFileNames.size() << std::endl;
		for(int i = 0; i < gFileNames.size(); i++)
			std::cout << "\tNumber of random variables in group " << i
					<< "(" << zFileNames.at(i) << ")"
					<< ": " << Z.at(i).get_cols() <<std::endl;
		clock_t init,final;
		init=clock();
		Glmm glmm(Y, X, eM, Z, G);
		glmm.setMethod("aireml"); // Use AI-REML iteration algorithm.
		glmm.setMaxIters(20); //Total number of iterations no greater than 20.
		glmm.setPriorEmIters(5); //First 5 iterations are EM algorithm.
		glmm.setDisp(true);//Display the iteration process.
		glmm.init();//Initialize parameters.
      glmm.calc();//Compute variance components, core method.
		glmm.calcFixedEffectPVal();

      QSMatrix<double> tmp1=glmm.calcVarCompSE();
		const std::vector<std::vector<double>> &covM = tmp1.getArrayRef();
		const std::vector<double> & bM 	= glmm.getB();
		const std::vector<double> & pvalM	= glmm.getBPval();
		const std::vector<std::vector<double>> & uM 	= glmm.getU();
		const std::vector<double> & varM 	= glmm.getVarComp();
		std::vector<double> varSEM = glmm.getVarCompSE();
		double logL 	= glmm.getLogL();
		double iters 	= glmm.getIters();
		std::cout << "Writing results into " << outputFileName;

//		nf.setMaximumFractionDigits(5);
		//Now write results from step 4 into a specified file.
		std::unique_ptr<std::ofstream> fout = FileIOHelper::openFileForWrite(outputFileName, false);
		(*fout) << "Number of iteras = " << iters<< std::endl;
		(*fout)<<"logL = " <<logL << std::endl;
		(*fout)<<"BLUEs\n";
		(*fout)<<"Id\tEst.\tPval\n";
		for(int i = 0; i < bM.size(); i++)
			*fout<< "b[" << i << "]\t" << bM[i]<< "\t" << pvalM[i] << "\n";
		(*fout)<<"Variance components\n";
		(*fout)<<"Id\tSource\tVar\tS.E\n";
		for(int i = 0; i < varM.size(); i++)
		{
			(*fout) << "u[" << i << "]\t";
			if(i < varM.size()- 1)
				(*fout) << gFileNames[i] << "\t";
			else
				(*fout) << RFielName << "\t";
			(*fout) << varM[i] << "\t" << varSEM[i] << "\n";
		}
		(*fout)<<"Cov\n";
		for(int i = 0; i < covM.size(); i++)
		{
			for(int j = 0; j < covM[i].size(); j++)
				*fout<<covM[i][j] << "\t";
			*fout<<"\n";
		}
		(*fout) <<"BLUPs\n";
		(*fout) <<"Id\tPred\n";
		for(int i = 0; i < uM.size(); i++)
		{
		   for(int j = 0; j < uM[i].size(); j++)
		   		(*fout) <<"u[" << i << "]["  << j << "]\t" << uM[i][j]<< "\n";
		}
		fout->close();//Finish output and close this file.

		final=clock()-init;
		std::cout << std::endl << (double)final/((double)CLOCKS_PER_SEC) << std::endl;


}
