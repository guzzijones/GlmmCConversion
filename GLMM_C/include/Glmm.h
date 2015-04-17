#ifndef GLMM_H
#define GLMM_H
#define _USE_MATH_DEFINES
#include "QSMatrix.h"
#include <iostream>
#include <string>
#include <vector>
#include "JArray.h"
#include <math.h>
#include "LDLT.h"
#include "solveInterface.h"
#include "FileIOHelper.h"
#include "Probability.h"

class Glmm {
	private:
		void calcV(QSMatrix<double> & vM);
//		const double M_PI = 3.14159265358979323846;
		double calcSymmetrixMatrixInverseByLDLT(QSMatrix<double> & vM);
		void em_reml();
		void calcPAtrace();
		double fastCalcTrace(const std::vector<std::vector<double>> & A, const std::vector<double> & B);
		double fastCalcTrace(const std::vector<std::vector<double>> & A, const std::vector<std::vector<double>> & B);
		void calcPY();
		void ai_reml(int iter);
		int constrain_varcmp(double y_Ssq);
		void calcRemlLogLik(double logDetV, double logDetxVx, double logLikCons);

		QSMatrix<double> Y;
		QSMatrix<double> X;
		std::vector<double> E;
		std::vector<QSMatrix<double>>  	Z; // r array with each of which is a n x nr matrix.
		std::vector<QSMatrix<double>>  	G; // r array with each of which is a nr x nr square matrix.
		std::vector<std::vector<std::vector<double>>> aM; // r x n x n for ZZ' (random)
		std::vector<double>bM; // p x 1 BLUEs, corresponding to X.
		std::vector<double>b_seM; // SE for bM.
		std::vector<double > b_pValueM; // p x 1 P values of BLUEs
		std::vector<std::vector<double>> uM; // n x r BLUPs
		double 				logL; // Log Likelihood.
		double  			dLogL;
		int 				r; // Number of random factors.
		int 				n; // Number of individuals.
		int					m;
		bool 			disp; // Display status of REML?
		QSMatrix<double> 				P;
		QSMatrix<double> 				Py;
		QSMatrix<double> 				H;
		QSMatrix<double> 				R;
		QSMatrix<double> 				s2M;
		QSMatrix<double> 				xvxM;
		std::vector<double>	s2_seM;
		QSMatrix<double>  			olds2M;
		QSMatrix<double>  			tr_PA;
		std::string 				method;
		int 				maxIters;
		int     			priorEmIterations;
		int 				iterations;
		std::string		EM ,	 AI_REML ;//AI-REML.

	protected:


	public:
		Glmm(const QSMatrix<double> & yM,
			const QSMatrix<double> & xM,
			const std::vector<double> & errM,
			const std::vector<QSMatrix<double> > & zM,
			const std::vector<QSMatrix<double> > & gM);
		void init();
		void setMethod(std::string inMethod);
		void setMaxIters(int iters);
		void setPriorEmIters(int EmIters);
		void setDisp(bool disp);
		void calc();
		void calcFixedEffectPVal();
		double calcPValueUsingTDistribution(double t);
   double calcStudentTValue(double estiPara, double expectedPara, double se);
   QSMatrix<double> calcVarCompSE();

   const std::vector<double> & getB();
   const std::vector<double> & getBPval();
   const std::vector<std::vector<double>> & getU();
   std::vector<double>  getVarComp();
   const std::vector<double> & getVarCompSE();
   double getDLog();
   double getLogL();
   int getIters();

};

#endif
