#include "Glmm.h"
Glmm::Glmm(const QSMatrix<double> & yM, // n x 1 2D array.
			const QSMatrix<double> & xM,
			const std::vector<double> & errM,
		const std::vector< QSMatrix<double> > & zM,
			const std::vector< QSMatrix<double> > & gM)
			{
				EM = "em",		//EM iteration.
				AI_REML = "aireml";//AI-REML.


				Y 					= yM;
				X 					= xM;
				Z 					= zM;
				G 					= gM;
				E 					= errM;
				aM 					;
				r 					= zM.size();
				n 					= Y.get_rows();
				m					= X.get_cols();
				bM					;
				uM 					;
				disp 				= true;
				b_pValueM 			;
				b_seM 				;

				tr_PA 				= QSMatrix<double>(r + 1, 1, 0);
				H 					= QSMatrix<double>(r + 1, r + 1,0);
				R 					= QSMatrix<double>(r + 1, 1, 0);
				s2M 				= QSMatrix<double>(r + 1, 1, 0);
				olds2M 				= QSMatrix<double>(r + 1, 1, 0);
				std::vector<double>sTmp;
				sTmp.resize(r+1);
				s2_seM 				= sTmp;
				xvxM				= QSMatrix<double>(m, m,0);
				P 					;
				Py 					= QSMatrix<double>(n, 1,0);
				logL 				= -100;
				method 				= EM;
				maxIters 			= 20;
				priorEmIterations 	= 3;
				iterations			= 0;

}//End of contractor.


void Glmm::init()
{
		//Initialize variance components.
//std::cout << "\nhere";
  // std::unique_ptr<std::ofstream> outx1= FileIOHelper::openFileForWrite("Testx.txt",false);
   //   FileIOHelper::writeLine(outx1,X.print()) ;
	//	std::cout << "x rows: " << X.get_rows()<< std::endl;
	//	std::cout << "x cols: " << X.get_cols()<< std::endl;


		QSMatrix<double> xt = X.transpose();
  //  std::unique_ptr<std::ofstream> outxt= FileIOHelper::openFileForWrite("xt.txt",false);
    //  FileIOHelper::writeLine(outxt,xt.print(1,5)) ;
      //std::unique_ptr<std::ofstream> outY= FileIOHelper::openFileForWrite("Y.txt",false);
      //FileIOHelper::writeLine(outY,Y.print(1,5));
	//	std::cout << "xt rows: " << xt.get_rows()<< std::endl;
	//	std::cout << "xt cols: " << xt.get_cols()<< std::endl;


   //   std::cout << "doing mult X" << std::endl;
	//	std::cout << "x rows: " << X.get_rows()<< std::endl;
	//	std::cout << "x cols: " << X.get_cols()<< std::endl;
		QSMatrix<double> xtx = xt*X;
    //  std::cout << "doing mult Y" << std::endl;
		QSMatrix<double> xty = xt*Y;
    //  std::cout << "done" << std::endl;
//        std::cout << xty.print(1,5);
        //xt and Y are the same

//         std::unique_ptr<std::ofstream> xtxout= FileIOHelper::openFileForWrite("xtx.txt",false);
 //       FileIOHelper::writeLine(xtxout,xtx.print(1,5)) ;

		solveInterface solvxtx(xtx) ;

		QSMatrix<double> B = solvxtx.solve(xty);//Use solve instead of inverting.
		//B is all zeros., should be one column of numbers
		 //solve is broken
	//	 std::unique_ptr<std::ofstream> B1= FileIOHelper::openFileForWrite("B.txt",false);
		// std::cout << "B" <<std::endl;
      //  std::cout << B.print(1,20) ;
		QSMatrix<double> yfit = Y-(X*B);
	//	std::unique_ptr<std::ofstream> yfito= FileIOHelper::openFileForWrite("yfit.txt",false);
     //   FileIOHelper::writeLine(yfito,yfit.print(1,10)) ;
		double pvar = JArray::var(yfit.getArrayRef(), 0);
	//	std::cout << "pvar" << pvar << std::endl;
		double gvar = 0.5 * pvar;
      //  std::cout <<"s2m before for: "<< std::endl<< s2M.print(1,20);
		for(int i = 0; i < r; i++)
		{
			s2M.set(i, 0, gvar / r);
			olds2M.set(i, 0, s2M.get(i, 0));
		}
//		std::cout << "s2m after:" << std::endl;
		s2M.set(r, 0, pvar - gvar);
		olds2M.set(r, 0, s2M.get(r, 0));
  //      std::unique_ptr<std::ofstream> outvm1= FileIOHelper::openFileForWrite("s2M.txt",false);
    //        FileIOHelper::writeLine(outvm1,s2M.print(1,20)) ;
      //      std::cout << s2M.print(1,20);
		//Initialize A.
		//aM;
		for(int i = 0; i < r; i++){
	/*	QSMatrix<double>  Zati=Z.at(i);
		QSMatrix<double>  Gati=G.at(i);
        QSMatrix<double>  temp;
        QSMatrix<double>  temp2;
       temp=Z.at(i)*G.at(i);
    temp2=Z.at(i).transpose();
        temp=temp*temp2;
       aM.push_back(temp.getArrayRef());*/
       aM.reserve(r);
		aM.push_back((Z.at(i)*
                (G.at(i))
                *(Z.at(i).transpose())).getArrayRef());
		}
		if(E.empty() )
		{
			E.clear() ;
			for(int i = 0; i < n; i++){
                E.reserve(n);
				E.push_back( 1.0);
				}
		}
}//End of this method.
void Glmm::setMethod(std::string inMethod){
	this->method=inMethod;
}

void Glmm::setMaxIters(int iters)
{
		maxIters = iters;
}//End of this method.
void Glmm::setPriorEmIters(int EmIters)
	{
		priorEmIterations = EmIters;
	}

void Glmm::setDisp(bool disp)
	{
		this->disp = disp;
	}//End of this method.

void Glmm::calc()
	{
		//Initialize local variables.
		QSMatrix<double> vM(n, n,0);
		QSMatrix<double> vxM ;
		double logDetV = 0;
		double logDetxVx = 0;
		double oldLogL = -100;

		if(disp)
		{
			std::cout << "\tIter\tAlgorithm\tLogL_Diff\tLogL(NumConstrains)\t";
			for(int i = 0; i < r + 1; i++)
				std::cout << "v" << i  << "\t";
			std::cout << "\n";
		}

//		NumberFormat nf = NumberFormat.getInstance();
//		nf.setMinimumFractionDigits(4);
//		nf.setMaximumFractionDigits(4);
		double crit = 0.0001;
		double y_Ssq = sqrt(Y.norm1())/(n - 1.0);
		double logLikCons = -0.5 * n * log(2 * M_PI);
		iterations = 0;
		for(int iter = 0; iter < maxIters; iter++)
		{
        //    std::unique_ptr<std::ofstream> outvm1= FileIOHelper::openFileForWrite("vMbeforeCaclv.txt",false);
          //  FileIOHelper::writeLine(outvm1,vM.print(1,20)) ;
			calcV(vM);// Get V.
            //std::unique_ptr<std::ofstream> outvm2= FileIOHelper::openFileForWrite("vMafterCaclv.txt",false);
           // FileIOHelper::writeLine(outvm2,vM.print(1,20)) ;
			//For the solver, please replace this part with your fast solution
			//rather than following LDLT code from this class.
			//vm is the same, problem in calcsysmm
			//need to check the value of vM and x
			//error is here - calcsymmetrixinverseby ldt changes vm
	//				  std::unique_ptr<std::ofstream> outvm= FileIOHelper::openFileForWrite("vMbefore.txt",false);
      //      FileIOHelper::writeLine(outvm,vM.print(1,20)) ;
			logDetV = calcSymmetrixMatrixInverseByLDLT(vM);//V -> V^-1
            //should be -4957.44
			// Calculate P matrix.

//                std::unique_ptr<std::ofstream> vxMf= FileIOHelper::openFileForWrite("vxM.txt",false);
 //            FileIOHelper::writeLine(vxMf,vxM.print(1,10)) ;
              std::unique_ptr<std::ofstream> vxMf2= FileIOHelper::openFileForWrite("X.txt",false);
  //           FileIOHelper::writeLine(vxMf2,X.print(1,10)) ;
			vxM = vM*X;

			//xvx is all zeros
			std::vector<std::vector<double>> & xvx= xvxM.getArrayRefEdit();
			//std::unique_ptr<std::ofstream> xvxf= FileIOHelper::openFileForWrite("xvx.txt",false);
			//QSMatrix<double> xvxM(xvx);
            //FileIOHelper::writeLine(xvxf,xvxM.print(1,10)) ;

			const std::vector<std::vector<double>> & x = X.getArrayRef();
			//std::unique_ptr<std::ofstream> xf= FileIOHelper::openFileForWrite("x.txt",false);
			//QSMatrix<double> xM(x);
            //FileIOHelper::writeLine(xf,xM.print(1,10)) ;

			const std::vector<std::vector<double>> & vx = vxM.getArrayRef();
            //std::unique_ptr<std::ofstream> vxf= FileIOHelper::openFileForWrite("vx.txt",false);
            //QSMatrix<double> vxM(vx);
            //FileIOHelper::writeLine(vxf,vxM.print(1,10)) ;
			//print each of the above and compare.
			//there is a problem here.
			//need to check the value of vm and x.  they must be diff.
			//or it is
			//vxm is incorrect.

			for(int i = 0; i < m; i++)
				for(int j = 0; j < m; j++)
				{
					double s = 0;
					for(int k = 0; k < n; k++)
						s += (x[k][i] * (double)vx[k][j]);
					xvxM(i,j) = s;
				}
			//xvxM = X.transpose().times(vxM);
//            std::unique_ptr<std::ofstream> xvxMf= FileIOHelper::openFileForWrite("xvxM.txt",false);
 //            FileIOHelper::writeLine(xvxMf,xvxM.print(1,10)) ;
			logDetxVx = this->calcSymmetrixMatrixInverseByLDLT(xvxM);//XV^-1V -> (XV^-1V)^-1
			//should be 507.44
			//xvxM = xvxM.inverse();//The time cost of this step depends on number of fixed effects.
			QSMatrix<double> tmpV=vxM*xvxM*(vxM.transpose());
			P = vM.minusInPlace(tmpV);//Time consuming.

			//Updating variance components.
			if(iter < priorEmIterations)//The first ten rounds are using EM.
				em_reml();//EM>
			else
			{
				if(method==AI_REML)
					ai_reml(iter - priorEmIterations);
				else if(method==EM)
					em_reml();
			}

			int constrain_num = constrain_varcmp(y_Ssq);

			// Calculate Log Likelihood.
			calcRemlLogLik(logDetV, logDetxVx, logLikCons);
			dLogL = logL - oldLogL;

			if(iter == 0 && dLogL < crit)
				dLogL = 100;
			oldLogL = logL;
			for(int i = 0; i <= r; i++)
				olds2M.set(i, 0, s2M.get(i, 0));//Update S2.
			if(disp)
			{
				std::cout << "\t" << (iter + 1) << "\t";
				if(iter < priorEmIterations)
					std::cout << "EM-Reml\t";
				else
					std::cout << "AI-Reml\t";
				std::cout << dLogL << "\t"
					<< logL << "(" << constrain_num << ")\t";
				for(int i = 0; i < s2M.get_rows(); i++)
					std::cout << s2M.get(i, 0) << "\t";
				std::cout << "\n";
			}
			if(dLogL < crit)
				break;
			iterations++;
		}//End of each iteration.

		// Calculate BLUE of fixed effects.
		//BLUE: B = (X' * V^-1 * X)^-1 * X' * V^-1 * Y.
		QSMatrix<double> B = xvxM*(vxM.transpose()*(Y));

		for(int i = 0; i < B.get_rows(); i++)
			bM.push_back(B.get(i, 0));//Save BLUE of b.

		// Calculate BLUPs of random effects.
		uM.resize(r);
		for(int i = 0; i < r; i++)
		{
			//Need further confirmation about u.
			//BLUP of U is U = G * vG * Z'* PY = G * vG * Z' * V-1 * (Y - XB).
			//U = G * vG * Z' * PY.
			QSMatrix<double> U = G.at(i)*(Z.at(i).transpose())*(Py)*(s2M.get(i, 0));
			uM[i] = U.getColumnPackedCopy();//Save BLUP of u.
		}
	}//End of this method.



void Glmm::calcV(QSMatrix<double> & vM)
	{

		std::vector<std::vector<double>> & V = vM.getArrayRefEdit();
		//Updating variance-covariance matrix vM.
		for(int i = 0; i < n; i++)
			for(int j = 0; j < n; j++)
				V[i][j] = 0.0;

		//V = E.
		double varE = s2M.get(r, 0);

		for(int i = 0; i < n; i++)
			V[i][i] = E[i] * varE;
   //     double test = aM[0][0][0];
		// V += A[0] * s0 + A[1] * s1 + A[2] * s2 + ....
		for(int i = 0; i < r; i++)
		{
			double varG = s2M.get(i, 0);
			for(int j = 0; j < n; j++)
				for(int k = 0; k < n; k++)
					V[j][k] += aM[i][j][k] * varG;
		}
	}//End of this method.

double Glmm::calcSymmetrixMatrixInverseByLDLT(QSMatrix<double> & vM)
	{
//			  std::unique_ptr<std::ofstream> outvm= FileIOHelper::openFileForWrite("vMinsidecalcsym.txt",false);
 //           FileIOHelper::writeLine(outvm,vM.print(1,5)) ;
		LDLT ldlt(vM.getArrayRefEdit());
//		 std::unique_ptr<std::ofstream> outvm2= FileIOHelper::openFileForWrite("vMafterconstructor.txt",false);
 //           FileIOHelper::writeLine(outvm2,vM.print(1,5)) ;
		ldlt.decompose();
//			 std::unique_ptr<std::ofstream> outvm3= FileIOHelper::openFileForWrite("vMafterdecompose.txt",false);
 //           FileIOHelper::writeLine(outvm3,vM.print(1,5)) ;
		ldlt.inverse();
//			 std::unique_ptr<std::ofstream> outvm4= FileIOHelper::openFileForWrite("vMafterinverse.txt",false);
 //           FileIOHelper::writeLine(outvm4,vM.print(1,20)) ;
		double logDetV = ldlt.getLodDet();//Note now vM was changed to inverse(vM).
		return logDetV;
	}

void Glmm::calcPAtrace()
	{
		const std::vector<std::vector<double>> & p = P.getArrayRef();//n x n matrix.
		std::vector<std::vector<double>>  a ;
		double t = 0;
		for(int i = 0; i < r; i++)
		{
			a = aM[i];//n x n matrix.
			t = fastCalcTrace(p, a);//The method only be used for two matrix with same size (n x n)
			tr_PA.set(i, 0, t);
		}//End of this loop.
		tr_PA.set(r, 0, fastCalcTrace(p, E));
	}//End of this method.

double Glmm::fastCalcTrace(const std::vector<std::vector<double>> & A, const std::vector<double> & B)
	{
		double t = 0;
		int n = A.size();
		for(int i = 0; i < n; i++)
			t += A[i][i] * B[i];
		return t;
	}

double Glmm::fastCalcTrace(const std::vector<std::vector<double>> & A, const std::vector<std::vector<double>> & B)
	{
		double t = 0;
		int n = A.size();
		for(int i = 0; i < n; i++)
			for(int j = 0; j < n; j++)
				t += A[i][j] * B[j][i];
		return t;
	}

void Glmm::calcPY()
	{
		//Get array from Matrix in order to
		//reducing the time for indexing.
	//	double p[][] = P.getArrayRef();
		//double y[][] = Y.getArrayRef();
	//	double py[][] = Py.getArrayRef();

		//Calculate P matrix (n x n) * Y matrix (n x 1).
		double s = 0;
		for(int j = 0; j < n; j++)
		{
			s = 0;
			for(int k = 0; k < n; k++)
				s += P.getArrayRef()[j][k] * Y.getArrayRef()[k][0];
			Py.getArrayRefEdit()[j][0] = s;
		}//End of each individual.
	}//End of this method.

void Glmm::em_reml()
	{
		// Calculate trace of PA matrix.
		this->calcPAtrace();
        std::unique_ptr<std::ofstream> tr_PAf= FileIOHelper::openFileForWrite("tr_PAf.txt",false);
        FileIOHelper::writeLine(tr_PAf,tr_PA.print(1,10)) ;
		this->calcPY();
       std::unique_ptr<std::ofstream> tr_PYf= FileIOHelper::openFileForWrite("Py.txt",false);
        FileIOHelper::writeLine(tr_PYf,Py.print(1,10)) ;
		//For individual variance component.

		//Avoid Matrix algebra to save memory.
		double s = 0;
		const std::vector<std::vector<double>> & py = Py.getArrayRef();
		for(int i = 0; i < r; i++)// Note aM[i] is ZGZ'.
		{
			double PytAPy = 0;//PytAPy = y'*P*A*P*y.
			const  std::vector<std::vector<double>> & a = aM[i];
			for(int j = 0; j < n; j++)
			{
				s = 0;
				for(int k = 0; k < n; k++)
					s += a[j][k] * py[k][0];
				s *= py[j][0];
				PytAPy += s;
			}//End of each individual.
			R.set(i, 0, PytAPy);
		}//End of each variance component.

		//For residual variance.
		//Avoid Matrix algebra to save memory.
		double YtPEPY = 0;//YtPEPY = y'*P*E*P*y.
		for(int i = 0; i < n; i++)
		{
			s = py[i][0];
			YtPEPY += s * s * E[i];
		}
		R.set(r, 0, YtPEPY);

		// Update variance components.
		double s2_updated = 0;
		for(int i = 0; i < r + 1; i++)
		{
			double s2_old = olds2M.get(i, 0);
			double s4_old = s2_old * s2_old;
			s2_updated = s2_old * n
					+ s4_old * R.get(i, 0)
					- s4_old * tr_PA.get(i, 0);
			s2_updated /= (double)n;
			s2M.set(i, 0, s2_updated);// Update variance components.
		}//End of each random factor.
	}//End of EM-REML method.

void Glmm::ai_reml(int iter)
	{
		int i = 0;
		int j = 0;
		calcPY();
        std::unique_ptr<std::ofstream> tr_PYf= FileIOHelper::openFileForWrite("Py.txt",false);
        FileIOHelper::writeLine(tr_PYf,Py.print(1,10)) ;
		std::vector<QSMatrix<double>> APy;
		for(i = 0; i < r; i++)
			APy.push_back(QSMatrix<double>(aM[i])*(Py));

		//R is the first derivative.
		for(i = 0; i < r; i++)
			R.set(i, 0, (Py.transpose()*(APy.at(i)))(0,0));

		QSMatrix<double> epy(n, 1,0);
		std::vector<std::vector<double>> & epyA = epy.getArrayRefEdit();
		const std::vector<std::vector<double>> & py = Py.getArrayRef();
		for(i = 0; i < n; i++)
			epyA[i][0] = E[i] * py[i][0];

		double PYtEPY = 0;
		for(i = 0; i < n; i++)
			PYtEPY += py[i][0] * epyA[i][0];
		R.set(r, 0, PYtEPY);

		calcPAtrace();
		R = tr_PA-R;
		R=R*(-0.5);

		//H is AI matrix.
		QSMatrix<double> PxPy ;
		for(i = 0; i < r; i++)
		{
			PxPy = P*(APy.at(i));
			//For diagonal entries of H.
			for(j = 0; j <= i; j++)
			{
				const std::vector<std::vector<double>> & apy = APy.at(j).getArrayRef();
				const std::vector<std::vector<double>> & pxpy = PxPy.getArrayRef();
				double YtPAPBPY = 0;
				for(int k = 0; k < n; k++)
					YtPAPBPY += apy[k][0] * pxpy[k][0]; //this number is correct
				H.set(j, i, YtPAPBPY);
				H.set(i, j, YtPAPBPY);//For non-diagonal entries of H.
			}//End of each column of H.
		}//End of each row of H.

		PxPy = P*(epy);//E).times(Py);
		for(j = 0; j < r; j++)
		{
			double YtPAPEPY = 0;
			const std::vector<std::vector<double>> & apy = APy.at(j).getArrayRef();
			const std::vector<std::vector<double>> & pxpy = PxPy.getArrayRef();
			for(int k = 0; k < n; k++)
				YtPAPEPY += apy[k][0] * pxpy[k][0];
			H.set(j, r, YtPAPEPY);
			H.set(r, j, YtPAPEPY);
		}
      //  std::unique_ptr<std::ofstream> hf1= FileIOHelper::openFileForWrite("h1.txt",false);
       // std::cout << "H1:";
       // std::cout << H.print(1,10) << std::endl;
		double YtPEPEPY = 0;
		const std::vector<std::vector<double>> &pxpy = PxPy.getArrayRef();
		for(i = 0; i < n; i++)
			YtPEPEPY += epyA[i][0] * pxpy[i][0];
		H.set(r, r, YtPEPEPY);
       // std::cout << "H2:";
        //std::cout << H.print(1,10) << std::endl;
		H=H*(0.5);//this is AI (average information).
		solveInterface Hsolve(H);
       // std::cout << "H3:";
        //std::cout << H.print(1,10) << std::endl;
		H = Hsolve.inverse();//Invert AI matrix. //blank after inverse.
		QSMatrix<double> deltaM = H*(R);
      //  std::unique_ptr<std::ofstream> hf= FileIOHelper::openFileForWrite("h.txt",false);
   //    std::cout << "H:";
     //   std::cout << H.print(1,10) << std::endl;
      //  FileIOHelper::writeLine(hf,H.print(1,10)) ;
		if(dLogL > 1.0 || iter < 3)
			s2M = olds2M+(deltaM*(0.316));
		else
			s2M = olds2M+(deltaM);
	}//End of AI-REML method.

int Glmm::constrain_varcmp(double y_Ssq)
	{
	    double delta=0.0;
	    int i = 0, num=0;
	    std::vector<double> constraine;
		constraine.resize(n+1);
	    for(i=0; i< r + 1; i++){
	        if(s2M.get(i, 0) < 0)
	        {
	            delta += y_Ssq * 1e-6 - s2M.get(i, 0);
	            s2M.set(i, 0, y_Ssq * 1e-6);
	            constraine[i] = 1;
	            num++;
	        }
	    }
	    delta /= (double)(r + 1 - num);
	    for(i = 0; i < r + 1; i++)
	    {
	        if(constraine[i] < 1 && s2M.get(i, 0)> delta)
	        {
	        	double s = s2M.get(i, 0) - delta;
	           	s2M.set(i, 0, s);
	        }
	    }
	    return num;
	}//End of this method.

void Glmm::calcRemlLogLik(double logDetV, double logDetxVx, double logLikCons)
	{
		// YPY = Y' * P * Y;
		double ypy = 0;

		const std::vector<std::vector<double>> & y = Y.getArrayRef();

		const std::vector<std::vector<double>> & py = Py.getArrayRef();

		for(int i = 0; i < n; i++)
			ypy += y[i][0] * py[i][0];
		logL = 	-0.5 * logDetV -0.5 * logDetxVx -0.5 * ypy + logLikCons;
	}

void Glmm::calcFixedEffectPVal()
	{
		// Calculate SE for b.
		b_seM.resize(bM.size());
		b_pValueM.resize(bM.size());
		for(int i = 0; i < bM.size(); i++)
			b_seM[i] = std::sqrt(xvxM.get(i, i));

		//Calculate P value for b.
		for(int i = 0; i < bM.size(); i++)
		{
			double t = calcStudentTValue(bM[i], 0, b_seM[i]);
			b_pValueM[i] = calcPValueUsingTDistribution(t);
		}//End of each explanatory variable.
	}

double Glmm::calcPValueUsingTDistribution(double t)
	{
		return 1 - Probability::studentT(1, t);
	}

double Glmm::calcStudentTValue(double estiPara, double expectedPara, double se)
   {
      double
         // Calculate the difference between the estimated and expected values.
         diff = std::abs(estiPara - expectedPara),
         t = diff / (double)se; // Obtain t value by standardizing DIFF by its SE.
      return t;
   }

QSMatrix<double> Glmm::calcVarCompSE()
	{
		//Calculate S.E. for estimates of variance components
		QSMatrix<double> hM(r + 1, r + 1,0);//hM, a squared information matrix.
		for(int i = 0; i < r; i++)
		{
         QSMatrix<double> tmp(aM[i]);
			QSMatrix<double> pa = P*(tmp);
			for(int j = 0; j < r; j++)
			{
            QSMatrix<double> tmp2(aM[j]);
				QSMatrix<double> pb = P*(tmp2);
				double tr = fastCalcTrace(pa.getArrayRef(), pb.getArrayRef()); //pa.times(pb).trace();
				hM.set(i, j, tr);
			}
		}//End of loops.

		//For diagonal elements of hM.
		QSMatrix<double> peM(n, n,0);
		std::vector<std::vector<double>> & pe = peM.getArrayRefEdit();
		const std::vector<std::vector<double>> & p = P.getArrayRef();
		for(int i = 0; i < n; i++)
		{
			for(int j = 0; j < n; j++)
			{
				double s = 0;
				s = p[i][j] * E[j];
				pe[i][j] = s;
			}
		}

		for(int i = 0; i < r; i++)
		{
         QSMatrix<double> tmp(aM[i]);
			QSMatrix<double> pb = P*(tmp);
			double tr = fastCalcTrace(pe, pb.getArrayRef()); //P.times(pb).trace();
			hM.set(i, r, tr);
			hM.set(r, i, tr);
		}

		hM.set(r, r, fastCalcTrace(pe, pe));//.getArray(), P.getArray()));//P.times(P).trace());
      solveInterface tmpSolve(hM);
		hM = tmpSolve.inverse();
		hM=hM*(2.0);
		for(int i = 0; i < r + 1; i++)
			s2_seM[i] = std::sqrt(hM.get(i, i));
		return hM;
	}//End of this method.


const std::vector<double> & Glmm::getB()
	{
		return bM;
	}//End of this method.

const std::vector<double> & Glmm::getBPval()
	{
		return b_pValueM;
	}//End of this method.

const std::vector<std::vector<double>> & Glmm::getU()
	{
		return uM;
	}//End of this method.
std::vector<double>  Glmm::getVarComp()
	{
		return s2M.getColumnPackedCopy();
	}//End of this method.

const std::vector<double> & Glmm::getVarCompSE()
	{
		return s2_seM;
	}

double Glmm::getDLog()
	{
		return dLogL;
	}
double Glmm::getLogL()
	{
		return logL;
	}//End of this method.
int Glmm::getIters()
	{
		return iterations;
	}

