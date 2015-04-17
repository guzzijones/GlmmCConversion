package LinearModel;
import java.io.PrintWriter;
import java.text.NumberFormat;
import java.util.ArrayList;

import FileIo.FileIOHelper;
import Stat.JArray;
import Stat.LDLT;
import Stat.Matrix;
import Stat.Probability;
/*
 * @Author: ZHIGANG GUO
 * @Reference: AI_REML paper.
 * @Status: Tested & validated by ZGG at 10/14/2014
 * @Reference: 
 */
public class Glmm 
{
	//Define properties for this class.
	private Matrix 				Y; // n x 1 observation matrix.
	private Matrix 				X; // n x p fixed design matrix, a full rank matrix.
	private double[]			E;
	private ArrayList<Matrix>  	Z; // r array with each of which is a n x nr matrix.
	private ArrayList<Matrix>  	G; // r array with each of which is a nr x nr square matrix.
	private double 				aM[][][]; // r x n x n for ZZ' (random)
	private double 				bM[]; // p x 1 BLUEs, corresponding to X.
	private double 				b_seM[]; // SE for bM.
	private double 				b_pValueM[]; // p x 1 P values of BLUEs
	private double 				uM[][]; // n x r BLUPs
	private double 				logL; // Log Likelihood.
	private double  			dLogL; 
	private int 				r; // Number of random factors.
	private int 				n; // Number of individuals.
	private int					m;
	private boolean 			disp; // Display status of REML?
	private Matrix 				P; 
	private Matrix 				Py; 			
	private Matrix 				H;		
	private Matrix 				R; 
	private Matrix 				s2M;
	private Matrix 				xvxM;
	private double  			s2_seM[];
	private Matrix  			olds2M;
	private Matrix  			tr_PA;
	private String 				method;
	private int 				maxIters;
	private int     			priorEmIterations;
	private int 				iterations;
	private static String		EM = "em",		//EM iteration. 
								AI_REML = "aireml";//AI-REML.
	/*
	 * Constructor to initialize variables. 
	 */
	public Glmm(Matrix yM, // n x 1 2D array. 
			Matrix xM,
			double[] errM,
			ArrayList<Matrix> zM,
			ArrayList<Matrix> gM)
	{
		Y 					= yM; 
		X 					= xM;
		Z 					= zM;
		G 					= gM;
		E 					= errM;
		aM 					= null; 
		r 					= zM.size();
		n 					= Y.getRowDimension();
		m					= X.getColumnDimension(); 	
		bM					= null;
		uM 					= null;
		disp 				= true;
		b_pValueM 			= null;
		b_seM 				= null;
		tr_PA 				= new Matrix(r + 1, 1, 0);
		H 					= new Matrix(r + 1, r + 1);
		R 					= new Matrix(r + 1, 1, 0);
		s2M 				= new Matrix(r + 1, 1, 0);
		olds2M 				= new Matrix(r + 1, 1, 0);
		s2_seM 				= new double[r + 1];
		xvxM				= new Matrix(m, m);
		P 					= null;
		Py 					= new Matrix(n, 1);
		logL 				= -100;
		method 				= EM;
		maxIters 			= 20;
		priorEmIterations 	= 3;
		iterations			= 0;
	}//End of contractor.
	
	/*
	 * Constructor to initialize variables. 
	 */
	public Glmm(Matrix yM, // n x 1 2D array. 
			Matrix xM,
			ArrayList<Matrix> zM,
			ArrayList<Matrix> gM)
	{
		Y 					= yM; 
		X 					= xM;
		Z 					= zM;
		G 					= gM;
		E 					= null;
		aM 					= null; 
		r 					= zM.size();
		m					= X.getColumnDimension(); 					
		n 					= Y.getRowDimension();
		bM					= null;
		uM 					= null;
		disp 				= true;
		b_pValueM 			= null;
		b_seM 				= null;
		tr_PA 				= new Matrix(r + 1, 1, 0);
		H 					= new Matrix(r + 1, r + 1);
		R 					= new Matrix(r + 1, 1, 0);
		s2M 				= new Matrix(r + 1, 1, 0);
		olds2M 				= new Matrix(r + 1, 1, 0);
		s2_seM 				= new double[r + 1];
		xvxM				= new Matrix(m, m);
		P 					= null;
		Py 					= new Matrix(n, 1);
		logL 				= -100;
		method 				= EM;
		maxIters 			= 20;
		priorEmIterations 	= 3;
		iterations          = 0;
	}//End of contractor.
	
	/*
	 * Initialization is critical to affect computational efficient. 
	 * The close prior of variance makes the program quickly go to 
	 * converge.
	 */
	public void init()
	{
		PrintWriter fout3 = FileIOHelper.openFileForWrite("ArrayJavax.txt", true);
		NumberFormat nf3 = NumberFormat.getInstance();
		nf3.setMinimumFractionDigits(0);
		nf3.setMaximumFractionDigits(0);
		X.print(fout3,nf3, 1);
		//Initialize variance components.
		Matrix xt = X.transpose();
		Matrix xtx = xt.times(X);
		Matrix xty = xt.times(Y);
		PrintWriter fout1 = FileIOHelper.openFileForWrite("ArrayJavaxt.txt", true);
		NumberFormat nf1 = NumberFormat.getInstance();
		nf1.setMinimumFractionDigits(5);
		nf1.setMaximumFractionDigits(5);
		xt.print(fout1,nf1, 1);
		
	/*	PrintWriter fout2 = FileIOHelper.openFileForWrite("ArrayJavax.txt", true);
		NumberFormat nf2 = NumberFormat.getInstance();
		nf1.setMinimumFractionDigits(5);
		nf1.setMaximumFractionDigits(5);
		X.print(fout2,nf2, 1);
		*/
		Matrix B = xtx.solve(xty);//Use solve instead of inverting.
		Matrix yfit = Y.minus(X.times(B));
		double pvar = JArray.var(yfit.getArray(), 0);
		double gvar = 0.5 * pvar;
		for(int i = 0; i < r; i++)
		{
			s2M.set(i, 0, gvar / r);
			olds2M.set(i, 0, s2M.get(i, 0));
		}
		s2M.set(r, 0, pvar - gvar);
		olds2M.set(r, 0, s2M.get(r, 0));
		
		//Initialize A.
		aM = new double[r][][];
		for(int i = 0; i < r; i++)
			aM[i] = Z.get(i).times(G.get(i)).times(Z.get(i).transpose()).getArray();
		
		if(E == null)
		{
			E = new double[n];
			for(int i = 0; i < n; i++)
				E[i] = 1.0;
		}
	}//End of this method.
	
	/*
	 * Perform REML analysis.
	 */
	public void calc()
	{
		//Initialize local variables.
		Matrix vM = new Matrix(n, n);
		Matrix vxM = null;
		double logDetV = 0;
		double logDetxVx = 0;
		double oldLogL = -100;
		
		if(disp)
		{
			System.out.print("\tIter\tAlgorithm\tLogL_Diff\tLogL(NumConstrains)\t");
			for(int i = 0; i < r + 1; i++)
				System.out.print("v" + i + "\t");
			System.out.print("\n");
		}
		
		NumberFormat nf = NumberFormat.getInstance();
		nf.setMinimumFractionDigits(4);
		nf.setMaximumFractionDigits(4);
		double crit = 0.0001;
		double y_Ssq = Math.sqrt(Y.norm1())/(n - 1.0);
		double logLikCons = -0.5 * n * Math.log(2 * Math.PI);
		iterations = 0;
		for(int iter = 0; iter < maxIters; iter++)
		{
			this.calcV(vM);// Get V.
			
			//For the solver, please replace this part with your fast solution
			//rather than following LDLT code from this class.
			logDetV = this.calcSymmetrixMatrixInverseByLDLT(vM);//V -> V^-1
							
			// Calculate P matrix.
			vxM = vM.times(X);
			double xvx[][] = xvxM.getArray();
			double x[][] = X.getArray();
			double vx[][] = vxM.getArray();
			for(int i = 0; i < m; i++)
				for(int j = 0; j < m; j++)
				{
					double s = 0;
					for(int k = 0; k < n; k++)
						s += x[k][i] * vx[k][j];
					xvx[i][j] = s;
				}
			//xvxM = X.transpose().times(vxM);
			logDetxVx = this.calcSymmetrixMatrixInverseByLDLT(xvxM);//XV^-1V -> (XV^-1V)^-1
			//xvxM = xvxM.inverse();//The time cost of this step depends on number of fixed effects.
			P = vM.minusInPlace(vxM.times(xvxM).times(vxM.transpose()));//Time consuming.
				
			//Updating variance components.
			if(iter < priorEmIterations)//The first ten rounds are using EM.
				em_reml();//EM>
			else
			{
				if(method.equalsIgnoreCase(AI_REML))
					ai_reml(iter - priorEmIterations);
				else if(method.equalsIgnoreCase(EM))
					em_reml();
			}
			
			int constrain_num = constrain_varcmp(y_Ssq);
				
			// Calculate Log Likelihood.
			this.calcRemlLogLik(logDetV, logDetxVx, logLikCons);
			dLogL = logL - oldLogL;
		
			if(iter == 0 && dLogL < crit)
				dLogL = 100;
			oldLogL = logL;
			for(int i = 0; i <= r; i++)
				olds2M.set(i, 0, s2M.get(i, 0));//Update S2.
			if(disp)
			{
				System.out.print("\t" + (iter + 1) + "\t");
				if(iter < priorEmIterations)
					System.out.print("EM-Reml\t");
				else
					System.out.print("AI-Reml\t");
				System.out.print(nf.format(dLogL) + "\t"  
					+ nf.format(logL) + "(" + constrain_num + ")\t");
				for(int i = 0; i < s2M.getRowDimension(); i++)
					System.out.print(nf.format(s2M.get(i, 0)) + "\t");
				System.out.print("\n");
			}
			if(dLogL < crit)
				break;
			iterations++;
		}//End of each iteration.
			
		// Calculate BLUE of fixed effects.
		//BLUE: B = (X' * V^-1 * X)^-1 * X' * V^-1 * Y. 
		Matrix B = xvxM.times(vxM.transpose().times(Y));
		bM = new double[B.getRowDimension()];
		for(int i = 0; i < bM.length; i++)
			bM[i] = B.get(i, 0);//Save BLUE of b.
		
		// Calculate BLUPs of random effects.
		uM = new double[r][];
		for(int i = 0; i < r; i++)
		{
			//Need further confirmation about u.
			//BLUP of U is U = G * vG * Z'* PY = G * vG * Z' * V-1 * (Y - XB).
			//U = G * vG * Z' * PY.
			Matrix U = G.get(i).times(Z.get(i).transpose()).times(Py).times(s2M.get(i, 0));
			uM[i] = U.getColumnPackedCopy();//Save BLUP of u.
		}
	}//End of this method.
	
	
	public void calcFixedEffectPVal()
	{
		// Calculate SE for b.
		b_seM = new double[bM.length];
		b_pValueM = new double[bM.length];
		for(int i = 0; i < bM.length; i++)
			b_seM[i] = Math.sqrt(xvxM.get(i, i));
		
		//Calculate P value for b.
		for(int i = 0; i < bM.length; i++)
		{
			double t = calcStudentTValue(bM[i], 0, b_seM[i]); 
			b_pValueM[i] = calcPValueUsingTDistribution(t);
		}//End of each explanatory variable.
	}
	
	public Matrix calcVarCompSE()
	{
		//Calculate S.E. for estimates of variance components
		Matrix hM = new Matrix(r + 1, r + 1);//hM, a squared information matrix.
		for(int i = 0; i < r; i++)
		{
			Matrix pa = P.times(new Matrix(aM[i]));
			for(int j = 0; j < r; j++)
			{
				Matrix pb = P.times(new Matrix(aM[j]));
				double tr = fastCalcTrace(pa.getArray(), pb.getArray()); //pa.times(pb).trace();
				hM.set(i, j, tr);
			}
		}//End of loops.
		
		//For diagonal elements of hM.
		Matrix peM = new Matrix(n, n);
		double pe[][] = peM.getArray();
		double p[][] = P.getArray();
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
			Matrix pb = P.times(new Matrix(aM[i]));
			double tr = fastCalcTrace(pe, pb.getArray()); //P.times(pb).trace();
			hM.set(i, r, tr);
			hM.set(r, i, tr);
		}
		
		hM.set(r, r, fastCalcTrace(pe, pe));//.getArray(), P.getArray()));//P.times(P).trace());
		hM = hM.inverse();
		hM.timesEquals(2.0);
		for(int i = 0; i < r + 1; i++)
			s2_seM[i] = Math.sqrt(hM.get(i, i));
		return hM;
	}//End of this method.
	
	
	/*
	 * This method is used to calculate phenotypic variance-covariance matrix.
	 */
	private void calcV(Matrix vM)
	{
		
		double V[][] = vM.getArray();
		//Updating variance-covariance matrix vM.
		for(int i = 0; i < n; i++)
			for(int j = 0; j < n; j++)
				V[i][j] = 0.0;
		
		//V = E.
		double varE = s2M.get(r, 0);
		for(int i = 0; i < n; i++)
			V[i][i] = E[i] * varE;
		
		// V += A[0] * s0 + A[1] * s1 + A[2] * s2 + .... 
		for(int i = 0; i < r; i++)
		{
			double varG = s2M.get(i, 0);
			for(int j = 0; j < n; j++)
				for(int k = 0; k < n; k++)
					V[j][k] += aM[i][j][k] * varG;
		}
	}//End of this method.
	
	/*
	 * This method is used to calculate the inverse of a symmetric matrix.
	 */
	private double calcSymmetrixMatrixInverseByLDLT(Matrix vM)
	{
		LDLT ldlt = new LDLT(vM.getArray());
		ldlt.decompose();
		ldlt.inverse();
		double logDetV = ldlt.getLodDet();//Note now vM was changed to inverse(vM).
		return logDetV;
	}
	
	/*
	 * This method is used to calculate REML log likelihood.
	 */
	private void calcRemlLogLik(double logDetV, 
			double logDetxVx, double logLikCons)
	{
		// YPY = Y' * P * Y;
		double ypy = 0;
		double y[][] = Y.getArray();
		double py[][] = Py.getArray();
		for(int i = 0; i < n; i++)
			ypy += y[i][0] * py[i][0];
		logL = 	-0.5 * logDetV -0.5 * logDetxVx -0.5 * ypy + logLikCons;
	}
	
	/*
	 * This method is used to perform REML analysis 
	 * based on AI-REML method.
	 * In our testing, AI-REML is faster than EM, but
	 * sometimes produces negative components. It has
	 * to been kept in mind to use EM when negative 
	 * variance components are observed.
	 * @Para: iter, index of iteration
	 * @Status: tested and validated by ZHIGANG GUO.
	 */
	private void ai_reml(int iter)
	{
		int i = 0; 
		int j = 0;
		this.calcPY();
		ArrayList<Matrix> APy = new ArrayList<Matrix>();
		for(i = 0; i < r; i++)
			APy.add(new Matrix(aM[i]).times(Py)); 
		
		//R is the first derivative.
		for(i = 0; i < r; i++)
			R.set(i, 0, (Py.transpose().times(APy.get(i))).get(0,0));
		
		Matrix epy = new Matrix(n, 1);
		double epyA[][] = epy.getArray();
		double py[][] = Py.getArray();
		for(i = 0; i < n; i++)
			epyA[i][0] = E[i] * py[i][0];
				
		double PYtEPY = 0;
		for(i = 0; i < n; i++)
			PYtEPY += py[i][0] * epyA[i][0]; 
		R.set(r, 0, PYtEPY);
				
		this.calcPAtrace();
		R = tr_PA.minus(R); 
		R.timesEquals(-0.5);
		
		//H is AI matrix.
		Matrix PxPy = null; 
		for(i = 0; i < r; i++)
		{
			PxPy = P.times(APy.get(i));
			//For diagonal entries of H.
			for(j = 0; j <= i; j++)
			{
				double apy[][] = APy.get(j).getArray();
				double pxpy[][] = PxPy.getArray();
				double YtPAPBPY = 0;
				for(int k = 0; k < n; k++)
					YtPAPBPY += apy[k][0] * pxpy[k][0];
				H.set(j, i, YtPAPBPY);
				H.set(i, j, YtPAPBPY);//For non-diagonal entries of H.
			}//End of each column of H.
		}//End of each row of H.
				
		PxPy = P.times(epy);//E).times(Py);
		for(j = 0; j < r; j++)
		{
			double YtPAPEPY = 0;
			double apy[][] = APy.get(j).getArray();
			double pxpy[][] = PxPy.getArray();
			for(int k = 0; k < n; k++)
				YtPAPEPY += apy[k][0] * pxpy[k][0]; 
			H.set(j, r, YtPAPEPY);
			H.set(r, j, YtPAPEPY);
		}
		
		double YtPEPEPY = 0;
		double pxpy[][] = PxPy.getArray();
		for(i = 0; i < n; i++)
			YtPEPEPY += epyA[i][0] * pxpy[i][0];
		H.set(r, r, YtPEPEPY);
		H.timesEquals(0.5);//this is AI (average information).
		H = H.inverse();//Invert AI matrix.
		Matrix deltaM = H.times(R);
		
		if(dLogL > 1.0 || iter < 3)
			s2M = olds2M.plus(deltaM.timesEquals(0.316)); 
		else
			s2M = olds2M.plus(deltaM);
	}//End of AI-REML method.
	
	private void calcPY()
	{
		//Get array from Matrix in order to 
		//reducing the time for indexing.
		double p[][] = P.getArray();
		double y[][] = Y.getArray();
		double py[][] = Py.getArray();
		
		//Calculate P matrix (n x n) * Y matrix (n x 1).
		double s = 0;
		for(int j = 0; j < n; j++)
		{
			s = 0;
			for(int k = 0; k < n; k++)
				s += p[j][k] * y[k][0];
			py[j][0] = s;
		}//End of each individual.
	}//End of this method.
	
	/*
	 * This method is used to perform REML analysis 
	 * based on EM algorithm or iteration.
	 * In general, more iterations are need for the algorithm to converge,
	 * compared with AI-REML.But EM-REML always give
	 * positive estimation of variance components comparing
	 * to negative variance components with REML or AI-REML
	 * in some cases. This may be considered as an advantage
	 * of EM-REML.
	 * @Author: ZHIGANG GUO.
	 * @Status: the code of this method has been tested at 10/17/2014.
	 */
	private void em_reml()
	{
		// Calculate trace of PA matrix.
		this.calcPAtrace();
		this.calcPY();
					
		//For individual variance component.
		//Avoid Matrix algebra to save memory.
		double s = 0;
		double py[][] = Py.getArray();
		for(int i = 0; i < r; i++)// Note aM[i] is ZGZ'.
		{
			double PytAPy = 0;//PytAPy = y'*P*A*P*y.
			double a[][] = aM[i];
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
			s2_updated /= n;
			s2M.set(i, 0, s2_updated);// Update variance components.
		}//End of each random factor.
	}//End of EM-REML method.
	
	/*
	 * This method is used to calculate t value given its mean and variance.
	 * @Para: estiPara, a point estimate of a parameter.
	 * @Para: expectedPara, expected value of a parameter.
	 * @Para: estimated variance of a parameter.
	 * @Return: t score.
	 */
	private double calcStudentTValue(double estiPara, 
			double expectedPara, double se)
	{
		double 
			// Calculate the difference between the estimated and expected values.
			diff = Math.abs(estiPara - expectedPara), 
			t = diff / se; // Obtain t value by standardizing DIFF by its SE.
		return t;
	}//End of this method.
	
	/*
	 * This method is used to calculate the p value for the 
	 * Hypothesis: MU = 0;
	 * Note that we are performing two-tail test.
	 * @Para: t score.
	 * @Return: p value.
	 */
	private double calcPValueUsingTDistribution(double t)
	{
		return 1 - Probability.studentT(1, t);
	}//End of this method.
	
	/*
	 * This method is used to calculate trace of PA. 
	 * @Author: Zhigang Guo.
	 * @This method has been tested by Zhigang Guo.
	 */
	private void calcPAtrace()
	{
		double[][] p = P.getArray();//n x n matrix.
		double[][] a = null;
		double t = 0;
		for(int i = 0; i < r; i++)
		{
			a = aM[i];//n x n matrix.
			t = fastCalcTrace(p, a);//The method only be used for two matrix with same size (n x n)
			tr_PA.set(i, 0, t);
		}//End of this loop.
		tr_PA.set(r, 0, fastCalcTrace(p, E));
	}//End of this method.
	
	private double fastCalcTrace(double A[][], double B[][])
	{
		double t = 0;
		int n = A.length;
		for(int i = 0; i < n; i++)
			for(int j = 0; j < n; j++)
				t += A[i][j] * B[j][i];
		return t;
	}
	
	private double fastCalcTrace(double A[][], double B[])
	{
		double t = 0;
		int n = A.length;
		for(int i = 0; i < n; i++)
			t += A[i][i] * B[i];
		return t;
	}

	private int constrain_varcmp(double y_Ssq)
	{
	    double delta=0.0;
	    int i = 0, num=0;
	    double constrain[] = new double[r+1];
	    for(i=0; i< r + 1; i++){
	        if(s2M.get(i, 0) < 0)
	        {
	            delta += y_Ssq * 1e-6 - s2M.get(i, 0); 
	            s2M.set(i, 0, y_Ssq * 1e-6);
	            constrain[i] = 1;
	            num++;
	        }
	    }
	    delta /= (r + 1 - num);
	    for(i = 0; i < r + 1; i++)
	    {
	        if(constrain[i] < 1 && s2M.get(i, 0)> delta)
	        {
	        	double s = s2M.get(i, 0) - delta;
	           	s2M.set(i, 0, s);
	        }
	    }
	    return num;
	}//End of this method.
	
	/*
	 * This method is used to set the status of displaying results.
	 * @Para: disp, a boolean variable.
	 */
	public void setDisp(boolean disp)
	{
		this.disp = disp;
	}//End of this method.
	
	/*
	 * This method is used to specify the maximum
	 * number of iterations used for REML iterations.
	 * @Para: iters, maximum number of iterations.
	 */
	public void setMaxIters(int iters)
	{
		maxIters = iters;
	}//End of this method.
	
	/*
	 * This method is used to return BLUEs 
	 * of fixed effects B.
	 * @Return: bM, (p x 1 ) 1D array.
	 */
	public double[] getB()
	{
		return bM;
	}//End of this method.
	
	/*
	 * This method is used to return P values of 
	 * each fixed effect bM.
	 * Note that this p value is calculated based on 
	 * a standard T test (two tails).
	 * @Return: (p x 1) 1D array.
	 */
	public double[] getBPval()
	{
		return b_pValueM;
	}//End of this method.
	
	public double[][] getU()
	{
		return uM;
	}//End of this method.
	
	/*
	 * This method is used to return 
	 * Log of likelihood of REML model.
	 * @Return: logL, logLikelihood of the model.
	 */
	public double getLogL()
	{
		return logL;
	}//End of this method.
	
	/*
	 * This method is used to return 
	 * a (r + 1) 1D array of variance
	 * components. 
	 * The last component is environmental 
	 * variance. 
	 * @Return: (r + 1) 1D array.
	 */
	public double[] getVarComp()
	{
		return s2M.getColumnPackedCopy();
	}//End of this method.
	
	
	public double[] getVarCompSE()
	{
		return s2_seM;
	}
	
	public double getDLog()
	{
		return dLogL;
	}
	
	/*
	 * This method is used to specify the method
	 * used for REML analysis. 
	 * There are three options:
	 * "EM" means solving REML by EM.
	 * "AIREML" means solving REML by AI method.
	 * @Para: method, must be one of above methods.
	 */
	public void setMethod(String method)
	{
		this.method = method;
	}//End of this method.
	
	public void setPriorEmIters(int EmIters)
	{
		priorEmIterations = EmIters;
	}
	
	public int getIters()
	{
		return iterations;
	}
	
}//End of this class.
