package LinearModel;
import java.io.PrintWriter;
import java.text.NumberFormat;
import java.util.ArrayList;
import java.util.Date;
import FileIo.FileIOHelper;
import FileIo.MatrixReader;
import Stat.Matrix;
/*
 * This class is used to perform GBLUP analysis.
 * @Author: ZGG
 * @Date: DEC 16, 2014.
 * @Status: Validated on DEC 18, 2014
 */
public class GlmmEngine 
{
	public static void doAnalysis(String YFielName, 
			String XFielName, 
			String RFielName, 
			ArrayList<String> zFileNames, 
			ArrayList<String> gFileNames, 
			String outputFileName) 
	{
		
		System.out.println("\nRunning REML to solve GLMM...");
		MatrixReader dataReader = new MatrixReader(YFielName);
		dataReader.read();
		double yM[][] = dataReader.getData();
		Matrix Y = new Matrix(yM);
		PrintWriter fout1 = FileIOHelper.openFileForWrite("ArrayJava.txt", true);
		NumberFormat nf1 = NumberFormat.getInstance();
		nf1.setMinimumFractionDigits(5);
		nf1.setMaximumFractionDigits(5);
		Y.print(fout1,nf1, 1);
		
		dataReader = new MatrixReader(XFielName);
		dataReader.read();
		double xM[][] = dataReader.getData();
		Matrix X = new Matrix(xM);
		
		dataReader = new MatrixReader(RFielName);
		dataReader.read();
		double rM[][] = dataReader.getData();
		double eM[] = new double[rM.length];
		for(int i = 0; i < eM.length; i++)
			eM[i] = rM[i][0];
		
		ArrayList<Matrix> Z = new ArrayList<Matrix>();
		for(int i = 0; i < zFileNames.size(); i++)
		{
			dataReader = new MatrixReader(zFileNames.get(i));
			dataReader.read();
			double zM[][] = dataReader.getData();
			Z.add(new Matrix(zM));
		}
		
		ArrayList<Matrix> G = new ArrayList<Matrix>();
		for(int i = 0; i < gFileNames.size(); i++)
		{
			dataReader = new MatrixReader(gFileNames.get(i));
			dataReader.read();
			double gM[][] = dataReader.getData();
			G.add(new Matrix(gM));
		}
		System.out.println("Total number of entries in Y: " + Y.getRowDimension());
		System.out.println("Number of fixed effects in X: " + X.getColumnDimension());
		System.out.println("Number of random groups: " + gFileNames.size());
		for(int i = 0; i < gFileNames.size(); i++)
			System.out.println("\tNumber of random variables in group " + i
					+ "(" + zFileNames.get(i) + ")"
					+ ": " + Z.get(i).getColumnDimension());
		System.out.println("\nComputing status with REML:");
		//Run GLMM (mixed model) to estimate variance components.
		long lStartTime = new Date().getTime(); 
		Glmm glmm = new Glmm(Y, X, eM, Z, G);
		glmm.setMethod("aireml"); // Use AI-REML iteration algorithm.
		glmm.setMaxIters(20); //Total number of iterations no greater than 20.
		glmm.setPriorEmIters(5); //First 5 iterations are EM algorithm.
		glmm.setDisp(true);//Display the iteration process.
		glmm.init();//Initialize parameters.
		glmm.calc();//Compute variance components, core method.
		glmm.calcFixedEffectPVal();
		double covM[][] = glmm.calcVarCompSE().getArray();
		double bM[] 	= glmm.getB();
		double pvalM[] 	= glmm.getBPval();
		double uM[][] 	= glmm.getU();
		double varM[] 	= glmm.getVarComp();
		double varSEM[] = glmm.getVarCompSE();
		double logL 	= glmm.getLogL();
		double iters 	= glmm.getIters();
		long lEndTime2 = new Date().getTime(); 
		long diff2 = lEndTime2 - lStartTime; 
		long diffSeconds2 = diff2 / (1000);
		System.out.println("Time spent on GLMM is " + diffSeconds2 + " seconds\n");
				
		System.out.println("Writing results into " + outputFileName);
		NumberFormat nf = NumberFormat.getInstance();
		nf.setMinimumFractionDigits(5);
		nf.setMaximumFractionDigits(5);
		//Now write results from step 4 into a specified file.
		PrintWriter fout = FileIOHelper.openFileForWrite(outputFileName, false);
		fout.println("Number of iteras = " + iters);
		fout.println("logL = " + nf.format(logL));
		fout.println("BLUEs");
		fout.println("Id\tEst.\tPval");
		for(int i = 0; i < bM.length; i++)
			fout.print("b[" + i + "]\t" + nf.format(bM[i]) + "\t" + pvalM[i] + "\n");
		fout.println("Variance components");
		fout.println("Id\tSource\tVar\tS.E");
		for(int i = 0; i < varM.length; i++)
		{
			fout.print("u[" + i + "]\t"); 
			if(i < varM.length - 1)
				fout.print(gFileNames.get(i) + "\t");
			else
				fout.print(RFielName + "\t");
			fout.print(nf.format(varM[i]) + "\t" + nf.format(varSEM[i]) + "\n");
		}
		fout.println("Cov");
		for(int i = 0; i < covM.length; i++)
		{
			for(int j = 0; j < covM[i].length; j++)
				fout.print(covM[i][j] + "\t");
			fout.print("\n");
		}
		fout.println("BLUPs");
		fout.println("Id\tPred");
		for(int i = 0; i < uM.length; i++)
		{
		   for(int j = 0; j < uM[i].length; j++)
		   		fout.print("u[" + i + "]["  + j + "]\t" + nf.format(uM[i][j]) + "\n");
		}
		fout.close();//Finish output and close this file.
	}//End of this method.
}//End of this class.