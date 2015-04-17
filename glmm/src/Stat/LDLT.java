package Stat;


/*
 * This class is used to perform LDLT decomposition for matrix A and 
 * is also used to calculate inversion of matrix A based on LDLT.
 * @Author: ZHIGANG GUO,
 * @Date: OCT-06-2014
 * @Status: the algorithm and code tested and validated on OCT-08-2014
 */
public class LDLT
{
	private double A[][];
	private double[] D;
	private int n; 
	/*
	 * @Para: aM, a symmetric and positive definite square matrix.
	 */
	public LDLT(double aM[][])
	{
		A = aM; //A is an n x n square matrix.
		n = A.length;//The order of matrix A.
		D = null;  //Prepare this array for D.
	}//End of this method.
	
	/*
	 * This method is used to perform LDLT decomposition
	 * for a symmetric and positive definite square matrix A.
	 * A = L * D * L'
	 * where L is lower triangle matrix with diagonal entry of ones.
	 * D is a diagonal matrix;
	 * L' is a transpose of square matrix L.
	 * After decomposition, L will be saved in lower triangle of A.
	 * The upper triangle of A is still the original A.
	 * The reason of doing this is to save memory. 
	 */
	public void decompose()
	{
		D = new double[n];//Initialize D.
		double v[] = new double[n];
		for(int j = 0; j < n; j++)
		{
			if(j > 0)
			{
				for(int k = 0; k <= j - 1; k++)
					v[k] = A[j][k] * D[k];
				
				double sum = 0;
				for(int k = 0; k <= j - 1; k++)
					sum += A[j][k] * v[k];
				v[j] = A[j][j] - sum;
				D[j] = v[j];
				
				if(j < n - 1)
					for(int k = j + 1; k < n; k++)
					{
						sum = 0;
						for(int w = 0; w <= j - 1; w++)
							sum += A[k][w] * v[w];
						A[k][j] = (A[j][k] - sum) / v[j];
					}
			}
			else
			{
				//Initialize...checked.
				v[0] = A[0][0];
				D[0] = v[0];
				for(int k = 1; k < n; k++)
					A[k][0] /= v[0];
			}
		}//End of this loop
		for(int i = 0; i < n; i++)
			A[i][i] = 1.0;
	}//End of this method.
	
	/*
	 * This method is used to compute the inverse of matrix A
	 * based on its LDLT decomposition.
	 */
	public double[][] inverse()
	{
		//Copy lower triangle data to upper triangle of matrix A.
		//So, the lower triangle part is L
		//The upper triangle part is L' which will be used later.
		//The diagonal of matrix A is 1.
		for (int j = 0; j < n; j++) 
			for (int i = j + 1; i < n; i++) 
			    A[j][i] = A[i][j];
		
		//Solve AZ = B with B is an n x n identify matrix.
		//This is equivalent to invert matrix A, but 
		//by using solving we save computational time.
		double b[] = new double[n];
		for(int i = 0; i < n; i++)
		{
			for(int j = 0; j < n; j++)
				b[j] = (j == i)?1.0:0.0;
			TriangleMatrix.solveLForLDLT(A, b, n);//Update b.
			//Note here we only update the lower triangle of A.
			for(int j = i; j < n; j++)
				A[j][i] = b[j];
		}
		
		//Solve DY = Z;
		//This is the fastest way of solving this equation given 
		//D is diagonal matrix.
	   for (int i = 0; i < n; i++) 
	      for(int j = 0; j <= i; j++)
	         A[i][j] /= D[i];//Only update the lower triangle of matrix A.
	  
	   //Solve L'X = Y with L' saved as upper triangle of A matrix.
	   	for(int i = 0; i < n; i++)
		{
		    for(int j = 0; j < n; j++)
			    b[j] = (j < i)? 0.0: A[j][i];
		   	TriangleMatrix.solveUForLDLT(A, b, n);//Update b.
			for(int j = i; j < n; j++)
				A[j][i] = b[j];//Update upper triangle of matrix A.
		}
	  
		//Make A a symmetric matrix.
		for (int j = 0; j < n; j++) 
		   for (int i = j + 1; i < n; i++) 
		     A[j][i] = A[i][j];
	   return A;
	}//End of this method.

	/*
	 * This method is used to get D after LDLT decomposition.
	 */
	public double[] getD()
	{
		return D;
	}
	
	/*
	 * This method is used to calculate Log(DET(A)).
	 */
	public double getLodDet()
	{
		double logD = 0; //Math.log(d);
	    for (int j = 0; j < n; j++)
	      if(D[j] > 0)
	    	  logD += Math.log(D[j]);
	    return logD;
	}
}//End of this class.
