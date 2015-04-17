package Stat;

public class TriangleMatrix 
{
	 /**
     * <p>
     * Solves for non-singular lower triangular matrices using forward substitution.
     * <br>
     * b = L<sup>-1</sup>b<br>
     * <br>
     * where b is a vector, L is an n by n matrix.<br>
     * </p>
     *
     * @param L An n by n non-singular lower triangular matrix. Not modified.
     * @param b A vector of length n. Modified.
     * @param n The size of the matrices.
     */
    public static void solveL(double L[][] , double []b , int n )
    {
    	double sum = 0;
    	for( int i = 0; i < n; i++ )
    	{
    	  sum = b[i];
    	  for( int k = 0; k < i; k++)
    	  	  sum -= L[i][k] * b[k];
    	  b[i] = sum / L[i][i];
    	}
    }//End of this method.
    
	 /**
     * <p>
     * Solves for non-singular lower triangular matrices using forward substitution.
     * <br>
     * b = L<sup>-1</sup>b<br>
     * <br>
     * where b is a vector, L is an n by n matrix.<br>
     * </p>
     *
     * @param L An n by n non-singular lower triangular matrix. Not modified.
     * @param b A vector of length n. Modified.
     * @param n The size of the matrices.
     */
    public static void solveLForLDLT(double L[][] , double []b , int n )
    {
    	double sum = 0;
    	for( int i = 0; i < n; i++ )
    	{
    	  sum = b[i];
    	  for( int k = 0; k < i; k++)
    	  	  sum -= L[i][k] * b[k];
    	  b[i] = sum;
    	}
    }//End of this method.
    
    /**
     * <p>
     * Solves for non-singular lower triangular matrices using forward substitution.
     * <br>
     * b = L<sup>-1</sup>b<br>
     * <br>
     * where b is a vector, L is an n by n matrix.<br>
     * </p>
     *
     * @param L An n by n non-singular lower triangular matrix. Not modified.
     * @param b A vector of length n. Modified.
     * @param n The size of the matrices.
     */
    public static void solveLForLDLT(float L[][] , float []b , int n )
    {
    	float sum = 0;
    	for( int i = 0; i < n; i++ )
    	{
    	  sum = b[i];
    	  for( int k = 0; k < i; k++)
    	  	  sum -= L[i][k] * b[k];
    	  b[i] = sum;
    	}
    }//End of this method.
    
     /**
     * <p>
     * This is a forward substitution solver for non-singular lower triangular matrices.
     * <br>
     * b = (L<sup>T</sup>)<sup>-1</sup>b<br>
     * <br>
     * where b is a vector, L is an n by n matrix.<br>
     * </p>
     * <p>
     * L is a lower triangular matrix, but it comes up with a solution as if it was
     * an upper triangular matrix that was computed by transposing L.
     * </p>
     *
     * @param L An n by n non-singular lower triangular matrix. Not modified.
     * @param b A vector of length n. Modified.
     * @param n The size of the matrices.
     */
    public static void solveTranL(double L[][], double []b, int n)
    {
        for( int i = n-1; i >= 0; i--)
        {
            double sum = b[i];
            for( int k = i + 1; k < n; k++)
                 sum -= L[k][i]* b[k];
            b[i] = sum/L[i][i];
        }
    }//end of this method.
    
    /**
     * <p>
     * This is a forward substitution solver for non-singular upper triangular matrices.
     * <br>
     * b = U<sup>-1</sup>b<br>
     * <br>
     * where b is a vector, U is an n by n matrix.<br>
     * </p>
     *
     * @param U An n by n non-singular upper triangular matrix. Not modified.
     * @param b A vector of length n. Modified.
     * @param n The size of the matrices.
     */
    public static void solveU( double U[][], double []b , int n )
    {
        for( int i = n - 1; i >= 0; i-- )
        {
            double sum = b[i];
            for( int j = i+1; j <n; j++ )
                 sum -= U[i][j] * b[j];
            b[i] = sum/U[i][i];
        }
    }
    
    /**
     * <p>
     * This is a forward substitution solver for non-singular upper triangular matrices.
     * <br>
     * b = U<sup>-1</sup>b<br>
     * <br>
     * where b is a vector, U is an n by n matrix.<br>
     * </p>
     *
     * @param U An n by n non-singular upper triangular matrix. Not modified.
     * @param b A vector of length n. Modified.
     * @param n The size of the matrices.
     */
    public static void solveUForLDLT( double U[][], double []b , int n )
    {
        for( int i = n - 1; i >= 0; i-- )
        {
            double sum = b[i];
            for( int j = i+1; j <n; j++ )
                 sum -= U[i][j] * b[j];
            b[i] = sum;
        }
    }
    
    /**
     * <p>
     * This is a forward substitution solver for non-singular upper triangular matrices.
     * <br>
     * b = U<sup>-1</sup>b<br>
     * <br>
     * where b is a vector, U is an n by n matrix.<br>
     * </p>
     *
     * @param U An n by n non-singular upper triangular matrix. Not modified.
     * @param b A vector of length n. Modified.
     * @param n The size of the matrices.
     */
    public static void solveUForLDLT( float U[][], float []b , int n )
    {
        for( int i = n - 1; i >= 0; i-- )
        {
            float sum = b[i];
            for( int j = i+1; j <n; j++ )
                 sum -= U[i][j] * b[j];
            b[i] = sum;
        }
    }
}