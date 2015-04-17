package Stat;


/**
 * @author u339591
 *
 */
public class PearsonCorrCoef
{
	private double 
		xM[],
		yM[],
		Rxy;
	private int
		N,
		DF;
	
	/*
	 * Constructor.
	 */
	public PearsonCorrCoef(double x[], double y[])
	{
		xM = x;
		yM = y;
		N = xM.length;
		DF = N - 1;
		Rxy = 0;
		calc();
	}
	
	/*
	 * Constructor.
	 */
	public PearsonCorrCoef(double x[][], double y[][])
	{
		int n = x.length;
		xM = new double[n];
		yM = new double[n];
		for(int i = 0; i < n; i++)
			xM[i] = x[i][0];
		for(int i = 0; i < n; i++)
			yM[i] = y[i][0];
		N = xM.length;
		DF = N - 1;
		Rxy = 0;
		calc();
	}
	
	/*
	 * Constructor.
	 */
	public PearsonCorrCoef(float x[], float y[])
	{
		N = x.length;
		xM = new double[N];
		yM = new double[N];
		for(int i = 0; i < N; i++)
		{
			xM[i] = x[i];
			yM[i] = y[i];
		}
		DF = N - 1;
		Rxy = 0;
		calc();
	}
	
	/*
	 * This mehtod is used to calculat pearson corr coef by
	 * Rxy = (1 / DF) * sum(normX * normY);
	 */
	private void calc()
	{
		double[] 
			Zx = calcZscore(xM),
			Zy = calcZscore(yM);
		Rxy = 0;
		for (int i = 0; i < N; i++)
		{
			Rxy += (Zx[i] * Zy[i]); 
		}
		Rxy /= DF;
	}
	/*
	 * This method is used to normalize obs of a sample 
	 * using its mean and std by
	 * ynorm = (y - mu) / std.
	 * @param: y, a (n x 1) array denoting obs of a sample.
	 * @return: ynorm, a (n x 1) array. 
	 */
	private double[] calcZscore(double y[])
	{
		int n = y.length;
		double ynorm[] = new double[n];
		double 
			mu = calcMean(y),
			std = calcSTD(y);
		for(int i = 0; i < n; i++)
		{
			ynorm[i] = (y[i] - mu) / std;
		}
		return ynorm;
	}
	/*
	 * This method is used to calculate standard deviation s
	 * @param: y, a (n x 1) array denoting obs of a sample
	 * @return: std of y.
	 */
	private double calcSTD(double y[])
	{
		double 
			s = 0,
			mu = calcMean(y);
		int 
			n = y.length,
			df = n - 1;
		for (int i = 0; i < n; i++)
		{
			double diff = y[i] - mu;
			s += diff * diff;
		}
		s /= df;
		s = Math.sqrt(s);
		return s;
	}
	/*
	 * This method is used to calculate the mean of a sample.
	 * @param: y, a (n x b) array denoting obs of a sample.
	 * @return: the mean of y.
	 */
	private double calcMean(double y[])
	{
		int n = y.length;
		double mu = 0;
		for (int i = 0; i < n; i++)
		{
			mu+= y[i];
		}
		mu /=n;
		return mu;
	}
	/*
	 * This method is used to access corr coef.
	 */
	public double getCorrCoef()
	{
		return Rxy;
	}
	
	public double calcPValue()
	{
		double t = Rxy / Math.sqrt((1 - (Rxy * Rxy))/ (N - 2));
		double p = Probability.studentT(N - 2, Math.abs(t));
		return p;
	}
	
	
	/*
	 * TEST
	 *
	public static void main(String args[])
	{
		double xM[] = {1, 3 ,5, 6};
		double yM[] = {8, 6 ,4 ,2};
		PearsonCorrCoef cc = new PearsonCorrCoef(xM, yM);
		System.out.println(cc.getCorrCoef());
	}
	*/
}
