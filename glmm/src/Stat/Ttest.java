package Stat;
/*
 * This class is used to provide a t-test to test the significance of 
 * the difference of two sample means from a population with equal 
 * variances.
 * @Status: this program has been tested and validated. Any change of code
 * should be cautious.
 */
public class Ttest 
{
	/*
	 * This method is used to test the significance of means
	 * from two populations xM and yM with equal variance.
	 * This is a two-tail test.
	 * @Para: xM[], 1D array.
	 * @Para: yM[], 2D array.
	 */
	public static double calc(double xM[], 
			double yM[], 
			String testType)
	{
		int s = xM.length;
		double zM[] = new double[s];
		for(int i = 0; i < s; i++)
			zM[i] = xM[i] - yM[i];
		double nx = zM.length;
		double mux = JArray.mean(zM);
		double sd = JArray.std(zM);
		double se = sd / Math.sqrt(nx);
		double t = mux/se;
		double pval = 1;
		if(testType.equalsIgnoreCase("both"))//Two tailed test.
			pval = 2 * Probability.studentT(nx - 1, -1 * Math.abs(t));
		else if(testType.equalsIgnoreCase("right"))//Right one tailed test.
			pval = Probability.studentT(nx - 1, -1 * t);
		else if(testType.equalsIgnoreCase("left"))//Left one tailed test.
			pval = Probability.studentT(nx - 1, t);
		return pval;
	}//End of this method.
	
	public static void main(String args[])
	{
		double p1[] = {30.02, 29.99,30.11, 29.97, 30.01, 29.99};
		double p2[] = {29.89, 29.93, 29.72, 29.98, 30.02, 29.98};
		double pval = Ttest.calc(p1, p2, "both");
		System.out.println("pval = " + pval);
	}

}
