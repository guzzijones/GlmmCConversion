package Stat;
public class Polynomial extends ConstantsForProb 
{
	protected Polynomial() {}
	public static double p1evl( double x, double coef[], int N) throws ArithmeticException
	{
		double ans;
		ans = x + coef[0];
		for(int i=1; i<N; i++) { ans = ans*x+coef[i]; }
		return ans;
	}
	public static double polevl( double x, double coef[], int N) throws ArithmeticException
	{
		double ans;
		ans = coef[0];
		for(int i=1; i<=N; i++) ans = ans*x+coef[i];
		return ans;
	}
}
