package Stat;

public class TruncatedNormalDistribution
{
	public static double calcUpperMean(double x, double mu, double variance)
	{
		double sigma = Math.sqrt(variance);
		double alpha = (x - mu) / sigma;
		double lamda_alpha = Normal.pdf(alpha) / (1 - Probability.normal(alpha));
		return mu + sigma * lamda_alpha;
	}
	
	public static double calcLowerMean(double x, double mu, double variance)
	{
		double sigma = Math.sqrt(variance);
		double alpha = (x - mu) / sigma;
		double lamda_alpha = -1 * Normal.pdf(alpha) / Probability.normal(alpha);
		return mu + sigma * lamda_alpha;
	}

}
