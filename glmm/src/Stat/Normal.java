package Stat;

import java.util.Random;

public class Normal
{
	/*
	 * Generates a randomly normal distributed number.
	 * @param: mu, the mean of the normal distribution.
	 * @param: sd, the standard deviation of the normal distribution.
	 * @return: the random number sampled from the normal distribution.
	 * @Author: ZG.
	 */
	public static float next(double mu, double sd)
	{
		Random generator = new Random();
		return (float)(mu + generator.nextGaussian() * sd);
	}
	
	public static double pdf(double x, double mu, double sigma2)
	{
		double prob = 0;
		prob = Math.exp(-0.5 * Math.pow(x - mu, 2) / sigma2 )/ Math.sqrt(2 * Math.PI * sigma2);
		return prob;
	}
	
	public static double pdf(double x)
	{
		return pdf(x, 0, 1); 
	}
}
