package Stat;

import java.util.Random;

public class NormRand
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
}
