package Stat;
/*
 * This class is used to calculate the likelihood of yM 
 * given the linear model
 * yM = xM * bM + eM.
 * where yM is observed phenotype;
 * xM are predictors. 
 * bM is the estimated predictor effects.
 * eM is normal residual error with mean zero
 * and variance sigma2.
 */
public class LikelihoodCalculator
{
	private double
		itsYM[],
		itsXM[][],
		itsMu,
		itsSigma2,
		itsBM[],
		itsLogLik,
		itsLik;
	private int
		itsN,
		itsM;
	public LikelihoodCalculator(
			double yM[], 
			double xM[][], 
			double bM[],
			double mu,
			double sigma2)
	{
		itsYM = yM;
		itsXM = xM;
		itsMu = mu;
		itsSigma2 = sigma2;
		itsBM = bM;
		itsN = itsYM.length;
		itsM = itsXM[0].length;
		itsLogLik = 0.0;
		itsLik = 0.0;
		this.calc();
	}
	/*
	 * This method will calculate the sum of log likelihood of YM
	 * given mu, bM, and sigma2.
	 */
	private void calc()
	{
		double sumResidualSquare = 0;
		for(int i = 0; i < itsN; i++)
		{
			// First calculate FITTED yM using predictor xM and bM.
			double fittedY = itsMu;
			for(int j = 0; j < itsM; j++)
				fittedY += itsXM[i][j] * itsBM[j];
			// Then calculate the RESIDUAL of each individual.
			double residual = itsYM[i] - fittedY;
			// Finally sum over the square of RESIDUAL.
			sumResidualSquare += Math.pow(residual, 2);
		}
		// Calculate the likelihood using NORMAL PDF.
		double cons = - 0.5 * (double)itsN * Math.log(2 * Math.PI * itsSigma2);
		itsLogLik = cons - sumResidualSquare / ( 2 * itsSigma2);
		itsLik = Math.exp(itsLogLik);
	}
	
	public double getLogLikelihood()
	{
		return itsLogLik;
	}
	public double getLikelihood()
	{
		return itsLik;
	}
}
