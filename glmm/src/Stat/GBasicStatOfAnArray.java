/**
 * 
 */
package Stat;

/**
 * @Use: this class provides the calculation of basic statistics based on 1-D array.
 * @author Zhigang Guo.
 */
public class GBasicStatOfAnArray
{
	// Member variables.
	private double 
		itsDataM[],
		itsMin,
		itsMax,
		itsMean,
		itsVar,
		itsSTD;
	private int
		itsMinIndex,
		itsMaxIndex,
		itsSampSize;
	private boolean
		isSampSizeCalculated,
		isMeanCalculated,
		isVarCalculated;
	// Constructor.
	public GBasicStatOfAnArray(double dataM[])
	{
		itsDataM = dataM;
		itsSampSize = 0;
		itsMinIndex = 0;
		itsMaxIndex = 0;
		itsMin = itsDataM[itsMinIndex];
		itsMax = itsDataM[itsMaxIndex];
		isSampSizeCalculated = false;
		isMeanCalculated = false;
		isVarCalculated = false;
	}
	
	public void sampleSize()
	{
		itsSampSize = itsDataM.length;
		isSampSizeCalculated = true;
	}

	// Member methods.
	public void min()
	{
		if (!isSampSizeCalculated)
		{
			sampleSize();
		}
		itsMinIndex = 0;
		itsMin = itsDataM[itsMinIndex];
		double 
			candidateMin;
		for(int i = 1; i < itsSampSize; i++)
		{
			candidateMin = itsDataM[i];
			if (candidateMin < itsMin) // Can we find the smaller one?
			{
				itsMinIndex = i;
				itsMin = itsDataM[itsMinIndex];
			}
		}//End Of all the individuals in the sample.
	}
	public void max()
	{
		if (!isSampSizeCalculated)
		{
			sampleSize();
		}
		itsMaxIndex = 0;
		itsMax = itsDataM[itsMaxIndex];
		double 
			candidateMin;
		for(int i = 1; i < itsSampSize; i++)
		{
			candidateMin = itsDataM[i];
			if (candidateMin > itsMax) // Can we find the bigger one?
			{
				itsMaxIndex = i;
				itsMax = itsDataM[itsMaxIndex];
			}
		}//End Of all the individuals in the sample.
	}
	public void mean()
	{
		if (!isSampSizeCalculated)
		{
			sampleSize();
		}
		double s = 0;
		for (int i = 0; i < itsSampSize; i++)
		{
			s += itsDataM[i];
		}
		s /= itsSampSize;
		itsMean = s;
		isMeanCalculated = true;
	}
	public void var()
	{
		if(!isMeanCalculated)
		{
			mean();
		}
		double
		temp,
		s = 0,
		df = itsSampSize - 1; // df: degree of freedom.
		for (int i = 0; i < itsSampSize; i++)
		{
			temp = (itsDataM[i] - itsMean);
			s += (temp * temp);
		}//End Of All the Individuals.
		s /= df;
		itsMean = s;
		isVarCalculated = true;
	}
	public void std()
	{
		if(!isVarCalculated)
		{
			var();
		}
		itsSTD = Math.sqrt(itsVar);
	}
	
	public double getSampSize()
	{
		return itsSampSize;
	}
	public double getMin()
	{
		return itsMin;
	}
	public int getMinIndex()
	{
		return itsMinIndex;
	}
	public double getMax()
	{
		return itsMax;
	}
	public int getMaxIndex()
	{
		return itsMaxIndex;
	}
	public double getMean()
	{
		return itsMean;
	}
	public double getVar()
	{
		return itsVar;
	}
	public double getStd()
	{
		return itsSTD;
	}
}
