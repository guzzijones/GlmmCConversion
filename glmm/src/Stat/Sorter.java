package Stat;

public class Sorter 
{
	private double[] itsDataM;
	private int[] itsRankM;
	private int n;
	
	public Sorter(double dataM[])
	{
		itsDataM = dataM;
		n = itsDataM.length;
		itsRankM = null;
	}
	
	public double[] get()
	{
		return itsDataM;
	}
	
	public int[] getRank()
	{
		return itsRankM;
	}
	
	public void calc()
	{
		double rawM[] = new double[n];
		for(int i = 0; i < n; i++)
			rawM[i] = itsDataM[i];
		
		int j = 0;
		double tmp = 0;
		for(int i=0;i<n;i++)
		{
		   j = i;
		   for(int k = i;k<n;k++)
		   {
		     if(itsDataM[j] > itsDataM[k])
		     {
		       j = k;
		     }
		   }
		   tmp = itsDataM[i];
		   itsDataM[i] = itsDataM[j];
		   itsDataM[j] = tmp;
		}
		
		int tempRankM[] = new int[n];
		int maxR = tempRankM[0];
		for(int i = 1; i < n; i++)
			if(itsDataM[i] > itsDataM[i - 1])
			{
				tempRankM[i] = tempRankM[i - 1] + 1;
				if(maxR < tempRankM[i])
					maxR = tempRankM[i];
			}
			else
				tempRankM[i] = tempRankM[i - 1];
				
		itsRankM = new int[n];
		for(int i = 0; i < n; i++)
		{
			int r = 0;
			for(int k = 0; k < n; k++)
			{
				if(rawM[i] == itsDataM[k])
				{
					r = tempRankM[k];
					break;
				}
			}
			itsRankM[i] = maxR - r;
		}
	}
	
	/*
	 * For testing purposes.
	 */
	public static void main(String args[])
	{
		
		double aM[] = {3};
		double bM[] = new double[aM.length];
		for(int i = 0; i < aM.length; i++)
			bM[i] = aM[i];
		Sorter sort = new Sorter(aM);
		sort.calc();
		//double baM = sort.get();
		int rM[] = sort.getRank();
		for(int i = 0; i < aM.length; i++)
			System.out.println(bM[i] + "\t" + rM[i]);
	}
}
