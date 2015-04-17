package Stat;

public class MinCalculator 
{
	private float 
		itsYM[],
		itsMinVal;
	private int 
		itsMinIndex;
	public MinCalculator(float yM[])
	{
		itsYM = yM;
		itsMinVal = 0;
		itsMinIndex = 0;
	}
	
	public void calc()
	{
		itsMinIndex = 0;
		itsMinVal = itsYM[itsMinIndex];
		int n = itsYM.length;
		if (n > 1)
		{
			for (int i = 1; i < itsYM.length; i++)
			{
				if (itsYM[i] < itsMinVal)
				{
					itsMinVal = itsYM[i];
					itsMinIndex = i;
				}
			}
		}
	}
	
	public float getMinValue()
	{
		return itsMinVal;
	}
	
	public int getMinIndex()
	{
		return itsMinIndex;
	}

}
