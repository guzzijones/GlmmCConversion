package Stat;

import java.util.ArrayList;

/**
 * @Use: this class provides the calculation of basic statistics based on 1-D array.
 * @author Zhigang Guo.
 */
public final class JArray
{
	
	/*
	 * Calculate mean of a sample.
	 */
	public static double mean(double yM[][], int columnIndex)
	{
	    double 
	    	mu = 0;
	    int 
	    	n = yM.length;
		for (int i = 0; i < n; i++)
			mu += yM[i][columnIndex];
		mu /= n;
		return mu;
	}
	
	/*
	 * Calculate mean of a sample.
	 */
	public static double mean(double yM[])
	{
	    double 
	    	mu = 0;
	    int 
	    	n = yM.length;
		for (int i = 0; i < n; i++)
			mu += yM[i];
		mu /= n;
		return mu;
	}
	
	
	/*
	 * Calculate varaince of a sample.
	 */
	public static double var(double yM[][],  int columnIndex)
	{
		double
			mu = mean(yM, columnIndex),
			var = 0,
			temp = 0;
		int 
			n = yM.length,
			df = n;
		for (int i = 0; i < n; i++)
		{
			temp = yM[i][columnIndex] - mu;
			var += (temp * temp);
		}//End Of All the Individuals.
		var /= df;
		return var;
	}
	
	/*
	 * Calculate varaince of a sample.
	 */
	public static double var(double yM[])
	{
		double
			mu = mean(yM),
			var = 0,
			temp = 0;
		int 
			n = yM.length,
			df = n - 1;
		for (int i = 0; i < n; i++)
		{
			temp = yM[i]- mu;
			var += (temp * temp);
		}//End Of All the Individuals.
		var /= df;
		return var;
	}
	
	/*
	 * Calculate std of a sample.
	 */
	public static double std(double yM[][], int columnIndex)
	{
		double v = var(yM, columnIndex);
		return (double)Math.sqrt(v);
	}
	
	/*
	 * Calculate std of a sample.
	 */
	public static double std(double yM[])
	{
		double v = var(yM);
		return (double)Math.sqrt(v);
	}
	
	
	/*
	 * Calculate the minimum value of an array
	 */
	public static double min(double yM[])
	{
		double minVal = yM[0];
		for (int i = 0; i < yM.length; i++)
		{
			if (yM[i] < minVal)
			{
				minVal = yM[i];
			}
		}
		return minVal;
	}
	
	/*
	 * Calculate the max value of an array
	 */
	public static int max(int yM[])
	{
		int maxVal = yM[0];
		for (int i = 0; i < yM.length; i++)
		{
			if (yM[i] > maxVal)
			{
				maxVal = yM[i];
			}
		}
		return maxVal;
	}
	
	public static double[][] convert(ArrayList<Double> y)
	{
		int n = y.size();
		double yM[][] = new double[n][1];
		for(int i = 0; i < n; i++)
			yM[i][0] = y.get(i);
		return yM;
	}
	
	/*
	 * This method is used to extract an unique list (STRINGN)
	 * for a given list.
	 * @Author: ZGG
	 * @Status: Validated on March 2, 2015.
	 */
	public static ArrayList<String> getUniq(ArrayList<String> nameList)
	{
		//Define an unique list.
		ArrayList<String> uniqList = new ArrayList<String>();
		for(int i = 0; i < nameList.size(); i++)
			if(!uniqList.contains(nameList.get(i)))
				uniqList.add(nameList.get(i));
		return uniqList;//return this uniqList.
	}//End of this method.
	
	public static void demo(double yM[][])
	{
		for(int i = 0; i < yM.length; i++)
		{
			for(int j = 0; j < yM[i].length; j++)
				System.out.print(yM[i][j] + "\t");
			System.out.println();
		}
	}
	
	public static void demo(ArrayList<String> aList)
	{
		for(int i = 0; i < aList.size(); i++)
			System.out.println(aList.get(i));
	}
}
