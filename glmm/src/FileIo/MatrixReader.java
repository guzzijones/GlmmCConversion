package FileIo;
import java.io.BufferedReader;
import java.util.ArrayList;
public class MatrixReader
{
	private String fileName;
	private double X[][];
	
	public MatrixReader(String fName)
	{
		fileName = fName;
		X = null;
	}//End of this method.
	
	public double[][] getData()
	{
		return X;
	}
	/*
	 * This method is used to write a Matrix into a specified file.
	 */
	public void read()
	{
		String strSplitter = ",";
		BufferedReader fin 	= FileIOHelper.openFileForRead(fileName);
		String sentence 	= null;
		String tokens[] 	= null;
		ArrayList<double[]> dataM = new ArrayList<double[]>();
		// Read data from the input file.
		while ((sentence = FileIOHelper.readLine(fin)) != null) 
		{
			tokens = sentence.trim().split(strSplitter);
			int n = tokens.length;
			double rowM[] = new double[n];
			for(int i = 0; i < n; i++)
				rowM[i] = Double.valueOf(tokens[i]);
			dataM.add(rowM);
		}
		FileIOHelper.closeFile(fin);
		int n = dataM.size();
		int m = dataM.get(0).length;
		X = new double[n][m];
		for(int i = 0; i < n; i++)
			X[i] = dataM.get(i);
	}//End of this method.
}//End of this class.
