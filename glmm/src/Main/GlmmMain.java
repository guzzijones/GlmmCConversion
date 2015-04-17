package Main;
import java.util.ArrayList;
import LinearModel.GlmmEngine;
/*
 * This program is developed to perform GLMM.  
 * @Author: ZGG.
 * @Date: March 17, 2015.
 * @Status: Validated by ZGG, March 16, 2015.
 */
public class GlmmMain
{
	public static void main(String args[])
	{
		System.out.println("*******************************");
		System.out.println("General linear mixed model tool");
		System.out.println("Developed by Zhigang Guo, Syngenta. ");
		System.out.println("Version: 1.0");
		System.out.println("Date: March 19, 2015");
		System.out.println("*******************************");
		
		System.out.println("\nUsage: ");
		System.out.print("java -jar glmm.jar -y yFileName.csv\n"); 
		System.out.print("\t\t   -f xFileName.csv\n");
		System.out.print("\t\t   -r zFileName1.csv:zFileName2.csv:zFileName3.csv:...:zFileNamep.csv\n");
		System.out.print("\t\t   -g gFileName1.csv:gFileName2.csv:gFileName3.csv:...:gFileNamep.csv\n");
		System.out.print("\t\t   -e rFileName.csv\n");
		System.out.print("\t\t   -o outputFileName.txt\n");
			
		String YFielName = null;
		String XFielName = null;
		String RFielName = null;
		String zFileNameStr = null;
		String gFileNameStr = null;
		String outputFileName = null;
		String splitter = ":";
		for(int i = 0; i < args.length; i++)
		{
			if(args[i].equalsIgnoreCase("-y"))
				YFielName = args[i + 1];
			else if(args[i].equalsIgnoreCase("-f"))
				XFielName = args[i + 1];
			else if(args[i].equalsIgnoreCase("-e"))
				RFielName = args[i + 1];
			else if(args[i].equalsIgnoreCase("-r"))
				zFileNameStr = args[i + 1];
			else if(args[i].equalsIgnoreCase("-g"))
				gFileNameStr = args[i + 1];
			else if(args[i].equalsIgnoreCase("-o"))
				outputFileName = args[i + 1];
		}//End of this loop.
		
		ArrayList<String> zFileNames = new ArrayList<String>();
		ArrayList<String> gFileNames = new ArrayList<String>();
		String termsZ[] = zFileNameStr.split(splitter);
		String termsG[] = gFileNameStr.split(splitter);
		int q = termsZ.length;
		for(int i = 0; i < q; i++)
		{
			zFileNames.add(termsZ[i]);
			gFileNames.add(termsG[i]);
		}
		System.out.println("\nYour inputs:");
		System.out.println("Y from the input file " +  YFielName);
		System.out.println("X (fixed design matrix) from the input file " +  XFielName);
		System.out.print("Z (design matrix) from " + q + " files: ");
		for(int i = 0; i < q; i++)
		{
			System.out.print(zFileNames.get(i));
			if(i < q - 1)
			 System.out.print(",");
		}
		System.out.println();
		System.out.print("G (covariance structure matrix) from " + q + " files: ");
		for(int i = 0; i < q; i++)
		{
			System.out.print(gFileNames.get(i));
			if(i < q - 1)
				System.out.print(",");
		}
		System.out.println();
		System.out.println("E (residual design matrix) from " + RFielName);
		System.out.println("Results will be output into " + outputFileName);
				
		GlmmEngine.doAnalysis(YFielName, XFielName, RFielName, 
				zFileNames, gFileNames, outputFileName);
		System.out.println("Results were written into the file " + outputFileName);
		System.out.println("The GLMM analysis is complete. ");
	}//End of this method.
}//End of this class.