package FileIo;
import java.io.BufferedReader;
import java.io.BufferedWriter;
import java.io.FileReader;
import java.io.FileWriter;
import java.io.IOException;
import java.io.PrintWriter;
/*
 * This class is used to open/close a file and read/write data
 * in the file. With this class, we can conveniently read or 
 * write data for a file without more considering about exceptions.
 */
public class FileIOHelper
{
	/*
	 * This method is used to open a text file.
	 * @Para: the name of a text file, a string.
	 * @Return: a handle of the file. It will return NULL
	 * if opening the file is failed.
	 * @Note: This method is a static method. This means that
	 * we can use this method as FileIOHelper.openTextFile(R) in 
	 * other programs.
	 */
	public static BufferedReader openFileForRead(String fileName)
	{
		BufferedReader reader = null;
		try
		{
			reader = new BufferedReader(new FileReader(fileName));
		}catch(IOException ex)
		{
			System.out.println("Error while trying to open [" + fileName + "]");
			System.out.println(ex.toString());
			reader = null;
		}
		return reader;
	}//End of OPENTEXTFILE method.
	
	/*
	 * This method is used to open a text file.
	 * @Para: the name of a text file, a string.
	 * @Return: a handle of the file. It will return NULL
	 * if opening the file is failed.
	 * @Note: This method is a static method. This means that
	 * we can use this method as FileIOHelper.openTextFile(R) in 
	 * other programs.
	 */
	public static PrintWriter openFileForWrite(String fileName, boolean keepwriting)
	{
		PrintWriter writer = null;
		try
		{
			FileWriter fout = new FileWriter(fileName, keepwriting);
            writer = new PrintWriter(fout, keepwriting);
		}
		catch(IOException ex)
		{
			System.out.println("Error while trying to write [" + fileName + "]");
			System.out.println(ex.toString());
			writer = null;
		}
		return writer;
	}//End of OPENTEXTFILE method.
	
	
	/*
	 * This method is used to read a line from an opened file.
	 * @Para: a handle of an opened file.
	 * @Return: a string of a line. It will return NULL if reading
	 * is failed.
	 * @Note: This method is a static method. This means that
	 * we can use this method as FileIOHelper.readLine(R) in 
	 * other programs.
	 */
	public static String readLine(BufferedReader reader)
	{
		String str = null;
		try
		{
			str = reader.readLine();
			
		}catch(IOException ex)
		{
			System.out.println("Error while trying to read from the file.");
			str = null;
		}
		return str;
	}//End of READLINE method.
	/*
	 * This method is used to open a file for writing data into it.
	 * @Para: FILENAME, the name of a file. It is a string.
	 * @Return: a handle of the file. It will return NULL if the 
	 * file doesn't exist or cannot be opened successfully.
	 * @Note: This method is a static method. This means that
	 * we can use this method as FileIOHelper.writeTextFile(W) in 
	 * other programs.
	 */
	public static PrintWriter writeTextFile(String fileName)
	{
		PrintWriter writer = null; // Default NULL class.
		try
		{
			writer = new PrintWriter(new FileWriter(fileName));
		}catch(IOException ex)
		{
			System.out.println("Error while trying to write data into [" + fileName + "]");
			System.out.println(ex.toString());
			writer = null;
		}
		return writer;
	}//End of WRITETEXTFILE method.
	/*
	 * This method is used to write a line in a specified file.
	 * @Para: a handle of a specified file.
	 * @Para: STR, a string needed to write.
	 * @Return: true if STR is written into the file,
	 * otherwise, it will return false.
	 * @Note: This method is a static method. This means that
	 * we can use this method as FileIOHelper.writeLine(W, STR) in 
	 * other programs.
	 */
	public static boolean writeLine(BufferedWriter writer, String str)
	{
		if (str == null) return false;
		try
		{
			writer.write(str);
			return true;
		}catch(IOException ex)
		{
			System.out.println("Error while trying to write data into the file.");
			return false;
		}
	}//End of WRITELINE method.
	
	/*
	 * This method is to test if a file specified by READER
	 * is closed successfully.
	 * @Para: a handle of a file opener.
	 * @Return: true if the specified file is closed successfully, 
	 * otherwise, it returns false.
	 * @Note: This method is a static method. This means that
	 * we can use this method as FileIOHelper.closeFile(R) in 
	 * other programs.
	 */
	public static boolean closeFile(BufferedReader reader)
	{
		if (reader == null) return false;
		try
		{
			reader.close();
			return true;
		}catch(IOException ex)
		{
			System.out.println("Error while trying to close the file.");
			return false;
		}
	}//End of CLOSEFILE method.
}//End of this class.
