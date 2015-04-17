#ifndef FILEIOHELPER_H
#define FILEIOHELPER_H
#include <iostream>
#include <fstream>
#include <exception>
#include <stdexcept>
#include <memory>
/*
 * This class is used to open/close a file and read/write data
 * in the file. With this class, we can conveniently read or
 * write data for a file without more considering about exceptions.
 */
class FileIOHelper{
      public:

      //static methods
   /*
    * This method is used to open a text file.
    * @Para: the name of a text file, a string.
    * @Return: a handle of the file. It will return NULL
    * if opening the file is failed.
    * @Note: This method is a static method. This means that
    * we can use this method as FileIOHelper.openTextFile(R) in
    * other programs.
    */
      static std::unique_ptr<std::ifstream> openFileForRead(std::string filename);
   /*
	 * This method is used to open a text file.
	 * @Para: the name of a text file, a string.
	 * @Return: a handle of the file. It will return NULL
	 * if opening the file is failed.
	 * @Note: This method is a static method. This means that
	 * we can use this method as FileIOHelper.openTextFile(R) in
	 * other programs.
	 */
      static std::unique_ptr<std::ofstream> openFileForWrite(std::string filename,bool keepwriting);
   /*
	 * This method is used to read a line from an opened file.
	 * @Para: a handle of an opened file.
	 * @Return: a string of a line. It will return NULL if reading
	 * is failed.
	 * @Note: This method is a static method. This means that
	 * we can use this method as FileIOHelper.readLine(R) in
	 * other programs.
	 */
      static std::string readLine(std::unique_ptr<std::ifstream> & in);
   /*
	 * This method is used to open a file for writing data into it.
	 * @Para: FILENAME, the name of a file. It is a string.
	 * @Return: a handle of the file. It will return NULL if the
	 * file doesn't exist or cannot be opened successfully.
	 * @Note: This method is a static method. This means that
	 * we can use this method as FileIOHelper.writeTextFile(W) in
	 * other programs.
	 */
      static std::unique_ptr<std::ofstream> writeTextFile(std::string filename);
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

	   static bool writeLine(std::unique_ptr<std::ofstream> &, const std::string & str);
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

	   static bool closeFile(std::unique_ptr<std::ofstream> &);
	   static bool closeFile(std::unique_ptr<std::ifstream> &);


};
#endif
