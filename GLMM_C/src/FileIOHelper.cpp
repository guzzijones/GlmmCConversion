#include "FileIOHelper.h"
using namespace std;
std::unique_ptr<ifstream> FileIOHelper::openFileForRead(std::string fileName){
   std::unique_ptr<ifstream> outPtr=nullptr;
		try
		{
		 std::unique_ptr<ifstream> outPtrTmp(new ifstream(fileName));
		 outPtr=std::move(outPtrTmp);
		 if(!outPtr->is_open()){
			throw runtime_error("could not open file");
		 }
		}catch(std::exception & e)
		{
			std::cout << "Error while trying to open [" <<  fileName << "] ";
			std::cout << e.what();
			outPtr= nullptr;
			throw;
		}
		return std::move(outPtr);

}


std::unique_ptr<ofstream> FileIOHelper::openFileForWrite(std::string fileName,bool keepwriting){
      std::unique_ptr<ofstream> writePtr=nullptr;
      try{
         if(keepwriting==true){
            std::unique_ptr<ofstream> writePtrTmp(new ofstream(fileName,std::ofstream::out|std::ofstream::app));
            writePtr=std::move(writePtrTmp);
			if(!writePtr->is_open()){
				throw runtime_error("could not open file");
			}
          }
          else{

            std::unique_ptr<ofstream> writePtrTmp(new ofstream(fileName,std::ofstream::out));
            writePtr=std::move(writePtrTmp);
			if(!writePtr->is_open()){
				throw runtime_error("could not open file");
			}
           }
      }catch(std::exception &e){
         std::cout << "Error while trying to write e [" << fileName << "] ";
         std::cout << e.what();
         writePtr=nullptr;
			throw;

      }
      return std::move(writePtr);
}

std::string FileIOHelper::readLine(std::unique_ptr<ifstream> & in){
		std::string str ;
		try
		{
        // std::cout << "reading" << std::endl;
			getline(*in,str);
         //std::cout << "done reading:"<< str << std::endl;

		}catch(std::exception &e)
		{
			std::cout << "Error while trying to read from the file. ";
			throw;
		}
		return str;

}

std::unique_ptr<ofstream> FileIOHelper::writeTextFile(std::string fileName){
      std::unique_ptr<ofstream> writePtr=nullptr;
      try{
            std::unique_ptr<ofstream> writePtrTmp(new ofstream(fileName,std::ofstream::out|std::ofstream::app));
            writePtr=std::move(writePtrTmp);

      }catch(std::exception &e){
         std::cout << "Error while trying to write e [" << fileName << "]";
         std::cout << e.what();
         writePtr=nullptr;
			throw;

      }
      return std::move(writePtr);
}

bool FileIOHelper::writeLine(std::unique_ptr<ofstream> &ofStream, const std::string & str){
	if (str == "") return false;
		try
		{
			(*ofStream) << str;
			return true;
		}catch(std::exception & e)
		{
			std::cout << "Error while trying to write data into the file.";
         std::cout << e.what();
			throw;
			return false;
		}

}
bool FileIOHelper::closeFile(std::unique_ptr<ifstream> &in){

      try{
         in->close();
			return true;
      }catch(std::exception & e){
         std::cout << "Error while closing file" ;
         std::cout << e.what();
			throw;
         return false;
      }
}

bool FileIOHelper::closeFile(std::unique_ptr<ofstream> &in){

      try{
         in->close();
			return true;
      }catch(std::exception & e){
         std::cout << "Error while closing file" ;
         std::cout << e.what();
			throw;
         return false;
      }
}
