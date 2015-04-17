#include <unistd.h>
#include "GlmmEngine.h"
#include "FileIOHelper.h"
#include "MatrixReader.h"
//#include "QSMatrix.h"
using namespace std;
int main(int argc, char * argv[]){
//////////////////////////////////
/// testing
//
////////////////////////////////
   //test open file for read
//   std::unique_ptr<ifstream> openFile=FileIOHelper::openFileForRead("./src/driver.cpp");   
   //test open file for write
 //  std::unique_ptr<ofstream> outputfile=FileIOHelper::openFileForWrite("t.txt",false);
  // FileIOHelper::writeLine(outputfile,"This is a new line");
  // *outputfile << "this is a message";
  // FileIOHelper::closeFile(outputfile);

   //test readline 
  // std::cout << "test readline: " << FileIOHelper::readLine(openFile) << std::endl;

  // openFile->close();

//	MatrixReader FileInput("g1.c");
//	FileInput.read();

/*
	MatrixReader readMatrix("./testdata/z1.csv");
	readMatrix.read();
	std::vector<std::vector<double>> myVectorDouble=readMatrix.getData();
	QSMatrix<double> MyQSMatrix(myVectorDouble);
	std::cout << "rows: " << MyQSMatrix.get_rows() << std::endl;
	std::cout << "columns: " << MyQSMatrix.get_cols() << std::endl;

*/

//////////////////////
/// end testing
//
///////////////////////



 if((argc <=1)||(argv[1])==NULL)
   {//there is nothing
		std::cerr<< "*******************************";
		std::cerr <<"General linear mixed model tool" << std::endl;
		std::cerr <<"Developed by Zhigang Guo, Syngenta. "<< std::endl;
		std::cerr <<"Version: 1.0"<< std::endl;
		std::cerr <<"Date: March 19, 2015"<< std::endl;
		std::cerr <<"*******************************"<< std::endl;
		
		std::cerr << "\nUsage: ";
		std::cerr <<"./driver.exe -y yFileName.csv\n"; 
		std::cerr <<"\t\t   -f xFileName.csv\n";
		std::cerr <<"\t\t   -r zFileName1.csv:zFileName2.csv:zFileName3.csv:...:zFileNamep.csv\n";
		std::cerr <<"\t\t   -g gFileName1.csv:gFileName2.csv:gFileName3.csv:...:gFileNamep.csv\n";
		std::cerr <<"\t\t   -e rFileName.csv\n";
		std::cerr <<"\t\t   -o outputFileName.txt\n";
			
		return 1;
	}	
		opterr=0;
		std::string value;
		int c;
		std::string yOpt;
		std::string fOpt;
		std::string rOpt;
		std::string gOpt;
		std::string eOpt;
		std::string oOpt;
		int yflag=0,fflag=0,rflag=0,gflag=0,eflag=0,oflag=0;
 //////////////////////////
   // get switched args
   //////////////////////////
   while ( (c = getopt(argc, argv, "y:f:r:g:e:o:")) != -1 ) {  // for each option...
        switch ( c ) {
	        case 'y':
                yOpt=optarg;
				yflag++;	
                break;
            case 'f':
                fOpt=optarg;
				fflag++;	
                break;
            case 'r':
                rOpt=optarg;
				rflag++;	
                break;
            case 'g':
                gOpt=optarg;
				gflag++;	
                break;
            case 'e':
                eOpt=optarg;
				eflag++;	
                break;
            case 'o':
                oOpt=optarg;
				oflag++;	
                break;
 
            case '?':  // unknown option...
               if (optopt == 'y'|| optopt == 'f'|| optopt == 'r'|| optopt == 'g'|| optopt == 'e'|| optopt == 'o')
                  fprintf (stderr, "Option -%c requires an argument.\n", optopt);
               else if (isprint (optopt))
                  fprintf (stderr, "Unknown option `-%c'.\n", optopt);
               else
                  fprintf (stderr, "Unknown option character `\\x%x'.\n",optopt);
               return 1;
            default:
               return 1;
        }
    }
		if(yflag+fflag+rflag+gflag+eflag+oflag!=6){
			std::cerr << "All Flags were not entered" << std::endl;
			return 1;
		}
		std::string YFielName =yOpt;
		std::string XFielName =fOpt;
		std::string RFielName =eOpt;
		std::string zFileNameStr =rOpt;
		std::string gFileNameStr =gOpt;
		std::string outputFileName =oOpt;
      
      std::vector<std::string> termsZ=nstring::split(zFileNameStr,':');
      std::vector<std::string> termsG=nstring::split(gFileNameStr,':');

		std::cout << "\nYour inputs:\n";
		std::cout << "Y from the input file " <<  YFielName << std::endl;
		std::cout << "X (fixed design matrix) from the input file " <<  XFielName <<std::endl;
      int q=termsZ.size();
		std::cout << "Z (design matrix) from " << q << " files: " ;
		for(int i = 0; i < q; i++)
		{
			std::cout << termsZ.at(i);
			if(i < q - 1)
			 std::cout << ",";
		}
		std::cout << std::endl;
		std::cout << "G (covariance structure matrix) from " << q << " files: ";
		for(int i = 0; i < q; i++)
		{
			std::cout << termsG.at(i);
			if(i < q - 1)
				std::cout << ",";
		}
      std::cout << std::endl;
		std::cout << "E (residual design matrix) from "  << RFielName <<std::endl;
		std::cout << "Results will be output into "  << outputFileName << std::endl;
		GlmmEngine::doAnalysis(YFielName, XFielName, RFielName, 
				termsZ, termsG, outputFileName);
	
		std::cout << "\nResults were written into the file " << outputFileName << std::endl;
		std::cout << "The GLMM analysis is complete. " << std::endl;
	
      
   return 0;
}
