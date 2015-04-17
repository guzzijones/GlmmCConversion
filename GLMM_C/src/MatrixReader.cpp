#include "MatrixReader.h"

using namespace std;
MatrixReader::MatrixReader(std::string in){
	_fileName=in;
}

const std::vector < std::vector<double> > & MatrixReader::getData()const {
	return _X;
}
void MatrixReader::read(){
	std::string strSplitter=",";
	std::unique_ptr<std::ifstream> fin=FileIOHelper::openFileForRead(_fileName);
	std::string sentence;
	std::vector< std::vector<double> > dataM;
	
	//initial row vector of doubles
	//preread
	sentence=FileIOHelper::readLine(fin);
	sentence=nstring::trim(sentence,",");
	sentence=nstring::trim(sentence);
	while(sentence.length()!=0){ //rows
		std::vector<double> tokensDoubles;
		std::stringstream ss(nstring::trim(sentence));
		if(ss){
			std::string token;
            while(std::getline(ss,token,',')){ //columns
				std::istringstream attrSS(token);
				double attrD;
				attrSS >> attrD;
				tokensDoubles.push_back(attrD);
            }
		}
		if(tokensDoubles.size()!=0)
			_X.push_back(tokensDoubles);//push row
		sentence=nstring::trim(FileIOHelper::readLine(fin),",");
		sentence=nstring::trim(sentence);
	}

	FileIOHelper::closeFile(fin);
}


