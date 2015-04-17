#ifndef NSTRING_H
#define NSTRING_H


#include <string>
#include <iostream>
#include <sstream>
#include <stack>
class nstring{
   private:
   protected:
      int bugFound;

   public:
      virtual ~nstring(){};
     std::string _s;
     nstring(std::string in) {

     _s=in;
     bugFound=setIsErrorCheckFileStack(_s);

     };

   ///////////////////////////////////
   //static
   //
   static inline std::vector<std::string> split(std::string input,const char splitChar){
      std::vector<std::string> returnVector;
		std::stringstream ss(trim(input));
		if(ss){
			std::string token;
            while(std::getline(ss,token,splitChar)){ //columns
				returnVector.push_back(token);
            }

		}
     
      return returnVector;
	}



   static inline std::string replaceWordCopy(const std::string inW,const std::string searchValue,std::string replace){
        std::string returnString=inW;
       while(returnString.find(searchValue)!=std::string::npos)
            returnString.replace(returnString.find(searchValue),searchValue.length(),replace);

        return returnString;

   }
   static inline bool setIsErrorCheckFileStack(std::string in){
      std::string inFile=removeAllIn(in);
      std::stack<char> stackOfChars;

      char checkString;
      int bugFound=false;
      unsigned int i;
      for(i=0;i<inFile.length() && bugFound !=true;i++){

         switch (inFile[i])
         {

            case '(':
               if(i>0 && inFile[i-1]!='\\'){
                     stackOfChars.push('(');
               }
					else
                     stackOfChars.push('(');
               break;
            case '{':
               if(i>0 && inFile[i-1]!='\\'){
                  stackOfChars.push('{');
               }
					else
                  stackOfChars.push('{');
               break;
            case ')':
               if(i>0 && inFile[i-1]!='\\'&& stackOfChars.size()!=0){
                  checkString=stackOfChars.top();
               }
               if(checkString!='('||stackOfChars.size()==0){
                    std::cout << "Found Errors in matching parenthesis" << std::endl;
                  bugFound=true;
               }
               if(stackOfChars.size()!=0)
                stackOfChars.pop();
               break;
            case '}':
               if(i>0 && inFile[i-1]!='\\' && stackOfChars.size()!=0){
                    
                  checkString=stackOfChars.top();
               }
               if(checkString!='{' || stackOfChars.size()==0){
                    std::cout << "Found Errors in matching curly braces" << std::endl;
                  bugFound=true;
               }
               if(stackOfChars.size()!=0)
                stackOfChars.pop();
               break;
        }

      }
      if(!stackOfChars.empty()){
         bugFound=true;
        std::cout << "Error in matching curly braces or parenthesis" << std::endl;
//     		throw(std::string( "Error in script\n "));
      }
      return bugFound;

   }
	static inline std::string removeAllIn(std::string in){

	   in=nstring::replaceIn(in,'\'');
	   in=replaceIn(in,'`');
	   in=replaceIn(in,'\"');
		return in;

	}
   static inline std::string removeComments(std::string &in){
      unsigned int intIteratorIt=0;
      std::string comments;

      int count =0;
      for(;intIteratorIt<in.length();intIteratorIt++){
       if(in[intIteratorIt]=='/'){
            if(in[intIteratorIt+1]=='/'){
               while (in[intIteratorIt]!='\n' && intIteratorIt<in.length()){
       //             std::cout << "iterator:" << intIteratorIt << std::endl;
        //            std::cout << "in.length(): " << in.length() << std::endl;
         //           std::cout << "atposition: " <<in[intIteratorIt] << std::endl;
                  comments=comments+in[intIteratorIt];
                  in[intIteratorIt]=' ';
                  intIteratorIt++;
               }
               comments=comments+'\n';
            }
         }

          count ++;
      }
      return comments;
   }

   static inline std::string ltrim(std::string ins,const char* t = " \t\n\r\f\v"){
       return(ins.erase(0, ins.find_first_not_of(t)));
   }

   static inline std::string rtrim(std::string ins,const char* t=" \t\n\r\f\v"){
       return(ins.erase(ins.find_last_not_of(t) + 1));
   }
   static inline std::string trim(std::string ins,const char* t = " \t\n\r\f\v")
   {
         return rtrim(ltrim(ins,t),t);

   }


   static inline std::string replaceIn(std::string ins,const char find ='\''){
      std::string::iterator it=ins.begin();
      for(;it!=ins.end();it++){
         char value=*it;
         if(value==find ){
            it++;
            char value2=*it;
			//	std::cout << "value2=" << value2 << std::endl;
            char value2sub1=*(it-1);
			//	std::cout << "value2sub1=" << value2sub1 << std::endl;
            while(value2!=find && value2sub1!='\\'){
   		//		std::cout << "changing=" << value2 << std::endl;
               *it=' ';
               it++;
               value2=*it;
             }
         }
      }
      return ins;
   }

   ///////////////////////////////////
   // trim from left
   inline void ltrim(const char* t = " \t\n\r\f\v")
   {
         _s=ltrim(_s,t);
//       _s.erase(0, _s.find_first_not_of(t));
   }

   // trim from right
   inline void rtrim( const char* t = " \t\n\r\f\v")
   {
         _s=rtrim(_s,t);
   }

   // trim from left & right
   inline void trim(const char* t = " \t\n\r\f\v")
   {
         _s=trim(_s,t);
   }

   inline void replaceIn(const char find ='\''){
      _s=replaceIn(_s,find);
   }

   inline std::string removeComments(){
     std::string comments= removeComments(_s);
      return comments;
  }

   inline bool setIsErrorCheckFileStack(){
            bugFound=setIsErrorCheckFileStack(_s);
           return (bugFound);
   }

};
#endif
