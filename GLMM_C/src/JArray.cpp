#include "JArray.h"
double JArray::mean(const std::vector<std::vector<double>> & yM, int columnIndex){
    double   	mu = 0;
	    int
	    	n = yM.size();
		for (int i = 0; i < n; i++)
			mu += yM[i][columnIndex];
		mu /= (double)n;
		return mu;

}

double JArray::mean(const std::vector<double> & yM){
	{
	    double
	    	mu = 0;
	    int
	    	n = yM.size();
		for (int i = 0; i < n; i++)
			mu += yM[i];
		mu /= (double)n;
		return mu;
	}


}

double JArray::var(const std::vector<std::vector<double>> & yM,  int columnIndex){
	double
			mu = mean(yM, columnIndex),
			var = 0,
			temp = 0;
		int
			n = yM.size(),
			df = n;
		for (int i = 0; i < n; i++)
		{
			temp = yM[i][columnIndex] - mu;
			var += (temp * temp);
		}//End Of All the Individuals.
		var /= (double)df;
		return var;

}
double JArray::var(const std::vector<double> & yM){
	double
			mu = mean(yM),
			var = 0,
			temp = 0;
		int
			n = yM.size(),
			df = n - 1;
		for (int i = 0; i < n; i++)
		{
			temp = yM[i]- mu;
			var += (temp * temp);
		}//End Of All the Individuals.
		var /= (double)df;
		return var;

}
double JArray::std(std::vector<std::vector<double>> yM, int columnIndex){

}
double JArray::std(std::vector<double> yM){

}

double JArray::min(std::vector<double> yM){

}
int JArray::max(std::vector<int> yM){

}
//	static double std::vector<std::vector<double>>  convert(ArrayList<Double> y)

std::vector<std::string> JArray::getUniq(std::vector<std::string> nameList){

}
void JArray::demo(std::vector<double> yM){

}
void JArray::demo(std::vector<std::string> aList){

}


