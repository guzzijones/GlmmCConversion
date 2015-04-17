#ifndef JARRAY_H
#define JARRAY_H
#include <vector>
#include <string>

class JArray{
	private:

	protected:

	public:

	static double mean(const std::vector<std::vector<double>>& yM, int columnIndex);
	static double mean(const std::vector<double> & yM);

	static double var(const std::vector<std::vector<double>> & yM,  int columnIndex);
	static double var(const std::vector<double> & yM);

	static double std(std::vector<std::vector<double>> yM, int columnIndex);
	static double std(std::vector<double> yM);

	static double min(std::vector<double> yM);
	static int max(std::vector<int> yM);
//	static double std::vector<std::vector<double>>  convert(ArrayList<Double> y)

	static std::vector<std::string> getUniq(std::vector<std::string> nameList);
	static void demo(std::vector<double> yM);
	static void demo(std::vector<std::string> aList);


};


#endif
