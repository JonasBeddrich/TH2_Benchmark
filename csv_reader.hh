#include <string>
#include <fstream>
#include <vector>
#include <utility>
#include <stdexcept>
#include <sstream>
#include <iostream>

class linear_interpolator {
private:
	std::vector<double> xData;
	std::vector<double> yData;
	int size;
	
	
public:
	linear_interpolator(){}
	
	linear_interpolator(std::vector<double> x, std::vector<double> y) {
		set_x_and_y(x,y); 
	}
	
	linear_interpolator(double point_distance, std::vector<double> y){
		set_x_and_y(point_distance, y); 
	}
	
	
	void set_x_and_y(std::vector<double> x, std::vector<double> y){
		xData = x;
		yData = y;
		size = xData.size();
	}

	void set_x_and_y(double point_distance, std::vector<double> y){
		yData = y; 
		size = yData.size(); 
		xData = std::vector<double>(size); 
		for (int i = 0; i < xData.size(); i++) {
			xData[i] = i * point_distance;
		}
	}

	double interpolate(double x) {
		int i = 0;
		if (x >= xData[size - 2]) {
			i = size - 2;
		} else {
			while (x > xData[i + 1])
				i++;
		}
		double xL = xData[i], yL = yData[i];
		double xR = xData[i + 1], yR = yData[i + 1];

		
		double dydx = (yR - yL) / (xR - xL);
		return yL + dydx * (x - xL);
	}
};


class csv_reader {
private:

public:
	csv_reader() {
	}
	// This one excepts a non double first line 
	std::vector<double> read_1D_axis(std::string filename) {
		std::vector<double> result;
		std::ifstream myFile(filename);
		std::string line;

		// just checkin 
		if (!myFile.is_open())
			throw std::runtime_error("Could not open file: " + filename);
		// okay everything is fine 
		if (myFile.good()) {
			int counter = -1;
			while (std::getline(myFile, line)) {
				if (counter == -1) {
					counter += 1;
				} else {
					double temp = ::atof(line.c_str());
					result.push_back(temp);
				}
			}
		}
		myFile.close();
		return result;
	}

	std::vector<std::vector<double>> read_2D_matrix(std::string filename) {
		std::vector<std::vector<double>> result;
		std::ifstream myFile(filename);
		std::string line; 

		if (!myFile.is_open())
			throw std::runtime_error("Could not open file: " + filename);
		if (!myFile.good())
			throw std::runtime_error("File could not be read: " + filename);
		
		while (std::getline(myFile, line)) {
			std::stringstream lineStream(line);
			std::string cell;
			std::vector<double> row;
			
			while (std::getline(lineStream, cell, ',')) {
				double cell_as_number = std::stod(cell); /*sTRING to dOUBLE*/ 
				row.push_back(cell_as_number);
			}
			result.push_back(row);
		}	
		myFile.close();
		return result;
	}
};