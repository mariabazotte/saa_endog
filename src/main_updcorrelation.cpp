#include <string>
#include <iostream>
#include <fstream>
#include <exception>

#include "input/input.hpp"
#include "tests/updcorrelation.hpp"

int main(int argc, char *argv[]) {
	/* Data information. */
	try{
		int mandatory = 0;
		std::string parameters_file;
		std::string network_file;
		int nb_samples = 0;
		int size_samples = 0;
		int nb_solutions = 0;
		int distribution = 0;
		for (int i = 1; i < argc; i += 2){
			if(std::string(argv[i]) == "-paramfile"){ // Mandatory parameters
				parameters_file = argv[i+1];
				mandatory += 1;
			}else if(std::string(argv[i]) == "-networkfile"){
				network_file = argv[i+1];
				mandatory += 1;
			}else if(std::string(argv[i]) == "-nbsamples"){
				nb_samples = std::stoi(argv[i+1]);
				mandatory += 1;
			}else if(std::string(argv[i]) == "-sizesamples"){
				size_samples = std::stoi(argv[i+1]);
				mandatory += 1;
			}else if(std::string(argv[i]) == "-d"){
				distribution = std::stoi(argv[i+1]);
				mandatory += 1;
			}else if(std::string(argv[i]) == "-nbsolutions"){
				nb_solutions = std::stoi(argv[i+1]);
				mandatory += 1;
			}
		}
		if(mandatory != 6){
			std::cerr << "ERROR: Not all mandatory arguments were defined." << std::endl;
			std::cerr << "You need to define: -paramfile -networkfile -nbsamples -sizesamples -d -nbsolutions -randomlygenerate";
			throw std::runtime_error(std::string("Incorrect line of command"));
		}
		UpdCorrelation correlation(parameters_file,network_file,distribution,nb_samples,size_samples,nb_solutions);
		correlation.solve();

	}catch (GRBException e){
    	std::cout << "Error code = " << e.getErrorCode() << std::endl;
    	std::cout << e.getMessage() << std::endl;
  	}
	catch (const std::string & e) { std::cout << "EXCEPTION | " << e << std::endl; }
	catch (const std::exception & e) { std::cout << "EXCEPTION | " << e.what() << std::endl; }

    return 0; 
}