#include <string>
#include <iostream>
#include <fstream>
#include <exception>

#include "input/input.hpp"
#include "solver/solverfactory.hpp"
#include "solver/abstractsolver.hpp"

int main(int argc, char *argv[]) {
	/* Data information. */
	try{
		Input input(argc,argv);
		input.writeHead();
		input.display();
		SolverFactory factory;
		AbstractSolver *solver = factory.createSolver(input);
		solver->solve();
		input.write(solver->write());
		delete solver;
	}catch (GRBException e){
    	std::cout << "Error code = " << e.getErrorCode() << std::endl;
    	std::cout << e.getMessage() << std::endl;
  	}
	catch (const std::string & e) { std::cout << "EXCEPTION | " << e << std::endl; }
	catch (const std::exception & e) { std::cout << "EXCEPTION | " << e.what() << std::endl; }

    return 0; 
}