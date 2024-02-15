#include "lshaped.hpp"

void LShaped::solve(AbstractMainProblem* main){
    start_time = time(NULL); 
    main->setStartTime(start_time); 

    /** Loop for Hot Start. */
    double LLB = -GRB_INFINITY;
    main->setLinearRelaxation();

    while(main->getNbHotStartIt() < input.getHotStartMaxIt()){
        main->increaseNbHotStartIt();
        main->solve();
        if(main->getStatus() == Status::Infeasible){ 
            std::cerr << "The relaxed problem is infeasible." << std::endl; 
            return; 
        }
        main->solveSubProblems(false); 
        if(input.getPapadakos() == true) { 
            if(main->getNbHotStartIt() == 1){
                main->initializeCorePoint();
            }else{ 
                main->updateCorePoint();
                main->solveSubProblems(true);
            }
        }
        if(input.getVerbose()>=1){
            printf ("i %3d | %12.4f | %3d | %3d | %12.4f \n", main->getNbHotStartIt(), main->getLB(), main->getNbViaCuts(), main->getNbOptCuts(), time(NULL)-start_time);
        }if (main->getNbHotStartIt() >= 60 && (main->getLB() - LLB) < input.getPrecision()){
            break;
        }else{
            LLB = main->getLB();
        }
    }
    if(main->getStatus() == Status::Infeasible){ return; }

    /** Integer Loop. */ 
    main->setInteger();
    if(input.getUseCallback() == true){
       main->solveCallback();
    }else{
        main->setSup(0.0);
        bool stop = false;
        
        main->increaseNbIt();
        if(input.getPapadakos() == true){
            main->solve(); 
            if(main->getStatus() == Status::Infeasible){ return; }
            main->initializeCorePoint(); 
        }
        if(main->getTimeLimit() <= 0.0){
            main->setStatus(Status::Time_Limit);
            stop = true;
        }
        while(stop == false){
            if(input.getPapadakos() == true) main->solveSubProblems(true);  
            main->solve(); 
            if(main->getStatus() == Status::Infeasible){ return; } 
            main->solveSubProblems(false);  
            if(input.getPapadakos() == true) main->updateCorePoint();
            
            if(input.getTrustRegion() == true && main->getNbIt() == 1){ 
                main->createTrustRegion();
                main->updateTrustRegion();
            }else if(input.getTrustRegion() == true && main->getNbIt() <= 10){
                main->updateTrustRegion();
            }else if(input.getTrustRegion() == true && main->getNbIt() == 11){
                main->removeTrustRegion();
            }

            if(main->getTimeLimit() <= 0.0){
                main->setStatus(Status::Time_Limit);
                stop = true;
            }else{
                main->updateUB();
                main->updateGAP();
            }
            if((main->getGap() < input.getPrecision()) && (main->getSup() < GRB_INFINITY)){
                main->setStatus(Status::Optimal);
                stop = true;
            }
            if(input.getVerbose() >= 1){
                if(main->getUB() >= GRB_INFINITY && main->getLB() <= -GRB_INFINITY){
                    printf ("i %3d | -infinity | infinity | %12.4f | %3d | %3d | %3d | %3d | %12.4f \n", main->getNbIt(), main->getGap(), main->getNbViaCuts(), main->getNbOptCuts(), main->getNbPpdViaCuts(), main->getNbPpdOptCuts(), time(NULL)-start_time);
                }else if(main->getUB() >= GRB_INFINITY ){
                    printf ("i %3d | %12.4f | infinity | %12.4f | %3d | %3d | %3d | %3d | %12.4f \n", main->getNbIt(), main->getLB(), main->getGap(), main->getNbViaCuts(), main->getNbOptCuts(), main->getNbPpdViaCuts(), main->getNbPpdOptCuts(), time(NULL)-start_time);
                }else if(main->getLB()<= -GRB_INFINITY){
                    printf ("i %3d | -infinity | %12.4f | %12.4f | %3d | %3d | %3d | %3d | %12.4f \n", main->getNbIt(), main->getUB(), main->getGap(), main->getNbViaCuts(), main->getNbOptCuts(), main->getNbPpdViaCuts(), main->getNbPpdOptCuts(), time(NULL)-start_time);
                } else{
                    printf ("i %3d | %12.6f | %12.6f | %12.4f | %3d | %3d | %3d | %3d | %12.4f \n", main->getNbIt(), main->getLB(), main->getUB(), main->getGap(), main->getNbViaCuts(), main->getNbOptCuts(), main->getNbPpdViaCuts(), main->getNbPpdOptCuts(), time(NULL)-start_time);
                }
            }
            main->increaseNbIt();
        }
    }
}

void LinearizationLShaped::create(){
    main = new LinearizationMainProblem(instance,input);
}

void LinearizationLShaped::solve(){
    double solve_start_time = time(NULL);
    lshaped->solve(main);
    lb = main->getLB();
    ub = main->getUB();
    gap = main->getGap();
    status = main->getStatus();
    setTime(time(NULL) - solve_start_time);
}