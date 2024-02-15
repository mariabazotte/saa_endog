#include "status.hpp"

Status statusFromGurobi(int gurobi_status) {
    if(gurobi_status == GRB_OPTIMAL){
        return Status::Optimal;
    }else if(gurobi_status == GRB_TIME_LIMIT){
        return Status::Time_Limit;
    }else if(gurobi_status == GRB_UNBOUNDED){
        return Status::Unbounded;
    }else if(gurobi_status == GRB_INF_OR_UNBD){
        return Status::Infeasible_or_Unbounded;
    }else if(gurobi_status == GRB_INFEASIBLE){
        return Status::Infeasible;
    }
    return Status::Unknown;
} 


std::ostream& operator<<(std::ostream& lhs, const Status & status) {
    switch(status) {
        case Status::Not_Solved: {
            lhs << "NotSolved";
            break;
        }
        case Status::Unknown: {
            lhs << "Unknown";
            break;
        }
        case Status::Optimal: {
            lhs << "Optimal";
            break;
        }
        case Status::Time_Limit: {
            lhs << "TimeLimit";
            break;
        }
        case Status::Unbounded: {
            lhs << "Unbounded";
            break;
        }
        case Status::Infeasible: {
            lhs << "Infeasible";
            break;
        }
        case Status::Infeasible_or_Unbounded: {
            lhs << "InfeasibleOrUnbounded";
            break;
        }
        default :{
            lhs << "";
            break;
        }
    }
    return lhs;
}