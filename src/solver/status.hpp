#ifndef STATUS
#define STATUS

#include "gurobi_c++.h"

enum Status{
    Not_Solved = 0,
    Unknown = 1,
    Optimal = 2,
    Time_Limit = 3,
    Unbounded = 4,
    Infeasible = 5,
    Infeasible_or_Unbounded = 6
};

Status statusFromGurobi(int gurobi_status);
std::ostream& operator<<(std::ostream& lhs, const Status & status);

#endif