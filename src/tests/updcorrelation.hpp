#ifndef UpdCorrelation_HPP
#define UpdCorrelation_HPP

#include "../input/input.hpp"
#include "../solver/abstractsolver.hpp"
#include "../solver/saa/saatestproblem.hpp"
#include "../data/instance.hpp"
#include "../data/saainstance.hpp"
#include "../data/allsceinstance.hpp"

class UpdCorrelation {
    protected:
        std::string params_file;           /* Parameters file. */
        std::string network_file;          /* Network file. */
        Input input;                       /* Input. */

        Instance * instance;
        SAAInstance *saainstance;                 /* SAA instance with sampled scenarios. */
        SAATest *saatest;                         /* SAA test problem to obtain the sample estimators value.*/

        int distribution;
        int nb_samples;
        int size_samples;
        int nb_solutions;

        std::vector<double*> x_;
        std::vector<double*> y_;
        std::vector<double*> mean_;
        std::vector<double*> stddev_;
        std::vector<double*> prob_;
        std::vector<double**> utility_;
        std::vector<std::vector<double>> objective;

        std::string solution_file;


        void binomialDistribution();
        void normalDistribution();
        void stdNormalization();
        void defineSolutionFile();
        void solveSAATest();
        void computeAndPrintEstimators();

        void initiateInput();
        void randomlyGenerateSolutions();
        void generateSolutions();
        void create();
        int nb_facilities;
        int nb_protec;
    
    public:
        UpdCorrelation(std::string, std::string, int, int, int, int);
        ~UpdCorrelation();

        void solve();
        std::vector<std::string> write() const { std::vector<std::string> out; return out; }
};


#endif
