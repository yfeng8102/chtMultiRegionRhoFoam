#include "createFluidFields.H"  
#include "createSolidFields.H"  

SolverPerformance<vector> solvPerfU;
SolverPerformance<scalar> solvPerfE;
SolverPerformance<scalar> solvPerfh;
//change p_rgh to p
//SolverPerformance<scalar> solvPerfp_rgh;
SolverPerformance<scalar> solvPerfp;
