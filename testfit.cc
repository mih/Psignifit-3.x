#include "psychometric.h"
#include <iostream>
#include <vector>

using namespace std;

int main ( int argc, char ** argv ) {
	vector<double> x (6);
	vector<int>    n (6,50);
	vector<int>    k (6);

	x[0] = 0;  x[1] = 2;  x[2] = 4;  x[3] = 6;  x[4] = 8;  x[5] = 10;
	k[0] = 33; k[1] = 39; k[2] = 39; k[3] = 44; k[4] = 48; k[5] = 48;

	PsiData *data = new PsiData (x,n,k,2);

	PsiPsychometric *pmf = new PsiPsychometric (2, new abCore, new PsiLogistic);
	pmf->setPrior( 2, new UniformPrior(0.,0.1));

	PsiOptimizer *opt = new PsiOptimizer (pmf, data);

	vector<double> solution (4);

	solution = opt->optimize(pmf,data);

	cerr << solution[0] << " " << solution[1] << " " << solution[2] << " " << solution[3] << "\n";
}
