#ifndef BOOTSTRAP_H
#define BOOTSTRAP_H

#include <vector>
#include <cmath>
#include "psychometric.h"
#include "mclist.h"

BootstrapList parametricbootstrap ( int B, const PsiData * data, const PsiPsychometric* model, std::vector<double> cuts, bool BCa=true );
JackKnifeList jackknifedata ( const PsiData * data, const PsiPsychometric* model );

#endif
