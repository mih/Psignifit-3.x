%module swignifit

%{
#define SWIG_FILE_WITH_INIT
#include "psipp.h"

%}

%include "std_vector.i"

namespace std {
    %template(vectord) vector<double>;
};

%include "psychometric.h"
%include "core.h"
%include "sigmoid.h"
    
