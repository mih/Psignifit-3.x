%module swignifit

%{
#define SWIG_FILE_WITH_INIT
#include "psipp.h"

%}

%include "std_vector.i"

namespace std {
    %template(vector_double) vector<double>;
    %template(vector_int) vector<int>;
};

%typemap(throws) BadArgumentError %{
      PyErr_SetString(PyExc_ValueError, "Contents");
      SWIG_fail;
%}

// we need to ignore the second constructor for PsiData since swig can't handle
// this type of overloading TODO write a factory method in python that
// implements this functionality
%ignore PsiData::PsiData (std::vector<double> x,
                          std::vector<int>    N,
                          std::vector<double> p,
                          int nAFC);

%include "data.h"
%include "psychometric.h"
%include "core.h"
%include "sigmoid.h"
    
