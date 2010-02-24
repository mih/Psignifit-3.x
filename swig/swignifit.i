%module swignifit

%{
#define SWIG_FILE_WITH_INIT
#include "psipp.h"

%}

%include "std_vector.i"

namespace std {
    %template(vector_double) vector<double>;
};

%typemap(throws) BadArgumentError %{
      PyErr_SetString(PyExc_ValueError, "Contents");
      SWIG_fail;
%}

%include "psychometric.h"
%include "core.h"
%include "sigmoid.h"
    
