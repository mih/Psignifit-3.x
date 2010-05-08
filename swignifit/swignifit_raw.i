/*
 *   See COPYING file distributed along with the psignifit package for
 *   the copyright and license terms
 */

/* This is the interface file for the swig wrapper to psignifit, swignifit
 */

%module swignifit_raw

%{
#define SWIG_FILE_WITH_INIT
#include "psipp.h"
%}

// handle all STL exceptions
%include "exception.i"
%exception {
    try {
        $action
    } catch (const std::exception& e) {
        SWIG_exception(SWIG_RuntimeError, e.what());
    }
}

// make the STL vectors available
%include "std_vector.i"
namespace std {
    %template(vector_double) vector<double>;
    %template(vector_int) vector<int>;
};

// This translates BadArgumentError (c++) -> ValueError (python)
// including the error message
%typemap(throws) BadArgumentError %{
      PyErr_SetString(PyExc_ValueError, $1.message);
      SWIG_fail;
%}

%include "std_string.i"

// we need to ignore the second constructor for PsiData since swig can't handle
// this type of overloading TODO write a factory method in python that
// implements this functionality
%ignore PsiData::PsiData (std::vector<double> x,
                          std::vector<int>    N,
                          std::vector<double> p,
                          int nAFC);

// We wrap the following headers
%include "data.h"
%include "psychometric.h"
%include "core.h"
%include "sigmoid.h"
%include "prior.h"
%include "mclist.h"
%include "bootstrap.h"
%include "mcmc.h"
%include "rng.h"
%include "optimizer.h"
