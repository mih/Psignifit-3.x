#ifndef GETSTART_H
#define GETSTART_H

#include "psychometric.h"
#include "data.h"

#include <list>
#include <vector>

/** \brief scalable and shiftable grid object used for the gridsearch */
class PsiGrid {
	private:
		std::list< std::vector<double> > grid1d;
		std::vector<double> lower_bounds;
		std::vector<double> upper_bounds;
	public:
		PsiGrid ( void ) {} ///< an empty grid
		PsiGrid (
				std::vector<double> xmin,  ///< lowest values for the parameters --- this should have the same length as the parameters we want to feed into the model
				std::vector<double> xmax,  ///< highest values for the parameres --- this should have the same length as the parameters we want to feed into the model
				unsigned int gridsize      ///< size of the grid in every dimension
				); ///< Generate a proper grid
		PsiGrid& shift ( std::vector<double> newposition ) const;  ///< shift the grid to be centered on newposition
		PsiGrid& shrink ( std::vector<double> newposition ) const; ///< shrink the grid around newposition
		PsiGrid& subgrid ( void ) const;                           ///< return a subgrid with the first dimension eliminated
		unsigned int get_gridsize ( void ) const { return grid1d.front().size(); } ///< return the size of the grid in every dimension
		bool empty ( void ) const { return grid1d.empty(); }   ///< check whether the grid is empty i.e. does not have any grid points
		const std::vector<double>& front ( void ) const { return grid1d.front(); }   ///< return the values along the first dimension of the grid
		unsigned int dimension ( void ) const { return grid1d.size(); }              ///< get dimension of the grid i.e. how many parameters are varied
		double get_lower ( unsigned int i ) const { return lower_bounds[i]; }        ///< get lower bound for parameter i
		double get_upper ( unsigned int i ) const { return upper_bounds[i]; }        ///< get upper bound for parameter i
};


std::vector<double> linspace (
		double xmin,
		double xmax,
		unsigned int n
		);   ///< create n linearly spaced values between xmin and xmax

void makegridpoints (
		const PsiGrid& grid,                           ///< PsiGrid object on from which the grid points should be generated
		std::vector<double> prm,                       ///< current setting of the parameter vector (irrelvant for initial call)
		unsigned int pos,                              ///< position of the currently modified parameter (should be 0 for initial call)
		std::list< std::vector<double> > *gridpoints   ///< gridpoints list to which new gridpoints are added
		);    ///< generate new grid points based on a PsiGrid object

void evalgridpoints (
		std::list< std::vector<double> > *grid,       ///< gridpoints on which the negative log posterior should be evaluated
		std::list< std::vector<double> > *bestprm,    ///< list to which the best parameters should be appended
		std::list< double > *L,                       ///< list to which the best negative log posteriors should be appended
		const PsiData* data,                          ///< data set on which the gridsearch should be performed
		const PsiPsychometric* pmf,                   ///< psychometric function for which the gridsearch should be performed
		unsigned int nbest,                           ///< how many "best" parameter constellations should be determined?
		);              ///< evaluate negative log posterior on a gridpoints and update bestprm and L to refer to the nbest parameter settings encountered so far

void updategridpoints (
		const PsiGrid& grid,                              ///< PsiGrid obejct from which the new grids and gridpoints should be generated
		const std::list< std::vector<double> >& bestprm,  ///< so far best parameter settings
		std::list< std::vector<double> > *newgridpoints,  ///< new gridpoints to be evaluated
		std::list< PsiGrid > *newgrids                    ///< grids corresponding to the new gridpoints
		);             ///< generate new gridpoints by shifting and shrinking the previous grid and store the newly generated PsiGrid objects

void a_range ( const PsiData* data, double *xmin, double *xmax );     ///< Heuristic to generate a first idea of the threshold parameter a
void b_range ( const PsiData* data, double *xmin, double *xmax );     ///< Heuristic to generate a first idea of the width parameter b
void lm_range ( const PsiData* data, double *xmin, double *xmax );    ///< Heuristic to generate a first idea of the lapse rate parameter lm
void gm_range ( const PsiData* data, double *xmin, double *xmax );    ///< Heuristic to generate a first idea of the guessing rate parameter gm

void parameter_range ( const PsiData* data, unsigned int prmindex, double *xmin, double *xmax );   ///< call the proper heuristic associated with the parameter at prmindex to generate a range for the respective parameter

std::vector<double> getstart (
		const PsiPsychometric* pmf,    ///< psychometric function model for which a starting value is desired
		const PsiData* data,           ///< data for which a starting value is desired
		unsigned int gridsize,         ///< number of grid points to be used
		unsigned int nneighborhoods,   ///< number of neighborhoods to be studied
		unsigned int niterations       ///< number of iterated neighborhood searched to be performed
		);    ///< Determine a good starting value using nested grid search

#endif
