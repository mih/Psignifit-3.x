#ifndef PYTOOLS_H
#define PYTOOLS_H

#include <string>

PsiData * create_dataset ( PyObject * pydata, int Nafc, int *nblocks ) {
	if ( !PySequence_Check ( pydata ) )
		throw std::string ( "data should be a sequence" );

	PyObject * pyblock;

	int Nblocks ( PySequence_Size ( pydata ) );
	*nblocks = Nblocks;
	std::vector<double> x ( Nblocks );
	std::vector<int> k ( Nblocks );
	std::vector<int> n ( Nblocks );

	for ( int i=0; i<Nblocks; i++ ) {
		pyblock = PySequence_GetItem ( pydata, i );
		if ( PySequence_Size ( pyblock ) != 3 ) {
			char msg[50];
			sprintf ( msg,"data in block %d do not have 3 entries", i );
			throw std::string ( msg );
		}
		x[i] = PyFloat_AsDouble ( PySequence_GetItem ( pyblock, 0 ) );
		k[i] = PyInt_AsLong ( PySequence_GetItem ( pyblock, 1 ) );
		n[i] = PyInt_AsLong ( PySequence_GetItem ( pyblock, 2 ) );
		std::cerr << i << " " << x[i] << " " << k[i] << " " << n[i] << "\n";
	}

	return new PsiData ( x, n, k, Nafc );
}

#endif
