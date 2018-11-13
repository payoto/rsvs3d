
/*
File for the implementation of the class template vectorarray
this .cpp file is INCLUDED as part of vectorarray.hpp and  cannot be 
compiled on its own.
*/


#ifndef VECTORARRAY_INCL_H_INCLUDED
#define VECTORARRAY_INCL_H_INCLUDED
#include "vectorarray.hpp" // include guarded does nothing (needed for the linter)

template<class T> void ArrayVec<T>::assign(int nR,int nC, T newelem)  
{
	vector<T> tempElems;
	tempElems.assign(nC,newelem);
	elems.assign(nR,tempElems);

	dim.clear(); 
	dim.push_back(nR); 
	dim.push_back(nC);
}

#endif // VECTORARRAY_INCL_H_INCLUDED