/*
File for the implementation of the class template SnakStruct
this .cpp file is INCLUDED as part of arraystructures.hpp and  cannot be 
compiled on its own.
This file adds the support for a second hashed variable called by KeyParent
*/

#ifndef SNAKSTRUCT_INCL_H_INCLUDED
#define SNAKSTRUCT_INCL_H_INCLUDED

#include "arraystructures.hpp" // never does anything but useful for linter


template<class T> bool SnakStruct<T>::checkready() 
{	
	readyforuse=ArrayStruct<T>::checkready();
	
	readyforuse=(readyforuse && isHashParent);

	return(readyforuse);
	
}

template<class T> void SnakStruct<T>::HashParent()
{
	// Generates a valid unordered_map for the current ArrayStruct
	// Function should not be called repeatedly 
	if(!hashParent.empty()){
		hashParent.clear();
	}
	hashParent.reserve(elems.size());
	for (int i = 0; i < int(elems.size()); ++i)
	{
		hashParent.emplace(elems[i].KeyParent(),i);
	}
	isHashParent=1;
	
	
   //cout << "Array Struct Succesfully Hashed" << endl;
}

template<class T> void SnakStruct<T>::DeHashParent(const int pos){
	
	int key=elems[pos].KeyParent();
	unordered_multimap<int,int>::iterator it = hashParent.find(key);
	while(it->second!=pos && it->first==key){++it;}

	if (it->second==pos && it->first==key){
		hashParent.erase (it);
	} else {

		cerr << "Error: Key value pair not found and could not be removed "<< endl;
		cerr << " key " << key << " pos " << pos << endl;
		cerr << "	in function:" <<  __PRETTY_FUNCTION__ << endl;
		RSVS3D_ERROR_ARGUMENT("Key value pair not found and could not be removed");
	}
}
template<class T> bool SnakStruct<T>::memberIsHashParent(const int pos) const{
	
	int key=elems[pos].KeyParent();
	unordered_multimap<int,int>::const_iterator it = hashParent.find(key);
	if(it != hashParent.end()){
		while(it != hashParent.end() 
			&& it->second!=pos && it->first==key){++it;}

		if (it != hashParent.end() 
			&& it->second==pos && it->first==key){
			return(true);
		} else {

			return(false);
		}
	}else {
		return(false);
	}
}


template<class T> void SnakStruct<T>::PrepareForUse()
{

	ArrayStruct<T>::PrepareForUse();

	if (isHashParent==0){
		this->HashParent();
	}

	readyforuse=true;
}

template<class T> void SnakStruct<T>::Concatenate(const SnakStruct<T> &other)
{
	ArrayStruct<T>::Concatenate(other);
	isHashParent=0;
}

template<class T> void SnakStruct<T>::remove(const vector<int> &sub)
{
	ArrayStruct<T>::remove(sub);
	isHashParent=0;
}

template<class T>
void SnakStruct<T>::ForceArrayReady()
{
	ArrayStruct<T>::ForceArrayReady();
	isHashParent=1;

}
template<class T> inline void SnakStruct <T>::clear()
{
	ArrayStruct<T>::clear();
	hashParent.clear();
	//isHashParent=0; 
}
template<class T> inline void SnakStruct <T>::Init(int n)
{
	this->clear();
	ArrayStruct<T>::Init(n);
 
}
template<class T>
int SnakStruct<T>::findparent(int key) const 
{
	if (isHashParent==0){
		cerr << "Warning: reading from potentially obsolete unordered_map " << endl;
		cerr << "          in snaxarray::findedge(int key)" << endl; 
	}
	auto search=hashParent.find(key);
	

	if (search==hashParent.end()){
		return(-1);
	}
	#ifdef SAFE_ACCESS
	int key2;
	key2=elems[search->second].KeyParent();
	if (key2!=key){
		RSVS3D_ERROR_ARGUMENT("FINDPARENT returned an invalid output ");
	}
	#endif //SAFE_ACCESS
	return(search->second);
}

template<class T>
void SnakStruct<T>::findsiblings(int key, vector<int> &siblings) const 
{
	if (isHashParent==0){
		cerr << "Warning: reading from potentially obsolete unordered_multimap " << endl;
		cerr << "          in " << __PRETTY_FUNCTION__ << endl; 
		cerr << "          To avoid this message perform read operations on ArrayStruct<T> using the () operator" << endl; 
	}
	ReturnDataEqualRange(key, hashParent, siblings);
}
template<class T> inline void SnakStruct <T>::push_back(T& newelem)  
{
	ArrayStruct<T>::push_back(newelem);
	
	hashParent.emplace(newelem.KeyParent(),elems.size()-1);
	//isHashParent=0;
}

// Surfarray methods
template<class T>
void ModiftrackArray<T>::SetNoModif(){
	int ii,n;
	n=ArrayStruct<T>::size();
	for(ii=0;ii<n;++ii){
		elems[ii].isModif=false;
	}
}

template<class T>
void ModiftrackArray<T>::ReturnModifInd(vector<int> &vecind){
	int ii,n;
	n=ArrayStruct<T>::size();
	vecind.clear();
	for(ii=0;ii<n;++ii){
		if(elems[ii].isModif){
			vecind.push_back(elems[ii].index);
		}
	}
}

template<class T>
void ModiftrackArray<T>::ReturnModifLog(vector<bool> &modiflog){
	int ii,n;
	n=ArrayStruct<T>::size();
	modiflog.assign(n,false);
	for(ii=0;ii<n;++ii){

		modiflog[ii]=(elems[ii].index);
	}
}

#endif // ARRAYSTRUCTS_INCL_H_INCLUDED 