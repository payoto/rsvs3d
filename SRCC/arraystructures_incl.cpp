/*
File for the implementation of the class template ArrayStruct
this .cpp file is INCLUDED as part of arraystructures.hpp and  cannot be 
compiled on its own.

*/

#ifndef ARRAYSTRUCTS_INCL_H_INCLUDED
#define ARRAYSTRUCTS_INCL_H_INCLUDED

template <class T> bool CompareDisp(T *mesh1,T *mesh2){
   bool compFlag;
   stringstream ss1,ss2;
   auto old_buf = cout.rdbuf(ss1.rdbuf()); 

   mesh1->disp();
   cout.rdbuf(ss2.rdbuf()); 
   mesh2->disp();
   std::cout.rdbuf(old_buf);

   compFlag=ss1.str().compare(ss2.str())==0;
   return(compFlag);
}

// member function definition template <class T> : "ArrayStruct"

template<class T> inline int ArrayStruct <T>::GetMaxIndex() const 
{
	return(maxIndex);
}

template<class T> inline int ArrayStruct <T>::find(int key) const 
{
	auto search=hashTable.find(key);
	if (search==hashTable.end()){
		return(-1);
	}
	return(search->second);
}

template<class T> inline int ArrayStruct <T>::size() const 
{
	return(int(elems.size()));
}



template<class T> void ArrayStruct <T>::disp() const 
{
	for (int ii=0 ; unsigned_int(ii)<elems.size();ii++){
		cout << "Array " << ii << " " ;
		elems[ii].disp();
	}
}


template<class T> inline void ArrayStruct <T>::Init(int n)
{
	T sT;
	elems.reserve(n);
	elems.assign(n,sT);
}


template<class T> void ArrayStruct <T>::HashArray()
{
	hashTable.reserve(elems.size());
	for (int i = 0; i < int(elems.size()); ++i)
	{
		hashTable.emplace(elems[i].Key(),i);
	}
   //cout << "Array Struct Succesfully Hashed" << endl;
}

template<class T> void ArrayStruct <T>::Concatenate(const ArrayStruct <T> &other)
{
	int nCurr,nNew,nTot,ii;
	nCurr=this->size();
	nNew=other.size();
	nTot=nCurr+nNew;
	elems.reserve(nTot);

	for (ii=0; ii<nNew;ii++){
		elems.push_back(other.elems[ii]);
	}
	this->HashArray();
}

template<class T> void ArrayStruct <T>::SetMaxIndex()
{
	int n=elems.size();
	maxIndex=-1;
	for(int ii=0;ii<n;++ii){
		maxIndex= (maxIndex<this->elems[ii].index) ? elems[ii].index : maxIndex;
	}
}

template<class T> void ArrayStruct <T>::PopulateIndices()
{
	int n=elems.size();
	maxIndex=-1;
	for(int ii=0;ii<n;++ii){
		elems[ii].index=ii+1;
	}
	maxIndex=n;
}

#endif // ARRAYSTRUCTS_INCL_H_INCLUDED