/*
File for the implementation of the class template ArrayStruct
this .cpp file is INCLUDED as part of arraystructures.hpp and  cannot be 
compiled on its own.

*/
#ifndef ARRAYSTRUCTS_INCL_H_INCLUDED
#define ARRAYSTRUCTS_INCL_H_INCLUDED

#include "arraystructures.hpp"
template<class T> void ConcatenateVector(vector<T> &vecRoot, const vector<T> &vecConcat)
{
	vecRoot.insert(vecRoot.end(),vecConcat.begin(),vecConcat.end());
}

template <class T> bool CompareDisp(T &mesh1,T &mesh2)
{
	bool compFlag;
	stringstream ss1,ss2;
	auto old_buf = cout.rdbuf(ss1.rdbuf()); 

	mesh1.disp();
	cout.rdbuf(ss2.rdbuf()); 
	mesh2.disp();
	std::cout.rdbuf(old_buf);

	compFlag=ss1.str().compare(ss2.str())==0;
	return(compFlag);
}

template<class T> void DisplayVector(vector<T> vec)
{
	cout << int(vec.size()) << " - "; 
	for (int i = 0; i < int(vec.size()); ++i)
	{
		cout << vec[i] << " ";
	}
	cout << " | " ;
}
template<class T, class R> 
R ConcatenateVectorField(const ArrayStruct<T> &arrayIn, R T::*mp, const vector<int> &subList)
{
	R surfInds;
	int ii;
	auto itVecInt=surfInds.begin();

	for(ii=0; ii<int(subList.size());++ii){
		surfInds.insert(itVecInt, (arrayIn(subList[ii])->*mp).begin(),
			(arrayIn(subList[ii])->*mp).end());
		itVecInt=surfInds.end();
		//vDisplayVector(arrayIn(subList[ii])->*mp);cout << endl;
	}
	return(surfInds);
}

template<class T, class R> 
vector<R> ConcatenateScalarField(const ArrayStruct<T> &arrayIn, R T::*mp, const vector<int> &subList)
{
	vector<R> surfInds;
	int ii;
	surfInds.reserve(subList.size());
	for(ii=0; ii<int(subList.size());++ii){
		surfInds.push_back((arrayIn(subList[ii])->*mp));
		//vDisplayVector(arrayIn(subList[ii])->*mp);cout << endl;
	}
	return(surfInds);
}

template<class T, class R> 
R ConcatenateVectorField(const ArrayStruct<T> &arrayIn, R T::*mp, int rStart,int rEnd)
{
	R surfInds;
	int ii;
	auto itVecInt=surfInds.begin();

	for(ii=rStart; ii<rEnd;++ii){
		surfInds.insert(itVecInt, (arrayIn(ii)->*mp).begin(),
			(arrayIn(ii)->*mp).end());
		itVecInt=surfInds.end();
		//vDisplayVector(arrayIn(subList[ii])->*mp);cout << endl;
	}
	return(surfInds);
}

template<class T, class R> 
vector<R> ConcatenateScalarField(const ArrayStruct<T> &arrayIn, R T::*mp, int rStart,int rEnd)
{
	vector<R> surfInds;
	int ii;
	surfInds.reserve((rStart-rEnd)>0?(rStart-rEnd):0);
	for(ii=rStart; ii<rEnd;++ii){
		surfInds.push_back((arrayIn(ii)->*mp));
		//vDisplayVector(arrayIn(subList[ii])->*mp);cout << endl;
	}
	return(surfInds);
}

template<class T, class R, class U, class  V> 
void OperArrayStructMethod(const ArrayStruct<T> &arrayIn,const vector<int> &subList,  R T::*mp , U &out , V oper)
{
	
	int ii;

	for(ii=0; ii<int(subList.size());++ii){
		out=oper(out,(arrayIn(subList[ii])->*mp)());
		
		//vDisplayVector(arrayIn(subList[ii])->*mp);cout << endl;
	}

}
template<class T, class R> vector<R> ReturnDataEqualRange(T key,const unordered_multimap<T,R> &hashTable)
{
	vector<R> subList;
	
	subList.reserve(5);
	auto range=hashTable.equal_range(key);
	for (auto it = range.first; it != range.second; ++it) {

		subList.push_back(it->second);
	}

	return(subList);
}
template<class T, class R> void ReturnDataEqualRange(T key,const unordered_multimap<T,R> &hashTable, vector<R> &subList)
{
	
	subList.clear();
	subList.reserve(5);
	auto range=hashTable.equal_range(key);
	for (auto it = range.first; it != range.second; ++it) {

		subList.push_back(it->second);
	}

	
}
template <typename T> inline void sort(vector<T> &vec)
{
	sort(vec.begin(),vec.end());
}
template <typename T> inline void unique(vector<T> &vec)
{
	auto itVecInt = std::unique (vec.begin(), vec.end());       
	vec.resize( std::distance(vec.begin(),itVecInt));
}
// template <typename T> inline void set_intersection(vector<T> &targVec,vector<T> &vec1,vector<T> &vec2,bool isSort)
// {
// 	T tempType;
// 	typename std::vector<T>::iterator it;
// 	targVec.assign(vec1.size(),tempType);
// 	if(!isSort){
// 		sort(vec1);
// 		sort(vec2);
// 	}
// 	it =set_intersection(vec1.begin(),vec1.end(),vec2.begin(),vec2.end(),targVec.begin());

// 	targVec.resize(it-targVec.begin());
// }

template <typename T> inline void set_intersection(vector<T> &targVec,const vector<T> &vec1,const vector<T> &vec2,bool isSort)
{
	T tempType;
	typename std::vector<T>::iterator it;
	targVec.assign(vec1.size(),tempType);
	if(!isSort){
		throw invalid_argument("Constant vectors are unsorted cannot intersect");
	}
	it =set_intersection(vec1.begin(),vec1.end(),vec2.begin(),vec2.end(),targVec.begin());

	targVec.resize(it-targVec.begin());
}

// Find option for vectors


// Templated test for all types that need to be derived from ArrayStruct<T>

//

template <class T> int TestTemplate_ArrayStruct()
{
	ArrayStruct<T> stackT,stackT2,stackT3;
	T singleT;
	vector<int> testSub = {2,5,10,7};
	vector<int> delInd={1,2};
	FILE *fidw,*fidr;
	int i=0,j=0;
	int errFlag=0;
	bool errTest;

	try {
		
		// Test Initialisation
		stackT.assign(5,singleT);
		stackT.PopulateIndices();
		stackT2.Init(5);
		stackT2.PopulateIndices();

		errTest=CompareDisp(stackT,stackT2);
		if (!errTest){
			cerr << "Error Displays were not the same (Stage 1)" << endl;
			errFlag++;
		} 
		// Test ASsignement
		singleT.index=10;
		stackT[2]=singleT;
		errTest=CompareDisp(stackT,stackT2);
		if (errTest){
			cerr << "Error Displays were not the same - assignement not working (Stage 2)" << endl;
			errFlag++;
		} 
		errFlag+=TestReadyness(stackT," after assignement",false);

		// Test prepare
		stackT.PrepareForUse();
		stackT2.PrepareForUse();
		errFlag+=TestReadyness(stackT," after PrepareForUse",true);

		// Test Find
		i=stackT.find(10);
		j=stackT.find(7);
		testSub=stackT.find_list(testSub);
		errFlag+=TestReadyness(stackT," after find",true);


		if (i!=2 || j!=-1){
			std::cerr << "FIND did not Succesfully identify the indices" << std::endl;
			errFlag++;
		}
		if (testSub[2]!=2 || testSub[3]!=-1){
			std::cerr << "FIND_LISTS did not Succesfully identify the indices" << std::endl;
			errFlag++;
		}

		try{
			stackT[2].index=7;
			errFlag+=TestReadyness(stackT," after [] assignement",false);

			i=stackT.find(10);
			// If this does not throw an error with have an issue
			#ifdef SAFE_ACCESS
			std::cerr << "FIND did not throw an error at the unsafe access" << std::endl;
			errFlag++;
			#endif //SAFE_ACCESS
		}  catch (exception const& ex) { 
			cout << "Previous Call should have thrown the following warning:" << endl;
			cout << "Warning: reading from potentially obsolete unordered_map" << endl;
			i=-1;
		}
		// Test Concatenation of arrays

		stackT.PrepareForUse();
		errFlag+=TestReadyness(stackT," after PrepareForUse (2nd Time)",true);

		stackT.Concatenate(stackT2);
		errFlag+=TestReadyness(stackT," after Concatenate",false);

		stackT.PrepareForUse();
		errFlag+=TestReadyness(stackT," after PrepareForUse (3nd Time)",true);

		stackT.PopulateIndices();
		errFlag+=TestReadyness(stackT," after PopulateIndices",false);

		stackT.PrepareForUse();
		errFlag+=TestReadyness(stackT," after PrepareForUse (4th Time)",true);

		// Test Read write:
		fidw=fopen("../TESTOUT/testarray.dat","w");
		if(fidw!=NULL){
			stackT.write(fidw);
			errFlag+=TestReadyness(stackT," Write out",true);
			fclose(fidw);
			fidr=fopen("../TESTOUT/testarray.dat","r");
			if(fidr!=NULL){
				stackT3.read(fidr);
				fclose(fidr);
				stackT3.PrepareForUse();
				errTest=CompareDisp(stackT,stackT3);
				if (!errTest){
					cerr << "Error Displays were not the same after write read" << endl;
					errFlag++;
				}
			}else{
				cout << "Error: Could not open file to read test arraystructures." << endl;
				errFlag++;  
			}
		}else{
			cout << "Error: Could not open file to write test arraystructures." << endl;
			errFlag++;  
		}


		stackT2=stackT;
		errTest=CompareDisp(stackT,stackT2);
		if (!errTest){
			cerr << "Error Displays were not the same after full assignement" << endl;
			errFlag++;
		}
		stackT.disp();
		errFlag+=TestReadyness(stackT," Disp",true);
		
		stackT.elems.erase(stackT.elems.begin(),++(++stackT.elems.begin()));
		stackT.isHash=0;
		stackT.isSetMI=0;
		stackT2.remove(delInd);

		errTest=CompareDisp(stackT,stackT2);
		if (!errTest){
			stackT.disp();
			stackT2.disp();
			cerr << "Error Displays were not the same erase" << endl;
			errFlag++;
		}

	} catch (exception const& ex) { 
		cerr << "Exception: " << ex.what() <<endl; 
		return -1;
	} 
	return(errFlag);
}

template<class T> int TestReadyness(T &stackT, const char* txt, bool errTarg)
{
	// Test if the behaviour of the readyness test is as expected
	// stackT is the class to be tested 
	// txt is a string to be displayed telling the user where the error is occuring
	// errTarg is 
	bool errTest;
	int errFlag=0;
	errTest=stackT.isready();
	if (!(errTest==errTarg)){
		cerr << "stackT wrongly marked as " << (errTarg ? "not" : "" ) <<" ready (isready()) " << txt << endl;
		errFlag++;
	} 
	errTest=stackT.checkready();
	if (!(errTest==errTarg)){
		cerr << "stackT wrongly marked as " << (errTarg ? "not" : "" ) << " ready (checkready())"  << txt  << endl;
		errFlag++;
	} 
	return(errFlag);
}

// member function definition template <class T> : "ArrayStruct"

template<class T> inline int ArrayStruct <T>::GetMaxIndex() const 
{
	#ifdef SAFE_ACCESS
	if (!isSetMI)
	{
		cerr << "warning: Potentially unsafe reading of max index - execute PrepareForUse() before access" << endl;
		cerr << "          in " << __PRETTY_FUNCTION__ << endl; 
	}
	#endif //SAFE_ACCESS
	return(maxIndex);
}

template<class T>  int ArrayStruct <T>::find(int key, bool noWarn) const 
{
	if (isHash==0 && !noWarn){
		cerr << "Warning: reading from potentially obsolete unordered_map " << endl;
		cerr << "          in " << __PRETTY_FUNCTION__ << endl; 
		cerr << "          To avoid this message perform read operations on ArrayStruct<T> using the () operator" << endl; 
	}
	auto search=hashTable.find(key);
	

	if (search==hashTable.end()){
		return(-1);
	}
	#ifdef SAFE_ACCESS
	int key2;
	key2=elems[search->second].index;
	if (key2!=key){
		cerr << "          Error in " << __PRETTY_FUNCTION__ << endl; 
		throw  invalid_argument ("FIND returned an invalid output ");
	}
	#endif //SAFE_ACCESS
	return(search->second);
}

template<class T> vector<int> ArrayStruct <T>::find_list(const vector<int> &key) const 
{
	vector<int> returnSub=key;
	int ii;
	for (ii=0;ii<int(key.size());++ii){
		returnSub[ii]=this->find(key[ii]);
	}
	return(returnSub);
}

template<class T> void ArrayStruct <T>::write(FILE *fid) const 
{
	int ii;
	fprintf(fid,"%i %i \n",int(this->size()),int(isInMesh));
	for (ii=0;ii<int(this->size());++ii){
		elems[ii].write(fid);
	}
}
template<class T> void ArrayStruct <T>::read(FILE *fid) 
{
	int ii,n;
	fscanf(fid,"%i %i ",&n,&ii);
	isInMesh=bool(ii);
	this->Init(n);
	for (ii=0;ii<n;++ii){
		elems[ii].read(fid);
	}
}

template<class T>  bool ArrayStruct <T>::checkready()  
{	
	readyforuse=true;
	int i=0;
	readyforuse=((isHash==1) & (isSetMI==1));
	if (readyforuse){
		while ( (i < this->size()) & readyforuse)
		{
			readyforuse=readyforuse & elems[i].isready(isInMesh);
			++i;
		}
	}

	return(readyforuse);
	
}
template<class T>  void ArrayStruct <T>::ForceArrayReady()
{
	isHash=1;
	isSetMI=1;
	readyforuse=true;
}

template<class T> void ArrayStruct <T>::disp() const 
{
	cout << "Array of size " << this->size() << endl;
	for (int ii=0 ; unsigned_int(ii)<elems.size();ii++){
		cout << "Array " << ii << " " ;
		elems[ii].disp();
	}
	cout << "Array Dat: isHash " << isHash << "; isSetMI " << isSetMI << "; isInMesh "  << isInMesh << endl;
}
template<class T> void ArrayStruct <T>::disp(int iStart, int iEnd) const 
{
	cout << "Array of size " << this->size() << endl;
	for (int ii=iStart ; iEnd<elems.size();ii++){
		cout << "Array " << ii << " " ;
		elems[ii].disp();
	}
	cout << "Array Dat: isHash " << isHash << "; isSetMI " << isSetMI << "; isInMesh "  << isInMesh << endl;
}
template<class T> void ArrayStruct <T>::disp(const vector<int> &subs) const
{
	int ii;
	cout << "Array of size " << this->size() << " displaying subset of size " << subs.size() << endl;
	for (ii=0 ; unsigned_int(ii)<subs.size();ii++){
		cout << "Array " << subs[ii] << " " ;
		elems[subs[ii]].disp();
	}
}

template<class T> inline void ArrayStruct <T>::Init(int n)
{
	T sT;
	this->clear();
	this->reserve(n);
	this->assign(n,sT);
}


template<class T> void ArrayStruct <T>::HashArray()
{
	// Generates a valid unordered_map for the current ArrayStruct
	// Function should not be called repeatedly 
	if(!hashTable.empty()){
		hashTable.clear();
	}
	hashTable.reserve(elems.size());
	for (int i = 0; i < int(elems.size()); ++i)
	{
		hashTable.emplace(elems[i].Key(),i);
	}
	isHash=1;
	
	
   //cout << "Array Struct Succesfully Hashed" << endl;
}

template<class T> void ArrayStruct <T>::TightenConnectivity()
{


	for (int i = 0; i < int(elems.size()); ++i)
	{
		elems[i].TightenConnectivity();
	}
	
	
   //cout << "Array Struct Succesfully Hashed" << endl;
}

template<class T> void ArrayStruct <T>::SetMaxIndex()
{
	// Sets the correct value of maxIndex
	int n=elems.size();
	maxIndex=0;
	for(int ii=n-1;-1<ii;--ii){
		if (maxIndex<this->elems[ii].index){
			maxIndex=elems[ii].index;
		}
	}
	isSetMI=1;
}
/*template<class T> void ArrayStruct <T>::SetMaxIndex()
{
	// Sets the correct value of maxIndex
	int n=elems.size();
	maxIndex=0;
	for(int ii=0;ii<n;++ii){
		maxIndex= (maxIndex<this->elems[ii].index) ? elems[ii].index : maxIndex;
	}
	isSetMI=1;
}*/
template<class T> void ArrayStruct <T>::Concatenate(const ArrayStruct<T> &other)
{
	int nCurr,nNew,nTot,ii;
	nCurr=this->size();
	nNew=other.size();
	nTot=nCurr+nNew;
	elems.reserve(nTot);

	for (ii=0; ii<nNew;ii++){
		elems.push_back(other.elems[ii]);
	}
	//this->HashArray();
	//this->SetMaxIndex();
	isHash=0;
	isSetMI=0;
	readyforuse=false;
}


template<class T> void ArrayStruct <T>::PopulateIndices()
{
	int n=elems.size();
	maxIndex=-1;
	for(int ii=0;ii<n;++ii){
		elems[ii].index=ii+1;
	}
	maxIndex=n;
	isSetMI=1;
	isHash=0;
	readyforuse=false;
}
template<class T> void ArrayStruct <T>::ChangeIndices(int nVert,int nEdge,int nSurf,int nVolu)
{
	int n=elems.size();
	for(int ii=0;ii<n;++ii){
		elems[ii].ChangeIndices( nVert, nEdge, nSurf, nVolu);
	}

	isSetMI=0;
	isHash=0;
	readyforuse=false;
}
template<class T> void ArrayStruct <T>::PrepareForUse()
{

	for (int ii = 0; ii < this->size(); ++ii)
	{
		this->elems[ii].PrepareForUse();
	}
	if (isSetMI==0){
		this->SetMaxIndex();
	}
	if (isHash==0){
		this->HashArray();
	}
	readyforuse=true;
	
}

template<class T> void ArrayStruct<T>::remove(vector<int> delInd)
{

	HashedVector<int,T> delHash;
	//int isHashCurr=isHash;

	delHash.vec=delInd;
	delHash.GenerateHash();

	elems.erase(std::remove_if(elems.begin(),elems.end(),delHash),elems.end());

	isHash=0;
	isSetMI=0;
	// This part of the code lets remove maintain the hashTable up to date.
	/*for (int i = 0; i < int(delInd.size()); ++i)
	{
		this->hashTable.erase(delInd[i]);
	}
	isHash=isHashCurr;*/

}

// Implementation of vector member functions into the base class

template<class T> inline int ArrayStruct <T>::size() const 
{
	return(int(elems.size()));
}
template<class T> inline int ArrayStruct <T>::capacity() const 
{
	return(int(elems.capacity()));
}
template<class T> inline void ArrayStruct <T>::assign(int n, T& newelem)  
{
	elems.assign(n,newelem);
	isHash=0;
	isSetMI=0;
	readyforuse=false;
}
template<class T> inline void ArrayStruct <T>::push_back(T& newelem)  
{
	elems.push_back(newelem);

	maxIndex=newelem.index>maxIndex ? newelem.index : maxIndex;
	hashTable.emplace(newelem.Key(),elems.size()-1);

}
template<class T> inline void ArrayStruct <T>::reserve(int n)  
{
	elems.reserve(n);
}
template<class T> inline void ArrayStruct <T>::clear()  
{
	elems.clear();
	hashTable.clear();
	// isHash=0;
	maxIndex=0;
	// readyforuse=false;
}
// Hashed Vector Template class implementations

template <class T,class Q,class R> inline void HashedVector<T,Q,R>::GenerateHash()
{
	hashTable.clear();
	HashVector(vec, hashTable);
	isHash=true;
}
template <class T,class Q,class R> inline int HashedVector<T,Q,R>::find(const T key) const
{
	return(FindSub(key, hashTable) );
}
template <class T,class Q,class R> inline vector<int> HashedVector<T,Q,R>::findall(const T key) const
{
	return(ReturnDataEqualRange(key,hashTable));
}

template <class T,class Q,class R> inline int HashedVector<T,Q,R>::count(const T key) const
{
	return(hashTable.count(key) );
}
template <class T,class Q,class R> vector<int> HashedVector<T,Q,R>::count(const vector<T> &key) const
{
	vector<int> subOut;
	subOut.reserve(key.size());
	for (int i = 0; i < int(key.size()); ++i)
	{
		subOut.push_back(hashTable.count(key[i]));
	}

	return(subOut);
}
template <class T,class Q,class R> inline vector<int> HashedVector<T,Q,R>::find_list(const vector<T> &key) const
{
	return(FindSubList(key, vec, hashTable) );
}
template <class T,class Q,class R> inline bool HashedVector<T,Q,R>::IsInVec(const Q &key) const
{
	return(FindSub(key.Key(), hashTable)!=-1);
}
template <class T,class Q,class R>  bool HashedVector<T,Q,R>::operator()(const Q &key) const
{
	return(FindSub(key.Key(), hashTable)!=-1);
}

template <class T,class Q,class R> inline void HashedMap<T,Q,R>::GenerateHash()
{

	hashTable.clear();
	HashVector(vec, hashTable, targ);
	isHash=true;
}
template<class T>  int FindSub(const T &key, const unordered_multimap<T,int> &hashTable)  
{
	auto search=hashTable.find(key);

	if (search==hashTable.end()){
		return(-1);
	}

	return(search->second);
}

template<class T> vector<int> FindSubList(const vector<T> &keyFind,const vector<T> &keyList, unordered_multimap<T,int> &hashTable) 
{
	vector<int> returnSub;
	int ii;

	returnSub.reserve(int(keyFind.size()));

	if (hashTable.empty()){
		HashVector(keyList,hashTable);
	}

	for (ii=0;ii<int(keyFind.size());++ii){
		returnSub.push_back(FindSub(keyFind[ii],hashTable));
      #ifdef SAFE_ACCESS
		if(returnSub[ii]>=0){
			if (keyList[returnSub[ii]]!=keyFind[ii]){
				throw  invalid_argument ("FIND returned an invalid output ");
			}
		}
      #endif //SAFE_ACCESS
	}
	return(returnSub);
}


template<class T> vector<int> FindSubList(const vector<T> &keyFind,const vector<T> &keyList, const unordered_multimap<T,int> &hashTable) 
{
	vector<int> returnSub;
	int ii;

	returnSub.reserve(int(keyFind.size()));

	if (hashTable.empty() && keyList.size()>0 && keyFind.size()>0){
		throw invalid_argument("Hash table needs to be full to find list");
	}

	for (ii=0;ii<int(keyFind.size());++ii){
		returnSub.push_back(FindSub(keyFind[ii],hashTable));
      #ifdef SAFE_ACCESS
		if(returnSub[ii]>=0){
			if (keyList[returnSub[ii]]!=keyFind[ii]){
				throw  invalid_argument ("FIND returned an invalid output ");
			}
		}
      #endif //SAFE_ACCESS
	}
	return(returnSub);
}


template<class T, class Q> void HashVector(const vector<T> &elems, unordered_multimap<T,Q> &hashTable,
	 const vector<Q> &targElems)
{
   // Generates a valid unordered_map for the current ArrayStruct
   // Function should not be called repeatedly 

	hashTable.reserve(elems.size());
	if(targElems.size()==0){
		for (int i = 0; i < int(elems.size()); ++i)
		{
			hashTable.emplace(elems[i],i);
		}
	} else if (targElems.size()==elems.size()){
		for (int i = 0; i < int(elems.size()); ++i)
		{
			hashTable.emplace(elems[i],targElems[i]);
		}
	} else {
		throw invalid_argument("HashVector failed as input vectors are of different sizes");
	}

}



template<template<class Q, class R> class T,class Q, class R>
	void EraseKeyPair(T<Q,R> hashTable, Q key, R pos){

	typename T<Q,R>::iterator it ;
	it= hashTable.find(key);
	while(it->second!=pos && it->first==key){++it;}

	if (it->second==pos && it->first==key){
		hashTable.erase (it);
	} else {

		cerr << "Error: Key value pair not found and could not be removed "<< endl;
		cerr << " key " << key << " pos " << pos << endl;
		cerr << "	in function:" <<  __PRETTY_FUNCTION__ << endl;
		throw invalid_argument ("Key value pair not found and could not be removed");
	}
}


#endif // ARRAYSTRUCTS_INCL_H_INCLUDED 