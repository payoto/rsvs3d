/*
File for the implementation of the class template ArrayStruct
this .cpp file is INCLUDED as part of arraystructures.hpp and  cannot be 
compiled on its own.

*/
#ifndef ARRAYSTRUCTS_INCL_H_INCLUDED
#define ARRAYSTRUCTS_INCL_H_INCLUDED

#include "arraystructures.hpp"
#include "warning.hpp"

template<class T> void ConcatenateVector(std::vector<T> &vecRoot,
	const std::vector<T> &vecConcat)
{
	vecRoot.insert(vecRoot.end(),vecConcat.begin(),vecConcat.end());
}

template <class T> bool CompareDisp(T &mesh1,T &mesh2)
{
	bool compFlag;
	std::stringstream ss1,ss2;
	auto old_buf = std::cout.rdbuf(ss1.rdbuf()); 

	mesh1.disp();
	std::cout.rdbuf(ss2.rdbuf()); 
	mesh2.disp();
	std::cout.rdbuf(old_buf);

	compFlag=ss1.str().compare(ss2.str())==0;
	return(compFlag);
}

template<class T> void DisplayVector(std::vector<T> vec)
{
	PrintVector(vec, std::cout);
}
template<class T> void PrintVector(std::vector<T> vec, std::ostream &streamout)
{
	streamout << int(vec.size()) << " - "; 
	for (int i = 0; i < int(vec.size()); ++i)
	{
		streamout << vec[i] << " ";
	}
	streamout << " | " ;
}

template<class T> void DisplayVectorStatistics(std::vector<T> vec)
{
	std::cout << int(vec.size()) << " - "; 
	std::cout << *min_element(vec.begin(), vec.end()) 
		<< " " << *max_element(vec.begin(), vec.end()) << " " ;
	T s = 0;
	for (auto& n : vec){
    	s += n;
	}
	std::cout << s << " " << double(s)/double(vec.size());
	std::cout << " | " ;
}

template<class T, class R> 
R ConcatenateVectorField(const ArrayStruct<T> &arrayIn, R T::*mp,
	const std::vector<int> &subList)
{
	R surfInds;
	int ii;
	auto itVecInt=surfInds.begin();

	for(ii=0; ii<int(subList.size());++ii){
		surfInds.insert(itVecInt, (arrayIn(subList[ii])->*mp).begin(),
			(arrayIn(subList[ii])->*mp).end());
		itVecInt=surfInds.end();
		//vDisplayVector(arrayIn(subList[ii])->*mp);std::cout << std::endl;
	}
	return(surfInds);
}

template<class T, class R> 
std::vector<R> ConcatenateScalarField(const ArrayStruct<T> &arrayIn, R T::*mp,
	const std::vector<int> &subList)
{
	std::vector<R> surfInds;
	int ii;
	surfInds.reserve(subList.size());
	for(ii=0; ii<int(subList.size());++ii){
		surfInds.push_back((arrayIn(subList[ii])->*mp));
		//vDisplayVector(arrayIn(subList[ii])->*mp);std::cout << std::endl;
	}
	return(surfInds);
}

template<class T, class R> 
R ConcatenateVectorField(const ArrayStruct<T> &arrayIn, R T::*mp,
	int rStart,int rEnd)
{
	R surfInds;
	int ii;
	auto itVecInt=surfInds.begin();

	for(ii=rStart; ii<rEnd;++ii){
		surfInds.insert(itVecInt, (arrayIn(ii)->*mp).begin(),
			(arrayIn(ii)->*mp).end());
		itVecInt=surfInds.end();
		//vDisplayVector(arrayIn(subList[ii])->*mp);std::cout << std::endl;
	}
	return(surfInds);
}

template<class T, class R> 
std::vector<R> ConcatenateScalarField(const ArrayStruct<T> &arrayIn, R T::*mp,
	int rStart,int rEnd)
{
	std::vector<R> surfInds;
	int ii;
	surfInds.reserve((rStart-rEnd)>0?(rStart-rEnd):0);
	for(ii=rStart; ii<rEnd;++ii){
		surfInds.push_back((arrayIn(ii)->*mp));
		//vDisplayVector(arrayIn(subList[ii])->*mp);std::cout << std::endl;
	}
	return(surfInds);
}

template<class T, class R, class U, class  V> 
void OperArrayStructMethod(const ArrayStruct<T> &arrayIn, 
	const std::vector<int> &subList,  R T::*mp , U &out , V oper)
{
	
	int ii;

	for(ii=0; ii<int(subList.size());++ii){
		out=oper(out,(arrayIn(subList[ii])->*mp)());
		
		//vDisplayVector(arrayIn(subList[ii])->*mp);std::cout << std::endl;
	}

}
template<class T, class R> std::vector<R> ReturnDataEqualRange(T key,
	const std::unordered_multimap<T,R> &hashTable)
{
	std::vector<R> subList;
	
	subList.reserve(5);
	auto range=hashTable.equal_range(key);
	for (auto it = range.first; it != range.second; ++it) {

		subList.push_back(it->second);
	}

	return(subList);
}
template<class T, class R> void ReturnDataEqualRange(T key, 
	const std::unordered_multimap<T,R> &hashTable, std::vector<R> &subList)
{
	
	subList.clear();
	subList.reserve(5);
	auto range=hashTable.equal_range(key);
	for (auto it = range.first; it != range.second; ++it) {

		subList.push_back(it->second);
	}

	
}
template <typename T> inline void sort(std::vector<T> &vec)
{
	sort(vec.begin(),vec.end());
}
template <typename T> inline void unique(std::vector<T> &vec)
{
	auto itVecInt = std::unique (vec.begin(), vec.end());       
	vec.resize( std::distance(vec.begin(),itVecInt));
}
// template <typename T> inline void set_intersection(std::vector<T> &targVec,std::vector<T> &vec1,std::vector<T> &vec2,bool isSort)
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

template <typename T> inline void set_intersection(std::vector<T> &targVec,
	const std::vector<T> &vec1,const std::vector<T> &vec2,bool isSort)
{
	T tempType;
	typename std::vector<T>::iterator it;
	targVec.assign(vec1.size(),tempType);
	if(!isSort){
		RSVS3D_ERROR_ARGUMENT("Constant vectors are unsorted cannot intersect");
	}
	it = set_intersection(vec1.begin(), vec1.end(), vec2.begin(),
		vec2.end(), targVec.begin());

	targVec.resize(it-targVec.begin());
}

// Find option for vectors


// Templated test for all types that need to be derived from ArrayStruct<T>

//

template <class T> int TestTemplate_ArrayStruct()
{
	ArrayStruct<T> stackT,stackT2,stackT3;
	T singleT;
	std::vector<int> testSub = {2,5,10,7};
	std::vector<int> delInd={1,2};
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
			std::cerr << "Error Displays were not the same (Stage 1)" << std::endl;
			errFlag++;
		} 
		// Test ASsignement
		singleT.index=10;
		stackT[2]=singleT;
		errTest=CompareDisp(stackT,stackT2);
		if (errTest){
			std::cerr << "Error Displays were not the same - assignement not "
				"working (Stage 2)" << std::endl;
			errFlag++;
		} 
		errFlag+=TestReadiness(stackT," after assignement",false,false);

		// Test prepare
		stackT.PrepareForUse();
		stackT2.PrepareForUse();
		errFlag+=TestReadiness(stackT," after PrepareForUse",true,true);

		// Test Find
		i=stackT.find(10);
		j=stackT.find(7);
		testSub=stackT.find_list(testSub);
		errFlag+=TestReadiness(stackT," after find",true,true);


		if (i!=2 || j!=rsvs3d::constants::__notfound){
			std::cerr << "FIND did not Succesfully identify the indices"
				<< std::endl;
			errFlag++;
		}
		if (testSub[2]!=2 || testSub[3]!=rsvs3d::constants::__notfound){
			std::cerr << "FIND_LISTS did not Succesfully identify the "
				"indices" << std::endl;
			errFlag++;
		}

		try{
			stackT[2].index=7;
			errFlag+=TestReadiness(stackT," after [] assignement",false,false);

			i=stackT.find(10);
			// If this does not throw an error with have an issue
			#ifdef SAFE_ACCESS
			std::cerr << "FIND did not throw an error at the unsafe access"
				<< std::endl;
			errFlag++;
			#endif //SAFE_ACCESS
		}  catch (std::exception const& ex) { 
			std::cout << "Previous Call should have thrown the following "
				"warning:" << std::endl;
			std::cout << "Warning: reading from potentially obsolete unordered_map"
				<< std::endl;
			i=-1;
		}
		// Test Concatenation of arrays

		stackT.PrepareForUse();
		errFlag+=TestReadiness(stackT," after PrepareForUse (2nd Time Time)",true,true);

		stackT.Concatenate(stackT2);
		errFlag+=TestReadiness(stackT," after Concatenate",false,false);

		stackT.PrepareForUse();
		errFlag+=TestReadiness(stackT," after PrepareForUse (3nd Time Time)",true,true);

		stackT.PopulateIndices();
		errFlag+=TestReadiness(stackT," after PopulateIndices",false,true);

		stackT.PrepareForUse();
		errFlag+=TestReadiness(stackT," after PrepareForUse (4th Time Time)",true,true);

		// Test Read write:
		fidw=fopen("../TESTOUT/testarray.dat","w");
		if(fidw!=NULL){
			stackT.write(fidw);
			errFlag+=TestReadiness(stackT," Write out",true,true);
			fclose(fidw);
			fidr=fopen("../TESTOUT/testarray.dat","r");
			if(fidr!=NULL){
				stackT3.read(fidr);
				fclose(fidr);
				stackT3.PrepareForUse();
				errTest=CompareDisp(stackT,stackT3);
				if (!errTest){
					std::cerr << "Error Displays were not the same after write read"
						<< std::endl;
					errFlag++;
				}
			}else{
				std::cout << "Error: Could not open file to read test arraystructures."
					<< std::endl;
				errFlag++;  
			}
		}else{
			std::cout << "Error: Could not open file to write test arraystructures." 
				<< std::endl;
			errFlag++;  
		}


		stackT2=stackT;
		errTest=CompareDisp(stackT,stackT2);
		if (!errTest){
			std::cerr << "Error Displays were not the same after full assignement"
				<< std::endl;
			errFlag++;
		}
		stackT.disp();
		errFlag+=TestReadiness(stackT," Disp",true,true);
		
		stackT.elems.erase(stackT.elems.begin(),++(++stackT.elems.begin()));
		stackT.isHash=0;
		stackT.isSetMI=0;
		stackT2.remove(delInd);

		errTest=CompareDisp(stackT,stackT2);
		if (!errTest){
			stackT.disp();
			stackT2.disp();
			std::cerr << "Error Displays were not the same erase" << std::endl;
			errFlag++;
		}

	} catch (std::exception const& ex) { 
		std::cerr << "Exception: " << ex.what() <<std::endl; 
		return -1;
	} 
	return(errFlag);
}

template<class T> int TestReadiness(T &stackT, const char* txt, bool errTarg,
	bool errTargCheck)
{
	// Test if the behaviour of the readyness test is as expected
	// stackT is the class to be tested 
	// txt is a string to be displayed telling the user where the error is occuring
	// errTarg is 
	bool errTest;
	int errFlag=0;
	errTest=stackT.isready();
	if (!(errTest==errTarg)){
		std::cerr << "stackT wrongly marked as " << (errTarg ? "not" : "" ) 
			<< " ready (isready()) " << txt << std::endl;
		errFlag++;
	} 
	errTest=stackT.checkready();
	if (!(errTest==errTargCheck)){
		std::cerr << "stackT wrongly marked as " << (errTarg ? "not" : "" ) 
			<< " ready (checkready())"  << txt  << std::endl;
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
		std::cerr << "warning: Potentially unsafe reading of max index"
			" - execute PrepareForUse() before access" << std::endl;
		std::cerr << "          in " << __PRETTY_FUNCTION__ << std::endl; 
	}
	#endif //SAFE_ACCESS
	return(maxIndex);
}

template<class T>  int ArrayStruct <T>::find(int key, bool noWarn) const 
{
	if (isHash==0 && !noWarn){
		std::cerr << "Warning: reading from potentially obsolete unordered_map ";
		std::cerr << std::endl << "          in " << __PRETTY_FUNCTION__ << std::endl; 
		std::cerr << "          To avoid this message perform read operations on"
			" ArrayStruct<T> using the () operator" << std::endl; 
	}
	auto search=hashTable.find(key);
	

	if (search==hashTable.end()){
		return(rsvs3d::constants::__notfound);
	}
	#ifdef SAFE_ACCESS
	int key2;
	key2=elems[search->second].index;
	if (key2!=key){
		std::cerr << "          Error in " << __PRETTY_FUNCTION__ << std::endl; 
		RSVS3D_ERROR_ARGUMENT("FIND returned an invalid output ");
	}
	#endif //SAFE_ACCESS
	return(search->second);
}

template<class T> std::vector<int> ArrayStruct <T>::find_list(const std::vector<int> &key,
	bool noWarn) const 
{
	std::vector<int> returnSub=key;
	int ii;
	for (ii=0;ii<int(key.size());++ii){
		returnSub[ii]=this->find(key[ii], noWarn);
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
	std::cout << "Array of size " << this->size() << std::endl;
	for (int ii=0 ; unsigned_int(ii)<elems.size();ii++){
		std::cout << "Array " << ii << " " ;
		elems[ii].disp();
	}
	std::cout << "Array Dat: isHash " << isHash << "; isSetMI " << isSetMI 
		<< "; isInMesh "  << isInMesh << std::endl;
}
template<class T> void ArrayStruct <T>::disp(int iStart, int iEnd) const 
{
	std::cout << "Array of size " << this->size() << std::endl;
	for (int ii=iStart ; iEnd<elems.size();ii++){
		std::cout << "Array " << ii << " " ;
		elems[ii].disp();
	}
	std::cout << "Array Dat: isHash " << isHash << "; isSetMI " << isSetMI 
		<< "; isInMesh "  << isInMesh << std::endl;
}
template<class T> void ArrayStruct <T>::disp(const std::vector<int> &subs) const
{
	int ii;
	std::cout << "Array of size " << this->size() << " displaying subset of size " 
		<< subs.size() << std::endl;
	for (ii=0 ; unsigned_int(ii)<subs.size();ii++){
		std::cout << "Array " << subs[ii] << " " ;
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
	
	
   //std::cout << "Array Struct Succesfully Hashed" << std::endl;
}

template<class T> void ArrayStruct <T>::TightenConnectivity()
{


	for (int i = 0; i < int(elems.size()); ++i)
	{
		elems[i].TightenConnectivity();
	}
	
	
   //std::cout << "Array Struct Succesfully Hashed" << std::endl;
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
	this->HashArray();
}
template<class T> void ArrayStruct <T>::ChangeIndices(int nVert, int nEdge,
	int nSurf, int nVolu)
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

template<class T> void ArrayStruct<T>::remove(std::vector<int> delInd)
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

// Implementation of std::vector member functions into the base class

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

template <class T,class Q,class R> 
inline void HashedVector<T,Q,R>::GenerateHash()
{
	hashTable.clear();
	HashVector(vec, hashTable);
	isHash=true;
}
template <class T,class Q,class R> 
inline int HashedVector<T,Q,R>::find(const T key) const
{
	return(FindSub(key, hashTable) );
}
template <class T,class Q,class R> 
inline std::vector<int> HashedVector<T,Q,R>::findall(const T key) const
{
	return(ReturnDataEqualRange(key,hashTable));
}

template <class T,class Q,class R> 
inline int HashedVector<T,Q,R>::count(const T key) const
{
	return(hashTable.count(key) );
}

template <class T,class Q,class R> 
std::vector<int> HashedVector<T,Q,R>::count(const std::vector<T> &key) const
{
	std::vector<int> subOut;
	subOut.reserve(key.size());
	for (int i = 0; i < int(key.size()); ++i)
	{
		subOut.push_back(hashTable.count(key[i]));
	}

	return(subOut);
}

template <class T,class Q,class R> 
inline std::vector<int> HashedVector<T,Q,R>::find_list(const std::vector<T> &key) const
{
	return(FindSubList(key, vec, hashTable) );
}

template <class T,class Q,class R> 
inline bool HashedVector<T,Q,R>::IsInVec(const Q &key) const
{
	return(FindSub(key.Key(), hashTable)!=rsvs3d::constants::__notfound);
}

template <class T,class Q,class R>
bool HashedVector<T,Q,R>::operator()(const Q &key) const
{
	return(FindSub(key.Key(), hashTable)!=rsvs3d::constants::__notfound);
}

template <class T,class Q,class R> inline void HashedMap<T,Q,R>::GenerateHash()
{
	hashTable.clear();
	HashVector(vec, hashTable, targ);
	isHash=true;
}

template<class T>  int FindSub(const T &key,
	const std::unordered_multimap<T,int> &hashTable)  
{
	auto search=hashTable.find(key);

	if (search==hashTable.end()){
		return(rsvs3d::constants::__notfound);
	}

	return(search->second);
}

template<class T> std::vector<int> FindSubList(const std::vector<T> &keyFind, 
	const std::vector<T> &keyList, std::unordered_multimap<T,int> &hashTable) 
{
	std::vector<int> returnSub;
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
				RSVS3D_ERROR_ARGUMENT("FIND returned an invalid output ");
			}
		}
      #endif //SAFE_ACCESS
	}
	return(returnSub);
}


template<class T> std::vector<int> FindSubList(const std::vector<T> &keyFind, 
	const std::vector<T> &keyList, const std::unordered_multimap<T,int> &hashTable) 
{
	std::vector<int> returnSub;
	int ii;

	returnSub.reserve(int(keyFind.size()));

	if (hashTable.empty() && keyList.size()>0 && keyFind.size()>0){
		RSVS3D_ERROR_ARGUMENT("Hash table needs to be full to find list");
	}

	for (ii=0;ii<int(keyFind.size());++ii){
		returnSub.push_back(FindSub(keyFind[ii],hashTable));
      #ifdef SAFE_ACCESS
		if(returnSub[ii]>=0){
			if (keyList[returnSub[ii]]!=keyFind[ii]){
				RSVS3D_ERROR_ARGUMENT("FIND returned an invalid output ");
			}
		}
      #endif //SAFE_ACCESS
	}
	return(returnSub);
}


template<class T, class Q> void HashVector(const std::vector<T> &elems, 
	std::unordered_multimap<T,Q> &hashTable, const std::vector<Q> &targElems)
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
		RSVS3D_ERROR_ARGUMENT("HashVector failed as input vectors are of "
			"different sizes");
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
		std::cerr << "Error: Key value std::pair not found and could not be removed "<< std::endl;
		std::cerr << " key " << key << " pos " << pos << std::endl;
		std::cerr << "	in function:" <<  __PRETTY_FUNCTION__ << std::endl;
		RSVS3D_ERROR_ARGUMENT("Key value std::pair not found and could not be removed");
	}
}


#endif // ARRAYSTRUCTS_INCL_H_INCLUDED 