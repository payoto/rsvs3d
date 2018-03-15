/*
File for the implementation of the class template ArrayStruct
this .cpp file is INCLUDED as part of arraystructures.hpp and  cannot be 
compiled on its own.

*/

#ifndef ARRAYSTRUCTS_INCL_H_INCLUDED
#define ARRAYSTRUCTS_INCL_H_INCLUDED

template <class T> bool CompareDisp(T &mesh1,T &mesh2){
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


// Templated test for all types that need to be derived from ArrayStruct<T>

//

template <class T> int TestTemplate_ArrayStruct()
{
	ArrayStruct<T> stackT,stackT2;
	T singleT;
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

		// Test prepare
		stackT.PrepareForUse();
		stackT2.PrepareForUse();

		// Test Find
		i=stackT.find(10);
		j=stackT.find(7);

		if (i!=2 || j!=-1){
			cerr << "FIND did not Succesfully identify the indices" << endl;
			errFlag++;
		}

		try{
			stackT[2].index=7;
			i=stackT.find(10);
			// If this does not throw an error with have an issue
			#ifdef SAFE_ACCESS
			cerr << "FIND did not throw an error at the unsafe access" << endl;
			errFlag++;
			#endif //SAFE_ACCESS
		}  catch (exception const& ex) { }
		// Test Concatenation of arrays

		stackT.Concatenate(stackT2);
		stackT.PopulateIndices();
		stackT.PrepareForUse();
		stackT.disp();

	} catch (exception const& ex) { 
		cerr << "Exception: " << ex.what() <<endl; 
		return -1;
	} 
	return(errFlag);
}


// member function definition template <class T> : "ArrayStruct"

template<class T> inline int ArrayStruct <T>::GetMaxIndex() const 
{
	return(maxIndex);
}

template<class T> inline int ArrayStruct <T>::find(int key) const 
{
	if (isHash==0){
		cerr << "Warning: reading from potentially obsolete unordered_map " << endl;
		cerr << "          in ArrayStruct <T>::find(int key)" << endl; 
		cerr << "          To avoid this message perform read operations on ArrayStruct<T> using the () operator" << endl; 
	}
	auto search=hashTable.find(key);
	if (search==hashTable.end()){
		return(-1);
	}
	#ifdef SAFE_ACCESS
	if (elems[search->second].index!=key){
		throw  invalid_argument ("FIND returned an invalid output ");
	}
	#endif //SAFE_ACCESS
	return(search->second);
}




template<class T> void ArrayStruct <T>::disp() const 
{

	for (int ii=0 ; unsigned_int(ii)<elems.size();ii++){
		cout << "Array " << ii << " " ;
		elems[ii].disp();
	}
	cout << "Array Dat: isHash " << isHash << "; isSetMI " << isSetMI << endl;
}


template<class T> inline void ArrayStruct <T>::Init(int n)
{
	T sT;
	elems.reserve(n);
	elems.assign(n,sT);
}


template<class T> void ArrayStruct <T>::HashArray()
{
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
	//this->HashArray();
	//this->SetMaxIndex();
	isHash=0;
	isSetMI=0;
}

template<class T> void ArrayStruct <T>::SetMaxIndex()
{
	int n=elems.size();
	maxIndex=-1;
	for(int ii=0;ii<n;++ii){
		maxIndex= (maxIndex<this->elems[ii].index) ? elems[ii].index : maxIndex;
	}
	isSetMI=1;
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
}
template<class T> inline void ArrayStruct <T>::push_back(T& newelem)  
{
	elems.push_back(newelem);
	isHash=0;
	isSetMI=0;
}
template<class T> inline void ArrayStruct <T>::reserve(int n)  
{
	elems.reserve(n);
}


#endif // ARRAYSTRUCTS_INCL_H_INCLUDED 