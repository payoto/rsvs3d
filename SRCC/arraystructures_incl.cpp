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
		errFlag+=TestReadyness(stackT," after assignement",false);

		// Test prepare
		stackT.PrepareForUse();
		stackT2.PrepareForUse();
		errFlag+=TestReadyness(stackT," after PrepareForUse",true);

		// Test Find
		i=stackT.find(10);
		j=stackT.find(7);
		errFlag+=TestReadyness(stackT," after find",true);


		if (i!=2 || j!=-1){
			std::cerr << "FIND did not Succesfully identify the indices" << std::endl;
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
		
		stackT.disp();
		errFlag+=TestReadyness(stackT," Disp",true);

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
	int key2;

	if (search==hashTable.end()){
		return(-1);
	}
	#ifdef SAFE_ACCESS
	key2=elems[search->second].index;
	if (key2!=key){
		throw  invalid_argument ("FIND returned an invalid output ");
	}
	#endif //SAFE_ACCESS
	return(search->second);
}

template<class T>  bool ArrayStruct <T>::checkready()  
{	
	readyforuse=true;
	int i;
	readyforuse=((isHash==1) & (isSetMI==1));
	if (readyforuse){
		while ( (i < this->size()) & readyforuse)
		{
			readyforuse=readyforuse & elems[i].isready();
			++i;
		}
	}

	return(readyforuse);
	
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
	this->reserve(n);
	this->assign(n,sT);
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
	readyforuse=false;
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
	isHash=0;
	isSetMI=0;
	readyforuse=false;
}
template<class T> inline void ArrayStruct <T>::reserve(int n)  
{
	elems.reserve(n);
}


#endif // ARRAYSTRUCTS_INCL_H_INCLUDED 