#include <iostream>
#include <cmath>
#include <ctime>
#include <vector>

#include "RSVSintegration.hpp"
#include "snake.hpp"
#include "snakeengine.hpp"

int SAFE_ALGO_TestConn(snake &snakein){
	int ret=0;

	if (snakein.Check3D()){
		#ifdef SAFE_ALGO
		ret = snakein.snakeconn.TestConnectivityBiDir();
		#endif //SAFE_ALGO
	}

	return(ret);
}

void SnakeConnectivityUpdate_legacy(snake &snakein,  vector<int> &isImpact){


	int start_s;

	start_s=clock();

	

	start_s=TimeStamp("position: ", start_s);

	snakein.SnaxImpactDetection(isImpact);
	MergeAllContactVertices(snakein, isImpact);
	snakein.PrepareForUse();


	start_s=TimeStamp("Merge: ", start_s);

	CleanupSnakeConnec(snakein);
	

	start_s=TimeStamp("Clean: ", start_s);
	snakein.SnaxImpactDetection(isImpact);
	SpawnArrivedSnaxels(snakein,isImpact);


	

	start_s=TimeStamp("Spawn: ", start_s);

	snakein.SnaxImpactDetection(isImpact);
	snakein.SnaxAlmostImpactDetection(isImpact, 0.01);
	snakein.PrepareForUse();
	SAFE_ALGO_TestConn(snakein);
	MergeAllContactVertices(snakein, isImpact);
	snakein.PrepareForUse();
	SAFE_ALGO_TestConn(snakein);
	

	start_s=TimeStamp("Impact: ", start_s);

	CleanupSnakeConnec(snakein);

	SAFE_ALGO_TestConn(snakein);
	snakein.OrientFaces();

	start_s=TimeStamp("Clean: ", start_s);

	

}

void SnakeConnectivityUpdate_robust(snake &snakein,  vector<int> &isImpact){
	/*
	Performs the snake step except the movement of the snake.

	This one performs it in two steps:
	 1) Impact Merge Clean
	 2) Impact Spawn Impact Merge Clean

	This function might be better in snakeengine.cpp
	*/
	double impactAlmostRange = 0.2;

	int start_s, start_f;
	start_f=clock();

	//===============================
	// Impact on edge
	// ======================
	// Impact
	SAFE_ALGO_TestConn(snakein);
	snakein.SnaxImpactDetection(isImpact);
	snakein.PrepareForUse();
	start_s=TimeStamp("Impact: ", start_f);
	// ======================
	// Merge
	MergeAllContactVertices(snakein, isImpact);
	snakein.PrepareForUse();
	start_s=TimeStamp("Merge: ", start_s);
	// ======================
	// Clean
	CleanupSnakeConnec(snakein);
	snakein.PrepareForUse();
	SAFE_ALGO_TestConn(snakein);
	snakein.OrientFaces();
	start_s=TimeStamp("Clean: ", start_s);


	//===============================
	// Spawn
	// ======================
	// Impact
	SAFE_ALGO_TestConn(snakein);
	snakein.SnaxImpactDetection(isImpact);
	snakein.SnaxAlmostImpactDetection(isImpact, impactAlmostRange);
	snakein.PrepareForUse();
	start_s=TimeStamp("Impact: ", start_s);
	// ======================
	// Spawn
	SAFE_ALGO_TestConn(snakein);
	snakein.SnaxImpactDetection(isImpact);
	SpawnArrivedSnaxels(snakein,isImpact);
	start_s=TimeStamp("Spawn: ", start_s);
	// ======================
	// Impact
	SAFE_ALGO_TestConn(snakein);
	snakein.SnaxImpactDetection(isImpact);
	snakein.PrepareForUse();
	start_s=TimeStamp("Impact: ", start_s);
	// ======================
	// Merge
	MergeAllContactVertices(snakein, isImpact);
	snakein.PrepareForUse();
	start_s=TimeStamp("Merge: ", start_s);
	// ======================
	// Clean
	CleanupSnakeConnec(snakein);
	snakein.PrepareForUse();
	SAFE_ALGO_TestConn(snakein);
	snakein.OrientFaces();
	start_s=TimeStamp("Clean: ", start_s);

	TimeStamp(" - Connec Update: ", start_f);
	
}

void SnakeConnectivityUpdate(snake &snakein,  vector<int> &isImpact){
	/*
	Performs the snake step except the movement of the snake.
	This one performs it in a 'single' step:
	 Impact Spawn Impact Merge Clean

	This function might be better in snakeengine.cpp
	*/
	double impactAlmostRange = 0.1;

	int start_s, start_f;
	start_f=clock();


	//===============================
	// Spawn
	// ======================
	// Impact
	SAFE_ALGO_TestConn(snakein);
	snakein.SnaxImpactDetection(isImpact);
	snakein.SnaxAlmostImpactDetection(isImpact, impactAlmostRange);
	snakein.PrepareForUse();
	start_s=TimeStamp("Impact: ", start_f);
	// ======================
	// Spawn
	SAFE_ALGO_TestConn(snakein);
	snakein.SnaxImpactDetection(isImpact);
	SpawnArrivedSnaxels(snakein,isImpact);
	start_s=TimeStamp("Spawn: ", start_s);
	// ======================
	// Impact
	SAFE_ALGO_TestConn(snakein);
	snakein.SnaxImpactDetection(isImpact);
	snakein.PrepareForUse();
	start_s=TimeStamp("Impact: ", start_s);
	// ======================
	// Merge
	MergeAllContactVertices(snakein, isImpact);
	snakein.PrepareForUse();
	start_s=TimeStamp("Merge: ", start_s);
	// ======================
	// Clean
	CleanupSnakeConnec(snakein);
	snakein.PrepareForUse();
	SAFE_ALGO_TestConn(snakein);
	snakein.OrientFaces();
	start_s=TimeStamp("Clean: ", start_s);

	TimeStamp(" - Connec Update: ", start_f);
	
}


int TimeStamp(const char* str,int start_s){
	int stop_s=clock();
	#ifdef TIME_EXEC
	cout << str << " " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << "ms; ";
	#endif
	return(stop_s);

}