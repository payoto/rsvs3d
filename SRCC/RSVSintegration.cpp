#include <iostream>
#include <cmath>
#include <ctime>
#include <vector>

#include "RSVSintegration.hpp"
#include "snake.hpp"
#include "snakeengine.hpp"

void SnakeConnectivityUpdate_legacy(snake &snakein,  vector<int> &isImpact){

	#ifdef TIME_EXEC
	int start_s;

	start_s=clock();
	#endif
	
	#ifdef TIME_EXEC
	start_s=TimeStamp("position: ", start_s);
	#endif //TIME_EXEC

	snakein.SnaxImpactDetection(isImpact);
	MergeAllContactVertices(snakein, isImpact);
	snakein.PrepareForUse();

	#ifdef TIME_EXEC
	start_s=TimeStamp("Merge: ", start_s);
	#endif //TIME_EXEC

	CleanupSnakeConnec(snakein);
	#ifdef SAFE_ALGO
	if (snakein.Check3D()){
		snakein.snakeconn.TestConnectivityBiDir();
	}
	#endif
	#ifdef TIME_EXEC
	start_s=TimeStamp("Clean: ", start_s);
	#endif //TIME_EXEC
	snakein.SnaxImpactDetection(isImpact);
	SpawnArrivedSnaxels(snakein,isImpact);


	
	#ifdef TIME_EXEC
	start_s=TimeStamp("Spawn: ", start_s);
	#endif //TIME_EXEC

	snakein.SnaxImpactDetection(isImpact);
	snakein.SnaxAlmostImpactDetection(isImpact, 0.01);
	snakein.PrepareForUse();
	#ifdef SAFE_ALGO
	if (snakein.Check3D()){
		snakein.snakeconn.TestConnectivityBiDir();
	}
	#endif
	MergeAllContactVertices(snakein, isImpact);
	snakein.PrepareForUse();
	#ifdef SAFE_ALGO
	if (snakein.Check3D()){
		snakein.snakeconn.TestConnectivityBiDir();
	}
	#endif
	
	#ifdef TIME_EXEC
	start_s=TimeStamp("Impact: ", start_s);
	#endif //TIME_EXEC

	CleanupSnakeConnec(snakein);

	#ifdef SAFE_ALGO
	if (snakein.Check3D()){
		snakein.snakeconn.TestConnectivityBiDir();
	}
	#endif
	snakein.OrientSurfaceVolume();
	#ifdef TIME_EXEC
	start_s=TimeStamp("Clean: ", start_s);
	#endif //TIME_EXEC

	

}

void SnakeConnectivityUpdate(snake &snakein,  vector<int> &isImpact){

	double impactAlmostRange = 0.01;
	#ifdef TIME_EXEC
	int start_s, start_f;
	start_f=clock();
	#endif
	// ======================
	// Spawn
	#ifdef SAFE_ALGO
	if (snakein.Check3D()){
		snakein.snakeconn.TestConnectivityBiDir();
	}
	#endif

	snakein.SnaxImpactDetection(isImpact);
	SpawnArrivedSnaxels(snakein,isImpact);
	
	#ifdef TIME_EXEC
	start_s=TimeStamp("Spawn: ", start_f);
	#endif //TIME_EXEC

	// ======================
	// Impact
	#ifdef SAFE_ALGO
	if (snakein.Check3D()){
		snakein.snakeconn.TestConnectivityBiDir();
	}
	#endif

	snakein.SnaxImpactDetection(isImpact);
	snakein.SnaxAlmostImpactDetection(isImpact, impactAlmostRange);
	snakein.PrepareForUse();

	#ifdef TIME_EXEC
	start_s=TimeStamp("Impact: ", start_s);
	#endif //TIME_EXEC

	// ======================
	// Merge

	MergeAllContactVertices(snakein, isImpact);
	snakein.PrepareForUse();
	

	#ifdef TIME_EXEC
	start_s=TimeStamp("Merge: ", start_s);
	#endif // TIME_EXEC

	// ======================
	// Clean


	CleanupSnakeConnec(snakein);

	#ifdef SAFE_ALGO
	if (snakein.Check3D()){
		snakein.snakeconn.TestConnectivityBiDir();
	}
	#endif
	snakein.PrepareForUse();
	snakein.OrientSurfaceVolume();

	#ifdef TIME_EXEC
	start_s=TimeStamp("Clean: ", start_s);
	TimeStamp(" - Connec Update: ", start_f);
	#endif // TIME_EXEC
}


int TimeStamp(const char* str,int start_s){
	int stop_s=clock();
	cout << str << " " << (stop_s-start_s)/double(CLOCKS_PER_SEC)*1000 << "ms; ";

	return(stop_s);

}