#ifndef RSVSINTEGRATION_H_INCLUDED 
#define RSVSINTEGRATION_H_INCLUDED 


//=================================
// forward declared dependencies
// 		class foo; //when you only need a pointer not the actual object
// 		and to avoid circular dependencies


//=================================
// included dependencies
#include <vector>

#include "snake.hpp"

// ================================
// declarations

using namespace std; 
 
void SnakeConnectivityUpdate(snake &testSnake,  vector<int> &isImpact);
void SnakeConnectivityUpdate_legacy(snake &snakein,  vector<int> &isImpact);
void SnakeConnectivityUpdate_robust(snake &snakein,  vector<int> &isImpact);
int TimeStamp(const char* str,int start_s);

//==================================
// Code
// NOTE: function in a class definition are IMPLICITELY INLINED 
//       ie replaced by their code at compile time
	

#endif //RSVSINTEGRATION_H_INCLUDED 