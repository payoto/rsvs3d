#include <stdio.h>
#include <stdlib.h>
#ifndef VOXEL_H_INCLUDED
#define VOXEL_H_INCLUDED
/* Type Definitions */
typedef struct {
	int index;
	double fill;
	int refineLvl;
	int *refineVec; /* This has to be allocated  */
	
} cellTemplate;

typedef struct {
	int index;
	int cellind[2];
	int vertex[2];
	int orientation; /* 1 is vertical 0 is Horizontal */
} edgeTemplate;

typedef struct {
	int index;
	double coord[2]; /* needs to be changed wiht dim */
	
} vertexTemplate;

/* Global Extern Variables */
/*extern int 2; */
/* Macros */
#if defined(__GNUC__) || defined(__GNUG__)
	/* GNU GCC/G++. --------------------------------------------- */
#define max(a,b) ({ __typeof__(a) _a = (a);  __typeof__(b) _b = (b);  _a > _b ? _a : _b; })
#define min(a,b) ({ __typeof__(a) _a = (a);  __typeof__(b) _b = (b);  _a < _b ? _a : _b; })

#elif defined(_MSC_VER)
	/* Microsoft Visual Studio. --------------------------------- */
	
#endif

#define dim() (2)
/* Function prototypes */


#endif

