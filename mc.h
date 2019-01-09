#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <limits.h> 
#include <string.h>
#include <time.h>  

#define NORM   (2.3283064365e-10)
#define RANDOM  ( (ira[ip++]=ira[ip1++]+ira[ip2++]) ^ira[ip3++] )        
#define RAND (NORM * RANDOM)

/********************************************************************
***                      Variable Declarations                    ***
********************************************************************/
unsigned zseed,ira[256];
unsigned char ip,ip1,ip2,ip3;
/********************************************************************
*            Random Number Generator by Parisi & Rapuano            *
*            RAND -> double in [0,1)                               *
********************************************************************/
unsigned rand4init(long unsigned seed)
{
     unsigned long long y;
     static long unsigned zseed=seed;
   
     y = (zseed*16807LL);
     zseed = (y&0x7fffffff) + (y>>31);

     if(zseed&0x80000000)
     	zseed = (zseed&0x7fffffff) + 1;
     return zseed;
}

void initRandom(unsigned long seed)
{
	if((seed%2)==0)seed++;

     ip=128;
     ip1=ip-24;
     ip2=ip-55;
     ip3=ip-61;
   
     for(long i=ip3; i<ip; i++)
	     ira[i] = rand4init(seed);
}
