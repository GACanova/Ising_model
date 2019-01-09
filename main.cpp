#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <time.h>
#include <math.h>
#include <memory>
#include "mc.h"
#include "types.h"
#include "parameters.h"

using namespace std;

/********************* Function Declarations *******************************/
void Init(void);
void InitialState(spin *s);
void SingleFlip(spin *s,spin *M,spin *E,uint **neighbors,const double P[]);
void WolffAlgorithm(spin *s,uint **neighbors,double P);
void SquareLattice(uint **neighbors);
spin Energy(spin *s,uint **neighbors);
spin Magnetization(spin *s);
bool CheckRange(double T);
void ChangeT(double &T);

ofstream file;
string name;

/*********************Main Program***************************************/
int main(void)
{
	Init();

	spin *s=new spin[N];
	uint **neighbors=new uint*[N];
    	for(uint i=0;i<N;i++)neighbors[i]=new uint[4];

	SquareLattice(neighbors);
	InitialState(s);
	spin E=Energy(s,neighbors);
	spin M=Magnetization(s);

	for(double T=Ti;CheckRange(T);ChangeT(T))
	{
#if(algorithm=='S')
		double P[3]={0.,exp(-4./T),exp(-8./T)};
#elif(algorithm=='W')
		double P=1.-exp(-2./T);
#endif
		float e=0.,m=0.;

		for(uint t=0;t<transient;t++)
		{	
#if(algorithm=='S')
			SingleFlip(s,&M,&E,neighbors,P);
#elif(algorithm=='W')
			WolffAlgorithm(s,neighbors,P);
#endif
		}
	
		e=0.; m=0.;

		for(uint t=0;t<samples;t++)
		{
			for(uint j=0;j<tau;j++)
			{
#if(algorithm=='S')
				SingleFlip(s,&M,&E,neighbors,P);
#elif(algorithm=='W')
				WolffAlgorithm(s,neighbors,P);
#endif			
			}
#if(algorithm=='W')
			E=Energy(s,neighbors);
			M=Magnetization(s);
#endif
			e+=E; m+=1.*M; 
		}

		file << T << " " << 1.*e/(samples*N) << " " << 1.*m/(samples*N) << "\n";
		file.flush();
	}

	for(uint i=0;i<N;i++)delete[] neighbors[i];
	delete[] neighbors;
	delete[] s;
	file.close();
	return 0;
}
/********************* Initial state ******************************************************/
void InitialState(spin *s)
{
	for(uint i=0;i<N;i++)
	{
#if(initialState==1)
		s[i]=2*(int)(2.*RAND)-1;
#else
		s[i]=1;
#endif
	}
}
/************************* Single-Flip **************************************************/
#if(algorithm=='S')
void SingleFlip(spin *s,spin *M,spin *E,uint **neighbors,const double P[])
{
	for(uint n=0;n<N;n++)
	{
		const uint i=(int)(N*RAND);
		spin Sj=0;

		for(uint k=0;k<4;k++)
		{
			Sj+=s[neighbors[i][k]];
		}

		const spin dE=2*s[i]*Sj; 

		if(dE<=0)
		{
			s[i]*=-1;
			(*E)+=dE;
			(*M)+=2*s[i];
		}
		else if(RAND<P[dE>>2])
		{
			s[i]*=-1;
			(*E)+=dE;
			(*M)+=2*s[i];
		}
	}
}
#elif(algorithm=='W')
/*************Wolf-Algorithm**************************************************/
void WolffAlgorithm(spin *s,uint **neighbors,double P)
{
	uint *added=new uint[N];
	uint *list =new uint[N];	

	for(uint i=0;i<N;i++)
	{
		added[i]=0;
		list[i]=N;
	}

	uint i=(int)(N*RAND);
	s[i]*=-1;
	added[i]=1;
	uint l=0,n=0;

	while(i!=N)
	{
		for(uint k=0;k<4;k++)
		{
        		const uint j=neighbors[i][k];

        		if(!added[j])
        		{
				
				if(s[j]==-s[i])
				{
					if(RAND<P)
                			{
                        			s[j]*=-1;
                        			list[l]=j;
                        			l++;
                        			added[j]=1;
                			}
				}
				
        		}
		}

		i=list[n];
		n++;
	}

	delete[] added;
	delete[] list;
}
#endif
/*******************Square Lattice**********************************************/
void SquareLattice(uint **neighbors)
{
	for(uint i=0;i<N;i++)
    	{
        	uint x=i%L;
        	uint y=i/L;

        	neighbors[i][0]=(x+1)%L + y*L;      //East
        	neighbors[i][1]=(x-1+L)%L + y*L;    //West
        	neighbors[i][2]=x + ((y-1+L)%L)*L;  //North
        	neighbors[i][3]=x + ((y+1+L)%L)*L;  //South
       }
}
/****************************************************************************/
void Init(void)
{
	initRandom((long unsigned)SEED);

	ostringstream nameTmp;
	nameTmp << "ising_model_" << ((algorithm=='S')?("SingleFlip"):("Wolff")) <<"_algorithm_L_" << L << "_Ti_" << Ti << "_Tf_" << Tf << "_" <<RAND<< ".dat";
        name=nameTmp.str();
	file.open(name.c_str());
	file << "#Temperature, Energy, Magnetization"<<endl;
}
/****************************************************************************/
spin Energy(spin *s,uint **neighbors)
{
	spin E=0;

	for(uint i=0;i<N;i++)
	{
		for(uint j=0;j<4;j++)
		{
			E+=-s[i]*s[neighbors[i][j]];
		}	
	}
	return E/2;
}
/****************************************************************************/
spin Magnetization(spin *s)
{
	spin M=0;

	for(uint i=0;i<N;i++)
        {
		M+=s[i];
	}
	return abs(M);
}
/**************************************************************************/
bool CheckRange(double T)
{
#if(initialState==0)
	return T<=Tf;
#else
	return T>=Tf;
#endif
}
/**************************************************************************/
void ChangeT(double &T)
{
#if(initialState==0)
        T+=dT;
#else
        T-=dT;
#endif
}
