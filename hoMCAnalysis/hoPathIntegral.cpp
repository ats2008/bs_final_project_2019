
#include <math.h>
#include<iostream>
#include<random>

#define RANDOM_ENGINE mt19937_64
#define RAND_IDX_MAX 512

using namespace std;
class lattice
{

	public:
    		lattice(int N=10,double dt=1, double m=2,double w=3,int randseed=0);
		~lattice() {} ;
		bool initialize(string type="zero");
		void takeSteps(int nSteps);
		void takeStride(int nStrides=1);

		void printLattice();
	private:
		int N;
		double dTau,mTilda,wTilda,xTildaFactor;
		vector<double> xVec;
		vector<double> actiondata;
	        double action;
		int skipDataCount;
		int initializationSeed;
		double alpha,beta;

		void findAction();

		RANDOM_ENGINE generator;	
		uniform_int_distribution<int> intDistribution;	
		uniform_real_distribution<double> dblDistribution;	
		int randIdCounter;
		int randIdx[RAND_IDX_MAX];
		double randVals[RAND_IDX_MAX];
		void populateRandomNumbers();
		
};


 lattice::lattice(int n,double dt, double m,double w ,int randseed)
		: N(n),dTau(dt),mTilda(m/dt), wTilda(w/dt),xTildaFactor(1.0/dt),
		  xVec(N,0.0),
		  alpha(mTilda),beta(mTilda*wTilda*wTilda),
		  initializationSeed(randseed),
		  intDistribution(0,N-1),
		  dblDistribution(-0.5,0.5)
{
	initialize();
}

bool lattice::initialize(string type)
{


	if(type=="zero")
	{
		for(int i=0;i<N;i++)
			xVec[i]=0.0;
	}
	if(type=="hot")
	{
		generator.seed(initializationSeed);
		for(int i=0;i<N;i++)
			xVec[i]=dblDistribution(generator);
	}
	populateRandomNumbers();
	findAction();
	actiondata.push_back(action);
}

void lattice::printLattice()
{
	cout<<"Printing Lattice with : "<<"\n";
	cout<<"              N  =  "<<N<<"\n";
	cout<<"             dt  =  "<<dTau<<"\n";
	cout<<"         mTilda  =  "<<mTilda<<"\n";
	cout<<"         wTilda  =  "<<wTilda<<"\n";
	cout<<"   xTildaFactor  =  "<<xTildaFactor<<"\n";
	cout<<"         action  =  "<<action<<"\n";

	cout<<" Lattice is "<<"\n";
	for(int i=0;i<N;i++)
		cout<<" i = " <<i<<" -> "<<xVec[i]<<"\n";
}

void lattice::findAction()
{
	double a(0.0),b(0.0);
	for(int i=0;i<N-1;i++)
	{
		a+=(xVec[i+1]-xVec[i])*(xVec[i+1]-xVec[i]);
		b+=xVec[i]*xVec[i];
	}
	
	a+=(xVec[N-1]-xVec[0])*(xVec[N-1]-xVec[0]);
	b+=xVec[N-1]*xVec[N-1];

	action=0.5*mTilda*a+0.5*wTilda*b;
}

void lattice::populateRandomNumbers()
{
	
	for(int i=0;i<RAND_IDX_MAX;i++)
	{
		randIdx[i]=intDistribution(generator);
		randVals[i]=dblDistribution(generator);
	}
	randIdCounter=0;
}

void lattice::takeSteps(int n)
{
	for(int i=0;i<n;i++)
	{
		
		auto idx=randIdx[randIdCounter];
		auto deltaX=randVals[randIdCounter]*xTildaFactor;
		randIdCounter++;
		if(randIdCounter==RAND_IDX_MAX)
			populateRandomNumbers();

		auto idxn = idx!=(N-1) ? idx+1 : 0;
		auto idxp = idx!=  0   ? idx-1 : N-1;
		auto deltaS= deltaX*(alpha*(2*xVec[idx]+deltaX-xVec[idxn]-xVec[idxp])+beta*(xVec[i]+deltaX/2));
		
		if(deltaS<0)
		       	xVec[idx]+=deltaX;
		else if (dblDistribution(generator) < exp(-1*deltaS) )
		       	xVec[idx]+=deltaX;
		else
			deltaS=0;
		cout<<" dS = "<<deltaS<<" dX = "<<deltaX<<"\n";
		action+=deltaS;
		cout<<" action = "<<action<<"\n";
		actiondata.push_back(action);
	
	}
}

void lattice::takeStride(int n)
{
	takeSteps(N*n);
}

int main()
{
	
	lattice alat;
	//alat.initialize("hot");
	alat.printLattice();
	alat.takeStride();
	alat.printLattice();

	return 0;

}
