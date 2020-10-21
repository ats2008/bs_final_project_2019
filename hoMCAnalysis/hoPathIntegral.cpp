
#include <math.h>
#include<iostream>
#include<random>
#include<fstream>


#define RANDOM_ENGINE mt19937_64
#define RAND_IDX_MAX 512

#define SWEEP_COUNT 10000

using namespace std;
class lattice
{

	public:
		lattice(double mT=1,double wT=1,int nT=120,double dx=1.0,int randseed=0,int skipSweepCount=-1,int writeEventCount=128);
		~lattice() {} ;
		void initialize(string type="zero");
		void takeSweep(int nSteps);

		void printToASCII(int n=-1);
		void printLattice(bool printLatticeSites= false);
		void clearBuff();
		void clearFile();

	private:
		int N;
		double mTilda,wTilda;
		vector<double> xVec;
	        double action;
		double h,idrate;
		int initializationSeed;
		int skipSweepCount;
		int skipCounter;
		int writeOutCounter;
		int writeOutCount;
		long int stepCount;
		long int sweepCount;
		vector<double> actiondata;
		vector<double> xVecBuffer;
		vector<long int> stepCountData;
		vector<long int> sweepCountData;
		double alpha,beta;
	
		void findAction();
		void fillBuff();
		RANDOM_ENGINE generator;	
		uniform_int_distribution<int> intDistribution;	
		uniform_real_distribution<double> dblDistribution;	
		int randIdCounter;
		int randIdx[RAND_IDX_MAX];
		double randVals[RAND_IDX_MAX];
		void populateRandomNumbers();
		

		string oFileName;
};

 lattice::lattice(double mT,double wT ,int nT,double dx,int randseed,int skipSweepCount,int writeEventCount)
		: N(nT),mTilda(mT), wTilda(wT),
		  h(dx),idrate(0.8), xVec(N,0.0),
		  skipSweepCount(skipSweepCount>-1 ? skipSweepCount : 0  ),writeOutCount(writeEventCount>-1 ? writeEventCount: 1),
		  alpha(mTilda),beta(mTilda*wTilda*wTilda),
		  initializationSeed(randseed),
		  intDistribution(0,N-1),
		  dblDistribution(-0.5,0.5)
	//	  oFileName(ofname)
{
	initialize();
	oFileName="N"+to_string(N)+"_Nm_"+to_string(int(N*mTilda*1000)/1000)+".txt";

}



void lattice::initialize(string type)
{
	generator.seed(initializationSeed);

	if(type=="zero")
	{
		for(int i=0;i<N;i++)
			xVec[i]=0.0;
	}
	if(type=="hot")
	{
		for(int i=0;i<N;i++)
			xVec[i]=dblDistribution(generator);
	}
	populateRandomNumbers();
	findAction();
	writeOutCounter=0;
	skipCounter=0;
	stepCount=0;
	sweepCount=0;
	clearBuff();
	fillBuff();
}

void lattice::printLattice(bool printLatticeSites)
{
	cout<<"Printing Lattice with : "<<"\n";
	cout<<"              N  =  "<<N<<"\n";
	cout<<"         mTilda  =  "<<mTilda<<"\n";
	cout<<"         wTilda  =  "<<wTilda<<"\n";
	cout<<"         action  =  "<<action<<"\n";
	cout<<"              h  =  "<<h<<"\n";
	cout<<"         idrate  =  "<<idrate<<"\n";
	cout<<"      stepCount  =  "<<stepCount<<"\n";
	cout<<"     sweepCount  =  "<<sweepCount<<"\n";
	cout<<" skipSweepCount  =  "<<skipSweepCount<<"\n";
	cout<<"  writeOutCount  =  "<<writeOutCount<<"\n";

	if(printLatticeSites)
	{

	cout<<" Lattice is "<<"\n";
	for(int i=0;i<N;i++)
		cout<<" i = " <<i<<" -> "<<xVec[i]<<"\n";
	}
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

void lattice::takeSweep(int n)
{
	double percentStep=0.25;
	int oneP=n*percentStep/100;
	double percent =0;
	cout<<"stepCount = "<<stepCount<<" , sweepCount = " <<sweepCount<<" , "<<" skipCounter = "<<skipCounter<<"/"<<skipSweepCount<<" , " ;
	cout<<"writeOutCounter = "<<writeOutCounter<<"/"<<writeOutCount<<"  [ "<<percent<<"%  ]    ";
	percent+=percentStep;	

	for(int k=0;k<n;k++)
	{
	 double accrete=0;
	//	int mm=0;
	
	for(int i=0;i<N;i++)
	{
		auto idx=randIdx[randIdCounter];
		auto deltaX=randVals[randIdCounter]*h;
		randIdCounter++;
		if(randIdCounter==RAND_IDX_MAX)
			populateRandomNumbers();

		auto idxn = idx!=(N-1) ? idx+1 : 0;
		auto idxp = idx!=  0   ? idx-1 : N-1;
		//cout<<idx<<" "<<idxn<<" "<<idxp<<"\n";
		auto deltaS= deltaX*(alpha*(2*xVec[idx]+deltaX-xVec[idxn]-xVec[idxp])+beta*(xVec[idx]+deltaX/2));
		//cout<< deltaS<<" = "<<deltaX<<"*"<<"("<<alpha<<"*"<<"(2*"<<xVec[idx]<<"+"<<deltaX<<"-"<<xVec[idxn]<<"-"<<xVec[idxp]<<")+"<<beta<<"*"<<"("<<xVec[idx]<<"+"<<deltaX<<"/2))";
		
		if(deltaS<0)
		{
			xVec[idx]+=deltaX;
			accrete+= double(1.0/N);
	//		mm+=1;	cout<<"m = "<<mm<<" , i =  "<<i<<" accrete+="<<" double(1.0/"<<N<<") "<< accrete<<" , "<<double(1.0/N)<<"\n";
	
		}
		else if ((dblDistribution(generator)+0.5) < exp(-1*deltaS) )
		{       
			xVec[idx]+=deltaX;
			accrete+= double(1.0/N);
	//		mm+=1;	cout<<"m = "<<mm<<" , i =  "<<i<<" accrete+="<<" double(1.0/"<<N<<") "<< accrete<<" , "<<double(1.0/N)<<"\n";

		}
		else
			deltaS=0;
		//cout<<" dS = "<<deltaS<<" dX = "<<deltaX;
		action+=deltaS;
		//cout<<" action = "<<action<<"\n";
		stepCount+=1;

	}
	
       	        sweepCount+=1;
		skipCounter+=1;
	        
		if(skipCounter>=skipSweepCount)
		{
			fillBuff();
			skipCounter=0;
			writeOutCounter+=1;
			
		}
		if(writeOutCounter==writeOutCount)
		{
			//cout<<"\nwriting out "<<"\n";
			printToASCII(writeOutCount);
			clearBuff();
			writeOutCounter=0;
		}
		
		if(sweepCount%oneP==0)
		{
		cout<<"\r stepCount = "<<stepCount<<" , sweepCount = " <<sweepCount<<" , "<<" skipCounter = "<<skipCounter<<"/"<<skipSweepCount<<" , " ;
	        cout<<"writeOutCounter = "<<writeOutCounter<<"/"<<writeOutCount<<"  [ "<<percent<<"%  ]    ";//<<"\n";
		percent+=percentStep;	
		}
		
		//cout<<" h = "<<h<<" accrete = "<<accrete;//<<" \n";
		h*= accrete/idrate;
	}
	cout<<"\rstepCount = "<<stepCount<<" , sweepCount = " <<sweepCount<<" , "<<" skipCounter = "<<skipCounter<<"/"<<skipSweepCount<<" , " ;
	cout<<"writeOutCounter = "<<writeOutCounter<<"/"<<writeOutCount<<"  [ "<<percent<<"%  ]    ";

}


void lattice::printToASCII(int n)
{
	if(n<0)
		n=actiondata.size();

	fstream file(oFileName.c_str(),ios::app);
	file<<"\n#N : "<<N<<"\n";
	file<<"#xVecSize : "<<xVecBuffer.size()<<"\n";
	file<<"#nWrite : "<<n<<"\n";
	auto adataIt =actiondata.begin();
	auto xBuffIt = xVecBuffer.begin();
	auto sweepCountIt =sweepCountData.begin();
	for(long int i=0;i < n; i++)
	{
		//cout<<"\n"<<"printing " << stepCount-(n-i);
		file<<*sweepCountIt<<","<<*adataIt<<",";
		file<<*(xBuffIt);
		for(int j=1;j<N;j++)
		{
			xBuffIt++;
			file<<","<<*(xBuffIt);
		}
		adataIt++;
		sweepCountIt++;
		file<<"\n";

	}

	file.close();
}

void lattice::fillBuff()
{
	for(int i=0;i< xVec.size();i++)
		xVecBuffer.push_back(xVec[i]);
	actiondata.push_back(action);
	sweepCountData.push_back(sweepCount);
	stepCountData.push_back(stepCount);
}


void lattice::clearBuff()
{
	actiondata.clear();
	xVecBuffer.clear();
	sweepCountData.clear();
	stepCountData.clear();
}
void lattice::clearFile()
{
	cout<<"Purging : "<<oFileName<<"\n";
	fstream file(oFileName.c_str(),ios::out);
	file.close();
}


int main(int argc,char *argv[])
{

	double m=1,mN=120;
	int randseed=0;
	int skipSweepCount =100;
	int writeEventCount=128;
	double h=1;
	long int eventsRequired=8;
	long int sweepMaxCount(eventsRequired*skipSweepCount);

	if(argc>1)
		m=atof(argv[1]);
	if(argc>2)
		mN=atof(argv[2]);
	if(argc>3)
		skipSweepCount=int(atof(argv[3]));
	if(argc>4)
		eventsRequired=atof(argv[4]);

	sweepMaxCount=eventsRequired*skipSweepCount;
	int N=int(mN/m);
	double w=m;


	//lattice(double mT=1,double wT=1,int nT=120,double dx=1.0,int randseed=0,int skipSweepCount=-1,int writeEventCount=128);
	lattice alat(m,w,N,h,randseed,skipSweepCount,writeEventCount);
	alat.initialize("hot");
	alat.printLattice();

	cout<<"\nDoing simulation for "<<sweepMaxCount<<" sweeps for storing "<<eventsRequired<<" configurations \n";
	cout<<"\n";
	alat.clearFile();
	cout<<"\n";
	alat.takeSweep(sweepMaxCount);
	cout<<"\n";
	cout<<"\n";
	alat.printLattice();
	alat.printToASCII();
	cout<<"\n";
	//alat.takeStride(10);
	//alat.printLattice();
	return 0;
}

