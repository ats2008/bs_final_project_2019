#include "stattools.h"

using namespace std;

#define BSTRAP_COUNT 100
#define BSTRAP_IID_COUNT 20000

int main(int argc,char ** argv)
{
	string BASE_FOLDER  = "simulation_data/raw/";
	string DESTN_FOLDER = "simulation_data/resampled/";
	string ofname="data_analysis.txt";
	if(argc>1)
	{
			BASE_FOLDER=argv[1];
			DESTN_FOLDER=BASE_FOLDER+"analysis/";
			mkdir(DESTN_FOLDER.c_str(),0777);
	}
	else
	{
		cout<<"Please provide *source_path ,ofname,**destn_path,**destn_folder_name\n";
		cout<<"files are to be saved under destn_path/destn_folder_name \n";
		return 0;
	}
	if(argc>2)
		ofname=argv[2];
	if(argc>3)
		{
			ofname=argv[2];
			DESTN_FOLDER=argv[3];
			DESTN_FOLDER+="/";
			mkdir(DESTN_FOLDER.c_str(),0777);
		}
	
	string temp_str;
	
	//Reading names & tau of the simulation files
	temp_str=BASE_FOLDER+"Hfnames.txt";
	fstream file(temp_str.c_str(),ios::in);
	if(!file)
	{
		cout<<temp_str<<"\n ERROR OPENING FILE \n";
		return 0;
	}
	lattice_O_n alat;
	vector<string> fnames;
	string l="";
	while( getline (file,l,'\n') )
	{
		fnames.push_back(l);
	}
	file.close();
	
	
	// Resampling  the main files
	int i=0;
	double rslt[10],H[5];
	fstream ofile;
	ofname=DESTN_FOLDER+ofname;
	ofile.open(ofname.c_str(),ios::out);
	cout<<"\n WRITING RESULTS TO : "<<(ofname)<<"\n";
	ofile<<"#dim,L,N,T,J,H0,E,E_err,M,M_err,SP_heat,SP_heat_err,SUCEP,SUCEP_err,M0,M0_err,M0_BiCu,M0_BiCu_err,stats\n";
	cout<<"_______________________________________________\n";
	vector<string>::iterator name = fnames.begin();
	int BSTRAP_IID_COUNT_inside;
	for(;name!=fnames.end();name++)
	{
		i++;
		cout<<i<<" / "<<fnames.size()<<"\n";
		cout<<"READING file\t:"<< i <<" : "<<*name<<"\n";
		alat.readFromBinary(BASE_FOLDER+*name+".dat");
		alat.get_ext_field(H);
		ofile<<alat.get_diamension()<<","<<alat.get_lattice_len()<<","<<alat.get_site_count()<<","<<alat.get_temparature()<<",";
		ofile<<alat.get_J()<<","<<H[0]<<",";
		BSTRAP_IID_COUNT_inside=alat.energy.size();
		alat.bootstarap_basic(rslt,BSTRAP_COUNT,BSTRAP_IID_COUNT_inside);
		cout<<"@ dim = "<<alat.get_diamension();
		cout<<"  L = "<<alat.get_lattice_len()<<", T "<<alat.get_temparature();
		cout<<", J "<<alat.get_J()<<", H = "<<H[0]<<"\n";
		cout<<"Energy  = "<<rslt[0]<<" +/- "<<rslt[1]<<"\n";
		cout<<"Mag     = "<<rslt[2]<<" +/- "<<rslt[3]<<"\n";
		cout<<"SpHeat  = "<<rslt[4]<<" +/- "<<rslt[5]<<"\n";
		cout<<"Sucepti = "<<rslt[6]<<" +/- "<<rslt[7]<<"\n";
		ofile<<rslt[0]<<","<<rslt[1]<<","<<rslt[2]<<","<<rslt[3]<<",";
		ofile<<rslt[4]<<","<<rslt[5]<<","<<rslt[6]<<","<<rslt[7]<<",";
		
		alat.bootstarap_mgn(rslt,0,BSTRAP_COUNT,BSTRAP_IID_COUNT);
		cout<<"Mag x = "<<rslt[0]<<" +/- "<<rslt[1]<<"\n";
		ofile<<rslt[0]<<","<<rslt[1]<<" , ";
		
//		alat.bootstarap_mgn(rslt,1,BSTRAP_COUNT,BSTRAP_IID_COUNT);
//		cout<<"Mag y = "<<rslt[0]<<" +/- "<<rslt[1]<<"\n";
//		ofile<<rslt[0]<<","<<rslt[1]<<"\n";

		alat.bootstarap_bindersC(rslt,0,BSTRAP_COUNT,BSTRAP_IID_COUNT);
		cout<<"BindersC for MgX = "<<rslt[0]<<" +/- "<<rslt[1]<<"\n";
		ofile<<rslt[0]<<","<<rslt[1];
		ofile<<","<<alat.energy.size()<<"\n ";
		
		alat.clear_lattice();
	}
	ofile.close();
return 0;
}
