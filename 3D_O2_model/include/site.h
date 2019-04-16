#ifndef SITE_H
#define SITE_H


typedef  double spin;
using namespace std;

class site{
	public:
		double * _components;
		int _id;
		vector <site *> _neighbour;
		site()
		{
			_components=NULL;
			_id=-1;
		}
};

#endif
