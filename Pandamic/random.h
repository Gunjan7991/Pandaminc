// Random Distribution
// created by Gunjan Basnet
// 04/04/2020

#include <map>
#include <vector>

using namespace std; 

class randoom{
	private:
		int infRate;
		int infDay;
		int totalInfected;
		int day;
		int infNum;
		vector<int> order,order1;
		vector<int> pop;
		vector<int> sus;
		vector<int> inf;
		vector<int> rec;
		vector<int> vac;
        vector <int> vaccDistribution;
		map <int,int> storeInfc;
		vector<int> storeMaxinfc;
		int storeTotInfc;

	public:
		void dispenceVaccine(int, int, map<int, int> &);
		void SIRVsim(int, map<int, int> &, int&, int&, int&,  map<int, int> &);
		void simulation(int , vector<vector<int>>);
		void print();
};