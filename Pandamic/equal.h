// Equal Distribution
// Created by Gunjan Basnet
// 04/04/2020

#include <map>
#include <vector>

using namespace std;

class equall
{
	private:
		int day;
		int infNum;
        int infRate;
		int infDay;
		int totalInfected;
        vector <int> vaccDistribution;
		vector<int> order,order1;
		vector<int> pop;
		vector<int> sus;
		vector<int> inf;
		vector<int> rec;
		vector<int> vac;
		map <int,int> storeInfc;
		vector<int> storeMaxinfc;
		int storeTotInfc;
		
	public:
		void dispenseVaccine(int, int, map<int, int> &);
		void SIRVsim(int, map<int, int> &, int&, int&, int&,  map<int, int> &);
		void simulation(int , vector<vector<int>>);
		void print();
};