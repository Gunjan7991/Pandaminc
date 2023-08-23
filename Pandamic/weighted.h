// Weighted Distribution
// WEIGHTED

#ifndef WEIGHTED
#define WEIGHTED

#include <map>
#include <vector>

using namespace std;

class weighted
{
private:
    int infRate;
    int infDay;
    int totalInfected;
    int day;
    int infNum;
    vector<int> vaccDistribution;
    vector<int> order, order1;
    vector<int> pop;
    vector<int> sus;
    vector<int> inf;
    vector<int> rec;
    vector<int> vac;
    map<int, int> storeInfc;
    vector<int> storeMaxinfc;
    int storeTotInfc;

public:
    void dispenVaccine(int, int, map<int, int> &);
    void SIRVsim(int, map<int, int> &, int &, int &, int &, map<int, int> &);
    void simulation(int, vector<vector<int>>);
    void print();
};

#endif