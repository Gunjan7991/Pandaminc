//Centrality Distribution

#ifndef DEGREE_H
#define DEGREE_H

#include <map>
#include <vector>

using namespace std;

class degree
{
private:
    int day;
    int infNum;
    int infRate;
    int infDay;
    int totalInfected;
    int storeTotInfc;
    vector<int> order, order1;
    vector<int> vaccDistribution;
    map<int, int> tempDeg;
    vector<int> pop;
    vector<int> sus;
    vector<int> inf;
    vector<int> rec;
    vector<int> vac;
    map<int, int> storeInfc;
    vector<int> storeMaxinfc;

public:
    void initialization(int);
    void disVaccine(int, map<int, int> &, map<int, int> &);
    void sorting(int, map<int, int> &);
    void SIRVsim(int, map<int, int> &, int &, int &, int &);
    void simulation(int, vector<vector<int>>);
    void print();
};

#endif