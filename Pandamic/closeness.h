// Closeness Centrality
#ifndef CLOSENESS_H
#define CLOSENESS_H

#include <map>
#include <vector>
#include <iostream>

using namespace std;

class closeness
{
private:
    int day;
    int infNum;
    int storeTotInfc;
    int infRate;
    int infDay;
    int totalInfected;
    vector<double> closeness;
    vector<int> vaccDistribution;
    map<int, int> tempDeg;
    vector<int> order, order1;
    vector<int> pop;
    vector<int> sus;
    vector<int> inf;
    vector<int> rec;
    vector<int> vac;
    map<int, int> storeInfc;
    vector<int> storeMaxinfc;

public:
    void initialization(int);
    void disVaccine(int, map<int, int> &, map<int, int> &, map<int, double> &);
    void sorting(int, map<int, double> &);
    void SIRVsim(int, map<int, int> &, int &, int &, int &, map<int, int> &);
    void simulation(int, vector<vector<int>>, map<int, double> &);
    void print();
};

#endif