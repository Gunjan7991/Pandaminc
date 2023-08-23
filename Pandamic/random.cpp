// Random Distribution
// Created by Gunjan Basnet
// 04/04/2020

#include <map>
#include <vector>
#include <iostream>
#include <algorithm>
#include <iomanip>
#include <math.h>

#include "random.h"

void randoom::dispenceVaccine(int numline, int numVac, map<int, int> &population)
{
    for (int i = 0; i < numline - 1; i++) //used to initialize the vacc distribution with 0;
    {
        vaccDistribution.push_back(0);
    }
    int vaccine = numVac;
    map<int, int>::iterator itr = population.begin();
    for (int i = 0; i < vaccDistribution.size(); i++) //loop runs once for all population.
    {
        if (vaccine >= 0) //untill vaccine runs out.
        {
            if (vaccine > itr->second) //is vaccine is more than the population of the particular region region
            {
                vaccine = vaccine - itr->second;
                vaccDistribution[i] = itr->second;
            }
            else //if vaccine is less than the population of the region.
            {
                vaccDistribution[i] = vaccine;
                vaccine = -1; //run population out after distributing vaccine.
            }
        }
        ++itr;
    }
}

void randoom::SIRVsim(int numline, map<int, int> &mapPopn, int &infArea, int &infPeriod, int &conRate, map<int, int> &degCount)
{
    day = 0;
    totalInfected = 0;
    infRate = conRate;
    infDay = infPeriod;

    map<int, int>::iterator ptr;
    for (ptr = degCount.begin(); ptr != degCount.end(); ++ptr)
    {
        order1.push_back(ptr->second);
    }

    for (int i = 0; i < (numline - 1); i++)
    {
        sus.push_back(0);
        inf.push_back(0);
        rec.push_back(0);
        vac.push_back(0);
    }

    for (int i = 0; i < (numline - 1); i++)
    {
        vac[i] = vaccDistribution[i];
    }

    map<int, int>::iterator itr;
    for (itr = mapPopn.begin(); itr != mapPopn.end(); ++itr)
    {
        pop.push_back(itr->second);
    }

    for (int i = 0; i < (numline - 1); i++)
    {
        sus[i] = pop[i] - vac[i];
    }
    if (sus[(infArea - 1)] > 0)
    {
        inf[(infArea - 1)] = 1;
        sus[(infArea - 1)] = sus[(infArea - 1)] - inf[(infArea - 1)];
    }

    for (int i = 0; i < numline - 1; i++)
    {
        totalInfected = totalInfected + inf[i];
    }

    infNum = 1;
    for (int j = 0; j < infDay; j++)
    {
        infNum = infNum + infNum * infRate;
    }

    //Print SIRV Model
    cout << endl
         << "RANDOM DISTRIBUTION" << endl;
    cout << "Day " << day << endl;
    for (int i = 0; i < (numline - 1); i++)
    {
        cout << fixed << left << setw(4) << i + 1 << setw(5) << "POP: " << setw(7) << pop[i] << setw(3) << "S: " << setw(7) << sus[i] << setw(3) << "I: " << setw(7) << inf[i] << setw(3) << "R: " << setw(7) << rec[i] << setw(3) << " V: " << setw(7) << vac[i] << endl;
    }
    cout << endl;
    storeTotInfc = 0;
}

void randoom::simulation(int numline, vector<vector<int>> adjList)
{
    vector<int> prevRec, prevRec2, prevInfc, l1, l2, l3, l4;
    vector<int> power, timer1, timer2, timer3;
    int rat = infRate + 1;
    for (int i = 0; i < numline - 1; i++)
    {
        prevRec.push_back(0);
        prevRec2.push_back(0);
        power.push_back(0);
        timer1.push_back(0);
        timer2.push_back(0);
        timer3.push_back(0);
        l1.push_back(1);
        l2.push_back(1);
        l3.push_back(1);
        l4.push_back(1);
    }
    int count = 0;

    while (totalInfected != 0)
    {
        vector<double> half, pophalf;
        totalInfected = 0;
        day++;
        count++;
        vector<int> adjNumbers;
        for (int i = 0; i < numline - 1; i++)
        {
            pophalf.push_back(pop[i] / 2);
            half.push_back(inf[i]);
            if (half[i] > pophalf[i])
            {
                for (int j = 0; j < order1[i]; j++)
                {
                    adjNumbers.push_back(adjList[i][j]);
                }
            }
        }

        for (int i = 0; i < numline - 1; i++)
        {
            prevRec[i] = rec[i];
            if ((inf[i] * (infRate + 1)) == infNum)
            {
                rec[i] = 1;
            }
            else if (rec[i] > 0)
            {
                rec[i] = (rec[i] * rat);
                if (power[i] > (infDay - 1) && rec[i] <= (pop[i] - vac[i]))
                {
                    rec[i] = rec[i] - (pow(rat, timer1[i]) * infRate);
                    if (timer1[i] >= (infDay + 1))
                    {
                        rec[i] = rec[i] + (pow(rat, timer2[i]) * (pow(infRate, (timer2[i] + 2))));
                        if (timer2[i] >= (infDay))
                        {
                            rec[i] = rec[i] - (pow(rat, (timer3[i] + 2)) * pow(infRate, (timer3[i] + 2)));
                        }
                        else if (rec[i] > (pop[i] - vac[i]))
                        {
                            rec[i] = (pop[i] - vac[i]);
                        }
                        timer2[i]++;
                    }
                    else if (rec[i] > (pop[i] - vac[i]))
                    {
                        rec[i] = (pop[i] - vac[i]);
                    }
                    timer1[i]++;
                }
                else if (rec[i] > (pop[i] - vac[i]))
                {
                    rec[i] = (pop[i] - vac[i]);
                }
                power[i]++;
            }
            if (((inf[i] * (infRate + 1)) - rec[i] + prevRec[i]) < (pop[i] - vac[i] - rec[i]))
            {
                inf[i] = (inf[i] * rat) - rec[i] + prevRec[i];
            }
            else
            {
                inf[i] = pop[i] - vac[i] - rec[i];
            }
            prevRec[i] = rec[i];
        }

        for (int i = 0; i < numline - 1; i++)
        {
            if (adjNumbers.size() > 0 && half[i] > pophalf[i])
            {
                for (int j = 0; j < adjNumbers.size(); j++)
                {
                    if (sus[adjNumbers[j]] > 0 && inf[adjNumbers[j]] == 0)
                    {
                        inf[(adjNumbers[j])] = 1;
                    }
                }
            }
        }

        for (int i = 0; i < numline - 1; i++)
        {
            if (sus[i] >= 0)
            {
                sus[i] = pop[i] - inf[i] - rec[i] - vac[i];
                if (sus[i] != pop[i] - inf[i] - rec[i] - vac[i])
                {
                }
            }
            else
            {
                sus[i] = 0;
            }
        }

        cout << "Day " << day << endl;
        for (int i = 0; i < (numline - 1); i++)
        {
            cout << fixed << left << setw(4) << i + 1 << setw(5) << "POP: " << setw(7) << pop[i] << setw(3) << "S: " << setw(7) << sus[i] << setw(3) << "I: " << setw(7) << inf[i] << setw(3) << "R: " << setw(7) << rec[i] << setw(3) << " V: " << setw(7) << vac[i] << endl;
        }
        cout << endl;

        for (int i = 0; i < numline - 1; i++)
        {
            totalInfected = totalInfected + inf[i];
        }
        storeMaxinfc.push_back(totalInfected);
        storeInfc.insert(pair<int, int>(day, totalInfected));
    }
    sort(storeMaxinfc.begin(), storeMaxinfc.end(), greater<int>());
    for (int i = 0; i < numline - 1; i++)
    {
        storeTotInfc = storeTotInfc + rec[i];
    }
}

void randoom::print()
{
    int maxDay, maxInfc;
    if (storeMaxinfc.size() > 0)
    {
        maxInfc = storeMaxinfc[0];
        map<int, int>::iterator itr;
        for (itr = storeInfc.begin(); itr != storeInfc.end(); ++itr)
        {
            if (maxInfc == itr->second)
            {
                maxDay = itr->first;
            }
        }
        cout << "Using the random distribution method, the peak number of infected was " << maxInfc << " on day " << maxDay << ". The outbreak ended on day " << day << " and the total number of infected individuals was " << storeTotInfc << "." << endl;
    }
    else
    {
        cout << "Using the random distribution method, the peak number of infected was " << 0 << " on day " << 0 << ". The outbreak ended on day " << 0 << " and the total number of infected individuals was " << 0 << "." << endl;
    }
}