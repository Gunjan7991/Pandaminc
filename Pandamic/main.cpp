#include <iostream>
#include <string>
#include <fstream>
#include <math.h>
#include <vector>
#include <map>
#include <iomanip>
#include <iterator>
#include <queue>

#include "closeness.h"
#include "degree.h"
#include "random.h"
#include "equal.h"
#include "weighted.h"

using namespace std;

bool readFile(int &, int &, int &, int &, int &, map<int, int> &, vector<vector<char>> &); //used to read the provided files and save them to use it later.
void printRegionalPopulation(map<int, int>);                                               // used to print Regional Population.
void printAdjacentMatrix(int, vector<vector<char>> &, vector<vector<int>> &, map<int, int> &);
void bfs(int &, map<int, int> &, vector<vector<int>> &, map<int, double> &);

int main()
{
    closeness close; //creating object of the class.
    degree deg;   // creating object of the class.
    randoom rant;     //creating object of the class.
    equall eq;       //creating object of the class.
    weighted weig;   //creating object of the class.

    int numLine;                                  //used to store the number of region.
    int infArea, infPeriod, conRate, numVaccines; // infected area, infectious period, contact rate, no. of vaccines;
    map<int, int> mapPopn;                        // Map to store the values of Population
    map<int, int> degSave;
    map<int, double> closeSave;                                                               // map container used to store the value of closeness centrality along with its corrosponding region.
    vector<vector<int>> adjList;                                                              //adjcent list : it is used to store the provided adjecent matrix but in a number system that is easy to read by user.
    vector<vector<char>> adjMat;                                                              //used to save the provided adjency matrix
    bool open = readFile(numLine, infArea, infPeriod, conRate, numVaccines, mapPopn, adjMat); //check if the given file was opened and  stores the boolen value returned by the function

    if (open == true)
    {
        printRegionalPopulation(mapPopn);                       //used to print the region and their corrosponting population.
        printAdjacentMatrix(numLine, adjMat, adjList, degSave); //function used to store and print the adjecent matrix in user freindly -readable way. degree of the region are also stored in this step.
        bfs(numLine, degSave, adjList, closeSave);              //functionused to breath first search method, and calculate the closeness of each reagion.
        //printRegionalPopulation(degSave); //need to comment out later.

        close.initialization(numLine);
        close.sorting(numLine, closeSave);
        close.disVaccine(numVaccines, mapPopn, degSave, closeSave);
        close.SIRVsim(numLine, mapPopn, infArea, infPeriod, conRate, degSave);
        close.simulation(numLine, adjList, closeSave);

        deg.initialization(numLine);
        deg.sorting(numLine, degSave);
        deg.disVaccine(numVaccines, mapPopn, degSave);
        deg.SIRVsim(numLine, mapPopn, infArea, infPeriod, conRate);
        deg.simulation(numLine, adjList);

        rant.dispenceVaccine(numLine, numVaccines, mapPopn);
        rant.SIRVsim(numLine, mapPopn, infArea, infPeriod, conRate, degSave);
        rant.simulation(numLine, adjList);

        eq.dispenseVaccine(numLine, numVaccines, mapPopn);
        eq.SIRVsim(numLine, mapPopn, infArea, infPeriod, conRate, degSave);
        eq.simulation(numLine, adjList);

        weig.dispenVaccine(numLine, numVaccines, mapPopn);
        weig.SIRVsim(numLine, mapPopn, infArea, infPeriod, conRate, degSave);
        weig.simulation(numLine, adjList);

        close.print();
        deg.print();
        rant.print();
        eq.print();
        weig.print();
    }
    return 0;
}

bool readFile(int &numLines, int &infArea, int &infPeriod, int &conRate, int &numVaccines, map<int, int> &mapPopn, vector<vector<char>> &onlyMat)
{
    int count = 1;
    vector<vector<char>> adjMat;
    string config_filename;   // to store configuration filename
    string popnFile, regFile; // to store population filename, region filename,
    string popnText, regText; // to store the contents of Population and Region Files
    int pText;                // stores the converted value of population to int
    string currentPopn;
    string temp; // to store a particular population from the Population file
    vector<char> adjlist;
    string infA, infP, conR, numV; //variables that store converted string to int

    cout << "Please enter the name of the configuration file:";
    getline(cin, config_filename);
    ifstream openConfigfile(config_filename);

    if (!openConfigfile.is_open()) //use to check if file is opened or not, and appropiate message is displayed if not opened.
    {
        cout << "Failed to open the file. Please try again!" << endl;
        return false;
    }

    else
    {
        // Get the Population Text-file
        openConfigfile.ignore(256, ':');   // Ignore until a colon appears
        getline(openConfigfile, popnFile); //first line opened. // name of population file stored in popnFile

        // popnFile stores the name of the Population Filename
        ifstream openPFile(popnFile);    // to find the number of lines in the Population file.
        ifstream openPopnFile(popnFile); //to get the region and its population.

        // Counting the number of lines in Population file
        numLines = 0;
        string linePopn;
        do
        {
            getline(openPFile, linePopn);
            numLines++;
        } while (!openPFile.eof());

        do
        {
            openPopnFile.ignore(256, ':');
            getline(openPopnFile, popnText);
            pText = stoi(popnText);                       //Converting the value read from file to int (for later manipulation)
            mapPopn.insert(pair<int, int>(count, pText)); //store region and its population.
            count++;
        } while (count != numLines);

        // Get the Region CSV-file
        openConfigfile.ignore(256, ':');  // Ignore until a colon appears
        getline(openConfigfile, regFile); //second line opened. name of the region file saved.
        // regFile stores the name of the Region Filename

        //Printing and Storing the Adjacency List
        ifstream openRegFile(regFile); // region file opened

        int maxChar = 0;
        char currentChar;
        getline(openRegFile, regText); //got inside the region file, skipied the first line.
        vector<int> linesize;          // used to save the region data linearly.
        while (getline(openRegFile, regText))
        {
            for (int i = 0; i < regText.size(); i++)
            {
                currentChar = regText.at(i);
                if (currentChar != ',' && currentChar != '\n') //used to skip ',' and newline
                {
                    adjlist.push_back(currentChar); //only stores elements without ',' and newline from the file to the vector
                }
            }
            maxChar = ceil(regText.size() / 2 + 1); // later used to count and skip the first char.
            linesize.push_back(maxChar);            //store the maxChar of each line. later used to count and skip the first char.
        }

        //Storing the elements in Vector
        int x = 0;
        for (int i = 0; i < (numLines - 1); i++) // used to store the linerar data in 2 d form
        {
            if (linesize[i] == numLines) //used if and else to store all the element. usefull in the case when the region number consist of 2 char.
            {
                vector<char> temp;
                for (int j = 0; j < (numLines); j++)
                {
                    temp.push_back(adjlist.at(x));
                    x++;
                }
                adjMat.push_back(temp);
            }
            else
            {
                vector<char> temp;
                for (int j = 0; j < (numLines + 1); j++)
                {
                    temp.push_back(adjlist.at(x));
                    x++;
                }
                adjMat.push_back(temp);
            }
        }

        //Storing the linear elements in 2D-matrix form
        for (int i = 0; i < numLines - 1; i++) // restored the previous matrix by removing all the first character.
        {
            if (linesize[i] == numLines)
            {
                vector<char> temp;
                for (int j = 1; j < (numLines); j++)
                {
                    temp.push_back(adjMat[i][j]);
                }
                onlyMat.push_back(temp);
            }
            else
            {
                vector<char> temp;
                for (int j = 2; j < (numLines + 1); j++)
                {
                    temp.push_back(adjMat[i][j]);
                }
                onlyMat.push_back(temp);
            }
        }

        // Get the Infected Area
        openConfigfile.ignore(256, ':'); // Ignore until a colon appears
        getline(openConfigfile, infA);
        infArea = stoi(infA);

        // Get the Infectious Period
        openConfigfile.ignore(256, ':'); // Ignore until a colon appears
        getline(openConfigfile, infP);
        infPeriod = stoi(infP);

        // Get the Contact Rate
        openConfigfile.ignore(256, ':'); // Ignore until a colon appears
        getline(openConfigfile, conR);
        conRate = stoi(conR);

        // Get the Number of Vaccines
        openConfigfile.ignore(256, ':'); // Ignore until a colon appears
        getline(openConfigfile, numV);
        numVaccines = stoi(numV);

        openConfigfile.close(); //Closing Configuration File
        openPFile.close();
        openPopnFile.close();
        openRegFile.close();
        return true;
    }
}

void printRegionalPopulation(map<int, int> mapPopn)
{
    //Printing Regional Population using Iterator for Map:
    cout << "Regional Population" << endl;
    map<int, int>::iterator itr;
    for (itr = mapPopn.begin(); itr != mapPopn.end(); ++itr)
    {
        cout << itr->first << " " << itr->second << '\n';
    }
    cout << endl;
}

void printAdjacentMatrix(int numLines, vector<vector<char>> &onlyMat, vector<vector<int>> &adjList, map<int, int> &degSave)
{
    int actualLines = numLines - 1;
    int remainingLines = 0;
    //Printing the Adjacency List
    cout << "Adjacency List" << endl;
    for (int i = 0; i < (numLines - 1); i++)
    {
        vector<int> temp;
        int degCount = 0;
        cout << i + 1 << ": ";

        for (int j = 0; j < (numLines - 1); j++)
        {
            if (onlyMat[i][j] == '1')
            {
                degCount++;
                temp.push_back(j);
                cout << j + 1 << " ";
            }
        }
        adjList.push_back(temp);                         // used to save adj list.
        degSave.insert(pair<int, int>(i + 1, degCount)); //store the region and its degree centrality.
        cout << endl;
    }
}

void bfs(int &numLine, map<int, int> &degSave, vector<vector<int>> &adjList, map<int, double> &hait)
{
    vector<int> checklist;
    for (int i = 0; i < numLine - 1; i++)
    {
        checklist.push_back(i); // saved a list of region.
    }
    vector<int> pathLength; //used to store the total path length of each region.

    for (int i = 0; i < checklist.size(); i++)
    {                         // loop run untill all the region finds connection with all other regions
        queue<int> lastPlace; // queue of last element
        double buffsize = 1, size2 = 1, count1 = 0, count2 = 0;
        double incSize1 = 0, incSize2 = 0;
        vector<int> visited; //used to store so can be used to check if the region is already visited or not.
        int place = 0;
        queue<int> vis;
        visited.push_back(i);
        vis.push(i);
        while (!vis.empty() && visited.size() < numLine)
        { //loop runs untill the all the region is visited.
            int row = vis.front();
            int pow = vis.back();
            vis.pop();
            map<int, int>::iterator itr = degSave.begin();
            int column;
            while (itr != degSave.end())
            {
                if ((row + 1) == itr->first)
                {
                    column = itr->second;
                }
                ++itr;
            }
            for (int j = 0; j < column; j++) //used to add all the element of the row.
            {
                bool cond = false;
                for (int k = 0; k < visited.size(); k++) //used to check if the element is already visted or not
                {
                    if (adjList[row][j] == visited[k])
                    {
                        cond = true;
                    }
                }
                if (cond == false)
                {
                    visited.push_back(adjList[row][j]);
                    vis.push(adjList[row][j]);
                    count1 = count1 + count2 + 1;
                }
            }

            //method used to increse the second count, in order to calculate closeness.
            if (incSize2 == 0)
            {
                size2 = visited.size();
                incSize2 = size2 - buffsize;
                buffsize = size2;
                count2++;
            }

            if (incSize2 > 0)
            {
                incSize2--;
            }
            int check = vis.front();
        }
        double why = checklist.size();
        double when = count1 / why;
        hait.insert(pair<int, double>(i + 1, when)); //region and its corrosponding closeness saved.
                                                     //while loop ends here.
    }
}