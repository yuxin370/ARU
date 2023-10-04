#ifndef DECISIONT_H
#define DECISIONT_H
/*
 *  Definition of decision table
 */
#include <fstream>
#include <vector>
#include <set>
#include <unordered_set>
#include <map>
#include <string>
#include <time.h>
using namespace std;


class DecisionT
{
public:
    DecisionT();
    DecisionT(ifstream& in,vector<int>& Cs, int Ds);
    DecisionT(
            const vector<vector<double> >& records,
            size_t len,
            vector<int>& Cs,
            int Ds);
    DecisionT(
            vector<int>& Cs,
            int Ds);
    DecisionT(DecisionT& d);

    DecisionT& getcopy() const;
    set<int>& getCs();
    vector<vector<double> >& getrecords();
    int getDs();
    void addrecords(vector<vector<double> >::const_iterator it,int);
    void delrecords(int size); 

private:
    vector<vector<double> > records;
    set<int> Cs;
    int Ds;

};

string getkey(const vector<double> &x, int i);
double Dependent(DecisionT& dt, set<int>& B,int holdpos,...); 
void POSplus(
        vector<vector<double> > records,
        vector<vector<double> >::const_iterator iold,
        size_t presize,
        vector<vector<double> >::const_iterator inew,
        size_t gsize,
        set<int>& B,
        int ds,
        vector<double>& prepos
        );
double pos_x(
        const vector<vector<double> >& records,
        size_t x,
        set<int> & B,
        int ds);
void POS(
        DecisionT& dt,
        set<int>& B,
        vector<double> & pos
        );
#endif
