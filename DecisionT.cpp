#include <iostream>
#include <fstream>
#include <vector>
#include <set>
#include <map>
#include <cstdarg>
#include <sstream>
#include <string>
#include <algorithm>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <float.h>
#include "DecisionT.h"
using namespace std;
const double abc = 1e-8;

/**
 * @brief get records
 * 
 * @return vector<vector<double> >& 
 */
vector<vector<double> >& DecisionT::getrecords()
{
    return this->records;
}

/**
 * @brief transfer a vector to set 
 *  
 * @param a dest set
 * @param b source vector
 */
static void trans_set(set<int>& a, vector<int> & b)
{
    for(vector<int>::iterator it = b.begin();
            it != b.end();
            ++it)
    {
        a.insert(*it);
    }
}

/**
 * @brief Construct a new Decision T:: Decision T object
 * 
 * @param in : input file
 * @param Cs : continious attributes location 
 * @param Ds : decision attributes location
 */
DecisionT::DecisionT(ifstream& in,vector<int>& Cs, int Ds)
{

    this->Ds = Ds;
    trans_set(this->Cs, Cs);

    string line;
    while(getline(in,line))
    {
        istringstream stream(line);

        vector<double> *record = new vector<double>();
        string field;
        while(getline(stream,field,','))
        {
            record->push_back(atof(field.c_str()));
        }
        this->records.push_back(*record);
        int val = (*record)[0];
    }
}

/**
 * @brief Construct a new Decision T:: Decision T object
 * 
 * @param records : instance records
 * @param len : size of instance records
 * @param Cs : continious attributes location 
 * @param Ds : decision attributes location
 */
DecisionT::DecisionT(
        const vector<vector<double> >& records,
        size_t len,
        vector<int>& Cs, 
        int Ds)
{
    this->Ds = Ds;
    trans_set(this->Cs, Cs);

    for(vector<vector<double> >::const_iterator record = records.begin();
            record != records.begin()+len;
            ++record)
    {
        this->records.push_back(*record);
        int val = (*record)[0];
    }
}

/**
 * @brief Construct a new Decision T:: Decision T object
 * 
 * @param Cs : continious attributes location 
 * @param Ds : decision attributes location
 */
DecisionT::DecisionT(
        vector<int>& Cs,
        int Ds)
{
    this->Ds = Ds;
    trans_set(this->Cs, Cs);
}

/**
 * @brief Construct a new Decision T:: Decision T object
 * 
 * @param d decision table
 */
DecisionT::DecisionT(DecisionT& d){
    this->records = d.records;
    this->Cs = d.getCs();
    this->Ds = d.getDs();
}

/**
 * @brief Construct a new Decision T:: Decision T object
 * 
 */
DecisionT::DecisionT()
{
}

/**
 * @brief get continious attribute location set
 * 
 * @return set<int>& 
 */
set<int>& DecisionT::getCs()
{
    return this->Cs;
}

/**
 * @brief get decision attribute location
 * 
 * @return int 
 */
int DecisionT::getDs()
{
    return this->Ds;
}

/**
 * @brief add records
 * 
 * @param iBegin :the iterator of the records to be added 
 * @param size : the size of the records to be added.
 */
void DecisionT::addrecords(vector<vector<double> >::const_iterator iBegin,int size)
{
    for(vector<vector<double> >::const_iterator record = iBegin;
            record !=iBegin+size;
            ++record)
    {
        this->records.push_back(*record);
    }

}

/**
 * @brief delete records at the end of the records
 * 
 * @param size 
 */
void DecisionT::delrecords(int size){
    for(int i = 0 ;i < size ; i ++)
    {
        if(!records.empty()){
            this->records.pop_back();
        }
    }
}


/**
 * @brief calculate the distance between to instance on the attribute
 * 
 * @param Xi : instance Xi
 * @param Xj : instance Xj
 * @param B : reduction
 * @return double 
 */
static double getdistance(
        const vector<double>& Xi,
        const vector<double>& Xj,
        set<int>& B)
{
    double dis= 0.0;
    for(set<int>::iterator it = B.begin();
            it != B.end();
            ++it)
    {
        double tmpd = fabs(Xi[*it]-Xj[*it]);

        if(dis<tmpd){
            dis = tmpd;
        }
    }
    return dis;
}

/**
 * @brief calculate the positive region on reduction B
 * 
 * @param records : all the instance records 
 * @param x : instance x
 * @param B : reduction
 * @param ds : decision attribute location
 * @return double 
 */
double pos_x(
       const vector<vector<double> >& records,
       size_t x,
       set<int> & B,
       int ds)
{
    double mindes = 1;

    for(vector<vector<double> >::const_iterator icur = records.begin();
            icur != records.end();
            ++icur)
    {
        if(records[x][ds] != (*icur)[ds]){
	    
            double tmpdes  = getdistance(records[x],*icur,B);
            if (tmpdes < mindes){
                mindes = tmpdes;
            }
        }
    }

    return mindes;
}


/**
 * @brief calculate the dependency
 * 
 * @param dt : the decision table
 * @param B : the reduction
 * @param holdpos : if hold the positive region result
 * @param ... 
 * @return double 
 */
double Dependent(
        DecisionT& dt,
        set<int>& B,
        int holdpos,
        ...
        ){
    double sum_pos = 0.0;
    vector<vector<double> >& records = dt.getrecords();

    if (holdpos){ 
        va_list  vlist;
        va_start(vlist,holdpos);
        void* ppos = va_arg(vlist,void*);
        vector<double>& pos = *(vector<double>*)ppos;
        for(size_t it = 0;
                it < records.size();
                ++it)
        {
            double posx = pos_x(records,it,B,dt.getDs());

            if(pos.size()>records.size())
                pos[it] = posx;
            else
                pos.push_back(posx);
        }
    }else{
        int k;int sum=records.size();
        for(size_t it = 0;
                it < sum;
                ++it)
        {
            double posx = pos_x(records,it,B,dt.getDs());
            sum_pos += posx;
        }
    }

    double dep = sum_pos/(double)records.size();
    return dep;
}

/**
 * @brief delet n(=size) instance records from the records 
 * 
 * @param records 
 * @param size 
 */
void delrcds(vector<vector<double>> &records, int size){
    for(int i = 0 ;i < size ; i ++)
    {
        if(!records.empty()){
            records.pop_back();
        }
    }
}

/**
 * @brief update the positive region after some instances are deleted
 * 
 * @param records 
 * @param iold : the iterator of all the instances 
 * @param presize : the presize
 * @param idel : the iterator of the instances to be deleted
 * @param dsize : the size of the records to be deleted
 * @param B : the reduction
 * @param ds : the decision attribute location
 * @param prepos : the previos positive region
 */
void POSplus(
        vector<vector<double> > records,
        vector<vector<double> >::const_iterator iold,
        size_t presize,
        vector<vector<double> >::const_iterator idel,
        size_t dsize,
        set<int>& B,
        int ds,
        vector<double>& prepos
        ){
    
    vector<vector<double>> deled_records = records;
    delrcds(deled_records,dsize);

    /**
     * For each instance in the remained records, update the positive region; 
     */
    for(size_t i = 0;i<presize-dsize; i++){
        double mindes = prepos[i];
        for(vector<vector<double> >::const_iterator icur = idel;//traverse all the instances to be deleted
                icur != idel+dsize;
                ++icur){
            if((*(iold+i))[ds]!=(*icur)[ds]){  //Ensure that the instances are not of the same category
                double tmpdes  = getdistance(iold[i],*icur,B);//calulate  the new distance
                if ((tmpdes-mindes)<= abc){
                    /**
                     * If the lower approximation of the positive example is provided by the deleted sample, 
                     * it indicates the lower approximation need to be recalculated.
                    */
                    mindes = pos_x(deled_records,i,B,ds);
                    break;
                }
            }
            
        }
        prepos[i] = mindes;
    }
    for(int i = 0; i < dsize ; i ++){
        prepos.pop_back();
    }
}

