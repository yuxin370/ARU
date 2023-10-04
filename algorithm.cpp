#include <iostream>
#include <sstream>
#include <algorithm>
#include <cmath>
#include <float.h>
#include <set>
#include <string>
#include <vector>
#include <queue>
#include <fstream>
#include "algorithm.h"
int atr_size = 0;
double eps = 0.2;  
double step = 0.06; 
const double abc = 1e-8;
double alpha = 0.2; // the threshold of ARU
ofstream save;
ofstream svTime;
double pos_x(
    const vector<vector<double>> &records,
    size_t x,
    set<int> &B,
    int ds);

/**
 * @brief Output a two-dimensional vector
 * 
 * @param r the two-dimensional vector to be print
 */
void print2d(vector<vector<int>> &r)
{

    for (int i = 0; i < r.size(); i++)
    {
        for (int j = 0; j < r[0].size(); j++)
        {
            cout << r[i][j] << " : ";
        }
        cout << endl;
    }
}

/**
 * @brief Check if two values are equal, allowing for an error margin of 'eps'
 * 
 * @param a 
 * @param b 
 * @return true 
 * @return false 
 */
bool dequal(double a, double b)
{
    return abs((a - b) / a < eps || fabs(a - b) < abc);
}

/**
 * @brief print vector in given format
 * 
 * @param r 
 */
void printv(set<int> &r)
{
    for (set<int>::iterator it = r.begin();
         it != r.end();
         ++it)
    {
        cout << *it << ",";
    }
    cout << endl
         << "     size: " << r.size() << endl;
}


/**
 * @brief output the values of a vector to a file "save" in CSV format
 * 
 * @param r 
 */
void saveredu(set<int> &r){
    for (set<int>::iterator it = r.begin();
         it != r.end();
         ++it)
    {
        // if((it+1) !=r.end())
            save << *it << ",";
        // else
            // save << *it <<endl;
    }
    save<<endl;
}

/**
 * @brief check if the attribute is redudant.
 * 
 * @param Bi : the attribute
 * @param redu : the reduction  
 * @param d : the decision table
 * @param prepos : the previous positive region
 * @return true 
 * @return false 
 */
bool candel(int Bi, set<int> &redu, DecisionT &d, vector<double> &prepos)
{
    redu.erase(Bi);
    for (size_t x = 0; x < prepos.size(); x++)
    {
        double newpos = pos_x( //If the positive region remains unchanged after deletion, then it is redundant
            d.getrecords(),
            x,
            redu,
            d.getDs());
        if (newpos + alpha < prepos[x])
        {
            redu.insert(Bi);
            return false;
        }
    }
    // redu.insert(Bi);

    return true;
}

/**
 * @brief clear the redudancy
 * 
 * @param d : the decision table
 * @param redu : the reduction 
 * @param prepos : the previous positive region
 */
void clearRedundancy(
    DecisionT &d,
    set<int> &redu,
    vector<double> &prepos)
{
    set<int> B = redu;
    
    set<int>::iterator it = B.begin();
    for (; it != B.end();)
    {

        if (candel(*it, redu, d, prepos))
        {
            // redu.erase(*it);
        }
        ++it;
    } 
}


/**
 * @brief Find the attribute that results in the maximum increase in the positive region for the ARU
 * 
 * @param d : the decision table
 * @param redu : the reduction
 * @param ds :the class label location
 * @param lef : the left attribute
 * @param prepos : the previous positive region
 * @param allpos : the positive region on all the attribute
 * @param a : the best attribute
 * @param conins : the unlearning set
 * @param t 
 */
void findbest(
    DecisionT d,
    set<int> &redu,
    int ds,
    set<int> &lef,
    vector<double> &prepos,
    vector<double> &allpos,
    int &a,
    vector<int> &conins,
    vector<int> &t)
{


    clock_t start0, finish0;
    double totaltime;
    double max = 0;
    set<int>::iterator i = lef.begin();
    a = *i;

    vector<int> tmp;
    bool flag = true;

    for (; i != lef.end(); ++i) //traverse all the left attribute
    {
        tmp.clear();
        int tmpa = *i;
        redu.insert(tmpa); //add tmpa into the reduction
        double count = 0;
        for (int j = 0; j < conins.size(); j++) //traverse the instance in the unlearning set
        {
            
            double posj = pos_x(
                d.getrecords(),
                conins[j],
                redu,
                ds);

            count += (posj-prepos[conins[j]]);
        }
        if (count > max)
        { 
            /*
             *  Compare, if this relative importance is greater than the current maximum, 
             *  then update it
            */
            max = count;
            a = tmpa;
        }
        redu.erase(tmpa); //rollback
    }
    redu.insert(a);
    lef.erase(a);
}

/**
 * @brief Find the attribute that results in the maximum increase in the positive region for the PAR
 * 
 * @param d : the decision table
 * @param redu : the reduction
 * @param ds :the class label location
 * @param lef : the left attribute
 * @param prepos : the previous positive region
 * @param allpos : the positive region on all the attribute
 * @param a : the best attribute
 * @param t 
 */
void findbest2(
    DecisionT d,
    set<int> &redu,
    int ds,
    set<int> &lef,
    vector<double> &prepos,
    vector<double> &allpos,
    int &a,
    // vector<int> &conins,
    vector<int> &t)
{
    clock_t start0, finish0;
    double totaltime;
    double max = 0;
    set<int>::iterator i = lef.begin();
    a = *i;


    for (; i != lef.end(); ++i) //traverse all the left attribute
    {
        int tmpa = *i;
        redu.insert(tmpa); //add tmpa into the reduction
        double count = 0;
        for (int j = 0; j < d.getrecords().size(); j++) //traverse all the atribute
        {
            
            double posj = pos_x(
                d.getrecords(),
                size_t(j),
                redu,
                ds);
            if(prepos[j] + alpha < allpos[j])
                count += (posj - prepos[j]);
        }
        if (count > max)
        {
            /*
             *  Compare, if this relative importance is greater than the current maximum, 
             *  then update it
            */
            max = count;
            a = tmpa;
        }
        redu.erase(tmpa); //rollback
    }
    redu.insert(a);
    lef.erase(a);
}

int a = 0;

/**
 * @brief check if the current attribute set is a redution(reach maxmun positive region) for PAR 
 * 
 * @param curpos : current positive region
 * @param allpos : previous positive region
 * @return true 
 * @return false 
 */
bool reach_maxmun(vector<double> curpos,vector<double> allpos){

    for (size_t i = 0; i < curpos.size(); i++)
    {
        if (curpos[i] + alpha < allpos[i])
        { 
            return false;
        }
    }
    return true;

}




/**
 * @brief the algorithm for PAR 
 * 
 * @param d the decision table
 * @param redu the reduction
 */
void car_dep_par(
    DecisionT &d,
    set<int> &redu)
{
    clock_t start0, finish0;
    double totaltime;
    set<int> lef;
    set<int> &Cs = d.getCs();

    set_difference(
        Cs.begin(),
        Cs.end(),
        redu.begin(),
        redu.end(),
        std::inserter(lef, lef.begin()));
        
    auto &records = d.getrecords();
    

    vector<double> allpos;
    Dependent(d, d.getCs(), 1, &allpos); //calculate the lower approximate for each instances on all the attributes
    vector<double> prepos;
    Dependent(d, redu, 1, &prepos); //calculate the lower approximate for each instances on the previous reduction
    vector<int> conins, t, tmp;

    while (!reach_maxmun(prepos,allpos))
    {
        if(redu.size()==atr_size)
            break;
        int i = -1;
        start0 = clock();
        findbest2(d, redu, d.getDs(), lef, prepos,allpos, i, t); //find a attribute with the maximum importance
        finish0 = clock();
        totaltime = (double)(finish0 - start0) / CLOCKS_PER_SEC; //calculate thr time of "findbest"
        prepos.clear();
        Dependent(d, redu, 1, &prepos); 

    }

    clearRedundancy(d, redu, allpos); //clear redudancy
}


/**
 * @brief the algorithm for ARU
 * 
 * @param d : the decision table
 * @param ib :
 * @param dsize : the loop size
 * @param redu : the reduction
 * @param prepos : the positive region on previous reduction
 * @param allpos : the positive region on all attribute
 */
void car_dep_plus(
    DecisionT &d,
    vector<vector<double>>::const_iterator ib,
    size_t dsize,
    set<int> &redu,
    vector<double> &prepos,
    vector<double> &allpos)
{
    clock_t start0, finish0;
    double totaltime;
    set<int> lef;
    set<int> &Cs = d.getCs();

    set_difference(
        Cs.begin(),
        Cs.end(),
        redu.begin(),
        redu.end(),
        std::inserter(lef, lef.begin()));
    auto &records = d.getrecords();

    start0 = clock();

    //update the positive region on all attribute
    POSplus(
        records,
        records.begin(),
        records.size(),
        ib,
        dsize,
        d.getCs(),
        d.getDs(),
        allpos); 


    //update the positive region on the previos reduction
    POSplus(
        records,
        records.begin(),
        records.size(),
        ib,
        dsize,
        redu,
        d.getDs(),
        prepos); 

    
    finish0 = clock();
    totaltime = (double)(finish0 - start0) / CLOCKS_PER_SEC; 
    d.delrecords(dsize);

    vector<int> conins, t, tmp;

    for (size_t i = 0; i < prepos.size(); i++)
    {
        if (prepos[i] + alpha < allpos[i])
        { 
            //if the postive region on the previous reduction < the postive region on all attribute
            // the instance is a member of unlearning set.
            conins.push_back(i);
        }
    }

    int s_size = conins.size(); //size of unlearning set

    while (s_size > 0)
    {
        if(redu.size()==atr_size)
            break;
        int i = -1;
        start0 = clock();
        findbest(d, redu, d.getDs(), lef, prepos, allpos, i, conins, t); //choose a attribute with the maximum importance
        finish0 = clock();
        totaltime = (double)(finish0 - start0) / CLOCKS_PER_SEC; //calculate the time for "findbest"

        //update the unlearning set.
        prepos.clear();
        Dependent(d, redu, 1, &prepos);
        conins.clear();
        for (size_t i = 0; i < prepos.size(); i++)
        {
            if (prepos[i] + alpha < allpos[i])
            { 
                conins.push_back(i);
            }
        }
        s_size = conins.size();
      

    }


    int size = redu.size();
    clearRedundancy(d, redu, allpos); //clear redudancy
    if(redu.size()<size) Dependent(d,redu,1,&prepos);
}

/**
 * @brief expiriment procedure for ARU
 * 
 * @param T : decision table
 * @param records : data records
 * @param initsize : initial size of the records
 * @param dsize : The number of instances deleted in each iteration
 * @param redu : reduction
 */
void icar_dep(
    DecisionT &T,
    vector<vector<double>> &records,
    size_t initsize,
    size_t dsize,
    set<int> &redu)
{

    clock_t start, finish;
    double totaltime;
    start = clock();
    vector<vector<double>>::const_iterator curit = records.begin();
    T.addrecords(curit, initsize);
    car_dep_par(T, redu); //initialize, calculate the initial reduction on the old reductoin


    // initialize, calculate the positive region on the previous reduction
    vector<double> curpos;
    Dependent(T, redu, 1, &curpos); 

    // initialize, calculate the positive region on all the attributes.
    vector<double> allpos;
    Dependent(T, T.getCs(), 1, &allpos); 
    finish = clock();
    totaltime = (double)(finish - start) / CLOCKS_PER_SEC; //calculate time of initialization

    printv(redu);
    saveredu(redu);

    cout << "loop 1, totaltime: " << totaltime << endl
         << endl;
    svTime<<totaltime<<",";

    curit = curit + initsize; 
    int index = 2;
    // loop, simulating the unlearning process.
    while (curit - records.begin() > 2 * dsize) 
    {

        start = clock();
        car_dep_plus(
            T,
            (curit - dsize > records.begin() ? curit - dsize : records.begin()),
            (curit - dsize > records.begin() ? dsize : curit - records.begin()),
            redu,
            curpos,
            allpos); 
        finish = clock();
        totaltime = (double)(finish - start) / CLOCKS_PER_SEC; 

        curit = curit - dsize > records.begin() ? curit - dsize : records.begin();

        printv(redu);
        saveredu(redu);
        cout << "loop " << index << " ,totaltime: " << totaltime << endl
             << endl;
        index++;
        svTime<<totaltime<<",";
    }
}


/**
 * @brief expiriment procedure for Naive-U
 * 
 * @param T : decision table
 * @param records : data records
 * @param initsize : initial size of the records
 * @param dsize : The number of instances deleted in each iteration
 * @param redu : reduction
 * @param finaldep 
 */
void _noncar_dep_unuse(
    DecisionT &T,
    vector<vector<double>> &records,
    size_t initsize,
    size_t dsize,
    set<int> &redu,
    double finaldep)
{
    clock_t start, finish;
    double totaltime;
    vector<vector<double>>::const_iterator curit = records.begin();
    T.addrecords(curit, initsize);
    start = clock();
    car_dep_par(T, redu);
    finish = clock();
    totaltime = (double)(finish - start) / CLOCKS_PER_SEC;

    printv(redu);
    saveredu(redu);
    cout << "loop 1, totaltime: " << totaltime << endl
         << endl;
    svTime<<totaltime<<",";
    curit = curit + initsize;
    int index = 2;

    

    while (curit - records.begin() > 2 * dsize)
    {
        // redu.clear();
       
        dsize = curit - dsize > records.begin() ? dsize : curit - records.begin();
        T.delrecords(dsize);
        
        
        start = clock();
        car_dep_par(T, redu);

        finish = clock();
        totaltime = (double)(finish - start) / CLOCKS_PER_SEC;

        curit = curit - dsize > records.begin() ? curit - dsize : records.begin();

        printv(redu);
        saveredu(redu);
        cout << "loop " << index << " ,totaltime: " << totaltime << endl
             << endl;
        index++;
        svTime<<totaltime<<",";
    }
}


/**
 * @brief expiriment procedure for PAR
 * 
 * @param T : decision table
 * @param records : data records
 * @param initsize : initial size of the records
 * @param dsize : The number of instances deleted in each iteration
 * @param redu : reduction
 * @param finaldep 
 */
void _icar_dep_unuse(
    DecisionT &T,
    vector<vector<double>> &records,
    size_t initsize,
    size_t dsize,
    set<int> &redu,
    double finaldep)
{
    clock_t start, finish;
    double totaltime;
    //initialize
    vector<vector<double>>::const_iterator curit = records.begin();
    T.addrecords(curit, initsize);
    start = clock();
    car_dep_par(T, redu);
    finish = clock();
    totaltime = (double)(finish - start) / CLOCKS_PER_SEC;

    printv(redu);
    saveredu(redu);
    cout << "loop 1, totaltime: " << totaltime << endl
         << endl;
    svTime<<totaltime<<",";
    
    curit = curit + initsize;
    int index = 2;


    while (curit - records.begin() > 2 * dsize)
    {


        redu.clear();
        dsize = curit - dsize > records.begin() ? dsize : curit - records.begin();
        T.delrecords(dsize);

       
        start = clock();
        car_dep_par(T, redu);

        finish = clock();
        totaltime = (double)(finish - start) / CLOCKS_PER_SEC;

        curit = curit - dsize > records.begin() ? curit - dsize : records.begin();

        printv(redu);
        saveredu(redu);
        cout << "loop " << index << " ,totaltime: " << totaltime << endl
             << endl;
        index++;
        svTime<<totaltime<<",";
    }
}

/**
 * @brief input data and revoke the algorithm ARU
 * 
 * @param datafile : datafile
 * @param cons : attribute number 
 * @param N : attribute size
 * @param ds : the decision attribute(class label) location
 * @param allsize : the size of all the attribute
 * @param initsize : the size of the initial records size
 * @param dsize : The number of instances deleted in each iteration
 * @param afa :
 * @param k : fold number
 * @param savepath : the path of result file
 * @param timepath  : the path of the execution time result file
 */
void test2(
    string datafile,
    int *cons,
    size_t N,
    int ds,
    size_t allsize,
    size_t initsize,
    size_t dsize,
    double afa,
    int k,
    string savepath,
    string timepath)
{
    atr_size = N;
    svTime.open(timepath,ios::out|ios::app);
    svTime<<"ARU"<<endl;
    cout << "##################################### ARU ####################################" << endl;
    alpha = afa;
    clock_t start, finish;
    double totaltime;
    size_t k_size = allsize - initsize;
    save.open(savepath,ios::out);

    ifstream in;
    in.open(datafile.c_str());
    if (!in.is_open())
    {
        cout << "open file fail!" << endl;
        return;
    }
    vector<vector<double>> init_records;
    string line;
    while (getline(in, line))
    {
        istringstream stream(line);

        vector<double> *record = new vector<double>();
        string field;
        while (getline(stream, field, ','))
        {
            record->push_back(atof(field.c_str())); 
        }
        init_records.push_back(*record);
    }
    in.close();

    for(int i = 0;i < k ;i++){

        cout<<endl<<"### Fold "<<(i+1)<<" ###"<<endl;
        save << "Fold "<<(i+1)<<endl;
        vector<vector<double>> records(init_records);
        if((i+1)*k_size <= allsize){
            records.erase(records.begin() + k_size*i,records.begin() + (i+1)*k_size);
        }else{
            initsize = allsize - k_size;
            records.erase(records.begin() + k_size*i,records.end());
        }
        set<int> redu;
        redu.clear();
        vector<int> cs(cons, cons + N);

        DecisionT _T(cs, ds);

        //excute the algorithm
        start = clock();
        icar_dep(_T, records, initsize, dsize, redu);
        finish = clock();

        //output the time cost
        totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
        cout << "--> " << totaltime << " <--" << endl;
        printv(redu);
        cout << "<--" << endl
            << endl;
        svTime<<endl;
    }
    save.close();
    svTime.close();
}

/**
 * @brief input data and revoke the algorithm PAR
 * 
 * @param datafile : datafile
 * @param cons : attribute number 
 * @param N : attribute size
 * @param ds : the decision attribute(class label) location
 * @param allsize : the size of all the attribute
 * @param initsize : the size of the initial records size
 * @param dsize : The number of instances deleted in each iteration
 * @param afa :
 * @param k : fold number
 * @param savepath : the path of result file
 * @param timepath  : the path of the execution time result file
 */
void test1(
    string datafile,
    int *cons,
    size_t N,
    int ds,
    size_t allsize,
    size_t initsize,
    size_t dsize,
    double afa,
    int k,
    string savepath,
    string timepath)
{
    atr_size = N;
    svTime.open(timepath,ios::out | ios::app);
    svTime<<"PAR"<<endl;
    cout << "##################################### PAR ####################################" << endl;

    clock_t start, finish;
    double totaltime;
    size_t k_size = allsize - initsize;
    save.open(savepath,ios::out);
    alpha = afa;
    // read data 
    ifstream in;
    in.open(datafile.c_str());
    if (!in.is_open())
    {
        cout << "open file fail!" << endl;
        return;
    }
    vector<vector<double>> init_records;
    string line;
    while (getline(in, line))
    {
        istringstream stream(line);

        vector<double> *record = new vector<double>();
        string field;
        while (getline(stream, field, ','))
        {
            record->push_back(atof(field.c_str()));
        }
        init_records.push_back(*record);
    }
    in.close();

    for(int i = 0; i < k ;i++){
        cout<<endl<<"### Fold "<<(i+1)<<" ###"<<endl;
        save << "Fold "<<(i+1)<<endl;
        
        vector<vector<double>> records(init_records);
        if((i+1)*k_size <= allsize){
            records.erase(records.begin() + k_size*i,records.begin() + (i+1)*k_size);
        }else{
            initsize = allsize - k_size;
            records.erase(records.begin() + k_size*i,records.end());
        }


        set<int> redu;
        vector<int> cs(cons, cons + N);
        DecisionT _T2(cs, ds);
        double finaldep = 0;

        start = clock();

        //execute the algorithm
        _icar_dep_unuse(_T2, records, initsize, dsize, redu, finaldep);
        finish = clock();

        //output the time cost
        totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
        cout << "--> " << totaltime << " <--" << endl;
        printv(redu);
        cout << "<--" << endl
        << endl;
        svTime<<endl;
    }
    save.close();
    svTime.close();

}

/**
 * @brief input data and revoke the algorithm Naive-U
 * 
 * @param datafile : datafile
 * @param cons : attribute number 
 * @param N : attribute size
 * @param ds : the decision attribute(class label) location
 * @param allsize : the size of all the attribute
 * @param initsize : the size of the initial records size
 * @param dsize : The number of instances deleted in each iteration
 * @param afa :
 * @param k : fold number
 * @param savepath : the path of result file
 * @param timepath  : the path of the execution time result file
 */
void test0(
    string datafile,
    int *cons,
    size_t N,
    int ds,
    size_t allsize,
    size_t initsize,
    size_t dsize,
    double afa,
    int k,
    string savepath,
    string timepath)
{
    atr_size = N;
    svTime.open(timepath,ios::out|ios::app);
    svTime<<"Naive-U"<<endl;
    cout << "##################################### Naive-U ####################################" << endl;

    clock_t start, finish;
    double totaltime;
    size_t k_size = allsize - initsize;
    save.open(savepath,ios::out);
    alpha = afa;
    // read data 
    ifstream in;
    in.open(datafile.c_str());
    if (!in.is_open())
    {
        cout << "open file fail!" << endl;
        return;
    }
    vector<vector<double>> init_records;
    string line;
    while (getline(in, line))
    {
        istringstream stream(line);

        vector<double> *record = new vector<double>();
        string field;
        while (getline(stream, field, ','))
        {
            record->push_back(atof(field.c_str()));
        }
        init_records.push_back(*record);
    }
    in.close();

    for(int i = 0; i < k ;i++){
        cout<<endl<<"### Fold "<<(i+1)<<" ###"<<endl;
        save << "Fold "<<(i+1)<<endl;
        
        vector<vector<double>> records(init_records);
        if((i+1)*k_size <= allsize){
            records.erase(records.begin() + k_size*i,records.begin() + (i+1)*k_size);
        }else{
            initsize = allsize - k_size;
            records.erase(records.begin() + k_size*i,records.end());
        }


        set<int> redu;
        vector<int> cs(cons, cons + N);
        DecisionT _T2(cs, ds);
        double finaldep = 0;

        start = clock();

        //excute the algorithm
        _noncar_dep_unuse(_T2, records, initsize, dsize, redu, finaldep);
        finish = clock();

        //output the time cost
        totaltime = (double)(finish - start) / CLOCKS_PER_SEC;
        cout << "--> " << totaltime << " <--" << endl;
        printv(redu);
        cout << "<--" << endl
        << endl;
        svTime<<endl;
    }
    save.close();
    svTime.close();

}
