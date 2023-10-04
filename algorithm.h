#ifndef ALGORITHM_H
#define ALGORITHM_H
#include <vector>
#include "DecisionT.h"
using namespace std;


void printv(set<int>& r);

bool dequal(double a,double b);
double car_dep(
        DecisionT &d,
        set<int>& redu
        );
void icar_dep(
        DecisionT& T,
        vector<vector<double> > & records,
        size_t initsize,
        size_t gsize,
        set<int>& redu
        );
void test2(
        string datafile,
        int* cons,
        size_t N,
        int ds,
        size_t allsize,
        size_t initsize,
        size_t dsize,
	double dlpha,
        int k,
        string savepath,
        string timepath
        );
void test1 (
        string datafile,
        int* cons,
        size_t N,
        int ds,
        size_t allsize,
        size_t initsize,
        size_t dsize,
        double alpha,
        int k,
        string savepath,
        string timepath
        );
void test0(
        string datafile,
        int* cons,
        size_t N,
        int ds,
        size_t allsize,
        size_t initsize,
        size_t dsize,
	double dlpha,
        int k,
        string savepath,
        string timepath
        );

#endif
