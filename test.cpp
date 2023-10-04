#include <iostream>
#include <string>
#include <sstream>
#include <queue>
#include <time.h>
#include <stdlib.h>
#include <iomanip>
#include "DecisionT.h"
#include "algorithm.h"

int main(int argc, char *argv[])
{
    //input the argvs
    int cs_start = atoi(argv[1]);
    int cs_end = atoi(argv[2]);
    int ds = atoi(argv[3]);
    string dataset = argv[4];
    int cons[cs_end - cs_start + 1];
    for (int i = cs_start; i <= cs_end; i++)
        cons[i - cs_start] = i;
    int cslen = cs_end - cs_start + 1;
    int allpart = atoi(argv[5]);
    double alpha = atof(argv[6]);
    int k = atoi(argv[7]);
    int initpart = (k - 1) * (allpart / k);
    int looppart = initpart / 6 + 10;

    string datafile = "standard//"+dataset+".csv";
    string timepath = "result\\" + dataset + "\\time_"+to_string(alpha)+".csv";

    for (int i = 2; i >= 0; i--)
    {
        if (i == 2)
        {
            // ARU
            string savepath = "result\\" + dataset + "\\aru_"+to_string(alpha)+".csv";

            test2(datafile,
                  cons,
                  cslen,
                  ds,
                  allpart,
                  initpart,
                  looppart,
                  alpha,
                  k,
                  savepath,
                  timepath);

        }
        else if (i == 1)
        {
            // PAR
            string savepath = "result\\" + dataset + "\\par_"+to_string(alpha)+".csv";

            test1(datafile,
                  cons,
                  cslen,
                  ds,
                  allpart,
                  initpart,
                  looppart,
                  alpha,
                  k,
                  savepath,
                  timepath);
        }else{
            //Naive-U
            string savepath = "result\\" + dataset + "\\naiveu_"+to_string(alpha)+".csv";

            test0(datafile,
                  cons,
                  cslen,
                  ds,
                  allpart,
                  initpart,
                  looppart,
                  alpha,
                  k,
                  savepath,
                  timepath);
	
        }
    }


    
    return 0;
}
