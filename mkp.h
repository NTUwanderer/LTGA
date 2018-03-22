#ifndef _mkp_h_
#define _mkp_h_

struct MKPinstance 
{
    int var;
    int constr;
    int opt;

    int** w;
    int* capacity, *profit;
};

double evaluateMKP(int*, MKPinstance*);
void loadMKP(char*, MKPinstance*);
void freeMKPinstance(MKPinstance*);

#endif
