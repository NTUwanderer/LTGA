#include "mkp.h"
#include <iostream>
#include <fstream>

using namespace std;

double evaluateMKP(int *x, MKPinstance *inst) 
{
    // Penalty function
    // ref: https://pdfs.semanticscholar.org/88de/57c535141d7540baa7488d5d907cbef57771.pdf
    int sum[inst->var];
    for(int i = 0; i < inst->var; ++i)
        sum[i] = 0.0;
    
    for(int i = 0; i < inst->constr; ++i)
    {   
        double temp = 0;
        for(int j = 0; j < inst->var; ++j)
            if (x[j])
                temp += inst->w[i][j];
        
        if(temp > inst->capacity[i])
            for (int j = 0; j < inst->var; ++j)
                if (x[j])
                    sum[j] += inst->w[i][j];
    }
    
    int penalty = 0;
    for(int i = 0; i < inst->var; ++i)
        penalty += inst->profit[i] * sum[i];

    double reward = 0.0;
    for(int i = 0; i < inst->var; ++i)
        if(x[i])
            reward += inst->profit[i];

    if (penalty != 0)
        reward /= penalty;

    return reward;
}


void repairMKPSolu(int* x, MKPinstance *inst)
{
    /*  NEED LP SOLVER */
    int r[inst->constr];
    for(int j = 0; j < inst->var; ++j)
        if (x[j])
        {   
            for(int i = 0; i < inst->constr; ++i)
                r[i] += inst->w[i][j];
        }
    
}


void loadMKP(char* name, MKPinstance *inst)
{
    ifstream input;
    string line;

    input.open(name);
    
    if(!input) {
        cout << "Fail to open file: " << name << endl;
        exit(0);
    }

    // First line n, m, optimal
    input >> inst->var >> inst->constr >> inst->opt;

    inst->profit = new int[inst->var];

    for(int i = 0; i < inst->var; ++i)
        input >> inst->profit[i];

    inst->w = new int*[inst->constr];

    for(int i = 0; i < inst->constr; ++i) {   
        inst->w[i] = new int[inst->var];
        for(int j = 0; j < inst->var; ++j) {
            input >> inst->w[i][j];
        }
    }

    inst->capacity = new int[inst->constr];

    for(int i = 0; i < inst->constr; ++i) 
        input >> inst->capacity[i];

    input.close();
}


void freeMKPinstance(MKPinstance *inst)
{
    delete[] inst->capacity;
    delete[] inst->profit;
    
    for(int i = 0; i < inst->constr; ++i)
        delete[] inst->w[i];
    delete[] inst->w;
}
