#include <iostream>

#include <string>
#include <typeinfo>
#include <stdexcept>
#include <memory>
#include <cstdlib>
#include <algorithm>
#include <cmath> 
#include <ctime>
#include <vector>
#include <armadillo>
#include <iomanip>
#include <random>
#include <fstream>
#include <sstream>

using namespace arma;

#ifndef __StateIndependent__
#define __StateIndependent__

class StateInd
{
    protected:
        //functions to create state/action space, viable actions, and arrivals
        imat list_states();
        imat list_actions();
        mat list_arrivals();

        //variables set on input
        int N_ED;
        int N_ICU;
        int N_A1;
        std::vector<double> lambda;
        int T;
        double MU;

        //matricies
        mat arrivals;

        //probability of actions
        double icup;
        double rrp;

        //Functions to ensure simulation running correctly
        void CheckProbRow(int s, mat P_action);


    public:
        int N_state;
        int N_action;

        imat states;
        imat actions;
        
        //functions:
        StateInd();
        ~StateInd();

        //set params
        virtual void setParams(int ed, int icu, int a1, int t, double l1, double mu, double icu_p, double reroute_p, double beta);
        void getParams();
        void getStateActionSpace();

        //simulation functions
        int nchoosek(int n, int k);
        int find_index(int n1, int n2);
        virtual bool viable_action(int s, int d);
        virtual vec value_iteration();
        virtual mat prob_action(vec value_function);
        int create_departures(int icu_s, int icu_a);
        virtual int choose_action(int s, double r, mat P_action);
        int update(int s, int t, int d, int dept, bool state_flag);
        virtual void write_outputs(int *s_indices, int *a_indicies, char *icup_c, char *rrp_c, char *rho, char *counter, char *beta);
};
                
#endif 
