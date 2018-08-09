#include "base_model.hpp"

#include <boost/tuple/tuple.hpp>
#include <boost/unordered_map.hpp>
#include <boost/math/special_functions/factorials.hpp>

#ifndef __DP__
#define __DP__

typedef boost::tuples::tuple<int, int> tuple_ind;

struct ihash
    : std::unary_function<tuple_ind, std::size_t>
{
    std::size_t operator()(tuple_ind const& e) const
    {
        std::size_t seed = 0;
        boost::hash_combine( seed, e.get<0>() );
        boost::hash_combine( seed, e.get<1>() );
        return seed;
    }
};

struct iequal_to
    : std::binary_function<tuple_ind, tuple_ind, bool>
{
    bool operator()(tuple_ind const& x, tuple_ind const& y) const
    {
        return ( x.get<0>()==y.get<0>() &&
                 x.get<1>()==y.get<1>());    }
};

typedef boost::unordered_map<tuple_ind, int, ihash, iequal_to> T_index;

class DP : public StateInd
{
    protected:
        //beta (discount factor)
        double BETA;

        //cost parameters
        double PHI1;
        double H1;
        double iter_tol;

        //vector to gather cost
        std::vector<double> stage_cost;

        vec condPoiPdf;

        T_index new_T_index;
        
    public:
        DP();
        
         //set params
         //NOTE: Instead of probabilities passed to set params, we now pass the cost parameters
         void setParams(int ed, int icu, int a1, int t, double l1, double mu, double phi1, double h1, double beta);
       
        bool viable_action(int s, int d);
        vec get_condPoiPdf(double arrival_rate, int arrival_cap);
        void transition_index();
        vec transition_prob(int s1, int d);
        mat prob_action(vec value_function);
        vec value_iteration();
        int choose_action(int s, double r, mat P_actions);
        void write_outputs(int *s_indices, int *a_indicies, char *icup_c, char *rrp_c, char *rho, char *counter, char *beta);

        
};
                
#endif 
