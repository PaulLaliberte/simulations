#include "dp.hpp"

DP::DP() : StateInd() {}

void DP::setParams(int ed, int icu, int a1, int t, double l1, double mu, double phi1, double h1, double beta)
{
    //NOTE: inherited function must match base class function definition to overload correctly
    N_ED = ed;
    N_ICU = icu;
    N_A1 = a1;
    N_state = (N_ICU + 1) * (N_ED + 1);
    T = t;
    MU = mu;

    //cost parameters
    BETA = beta;
    PHI1 = phi1;
    H1 = h1;
    iter_tol = 1e-8;

    //add elements to lambda array
    lambda.push_back(l1);

    //state/action space ane viable actions
    states = list_states();
    actions = list_actions();

    //arrivals
    arrivals = list_arrivals();

    //construct condition poisson pdf vector before-hand
    condPoiPdf = get_condPoiPdf(lambda.at(0), N_A1);

    transition_index();
}

bool DP::viable_action(int s, int d)
{
    //Last condition is to keep the state space from overflowing
    if (actions(d, 0) <= N_ICU - states(s, 1) && actions(d, 0) + actions(d, 1) <= states(s, 0))
    {
		return true;
	}
    else
    {
        return false;
    }
}

vec DP::get_condPoiPdf(double arrival_rate, int arrival_cap)
{
    vec pdfs = zeros<vec>(arrival_cap+1);

    //unrolled recursion (arrival cap small --- will not having any floating point errors)
    for (int i = 0; i < arrival_cap; i++)
    {
        pdfs[i] = (std::pow(arrival_rate, i) * std::exp(-1 * arrival_rate)) / boost::math::factorial<double>((double) i);
    }

    double tmp_sum = 0.0;
    for (auto &x : pdfs)
    {
        tmp_sum += x;
    }

    //1 - sum_{i=0}^{A-1}
    pdfs[arrival_cap] = 1 - tmp_sum;


    //check that probability sums to 1
    double total_sum = tmp_sum + pdfs[arrival_cap];
    if (total_sum != 1.0)
    {
        std::cout << total_sum << std::endl;
        throw std::invalid_argument("Conditional Poisson Pdf does not sum to 1.");
    }

    return pdfs;
}

void DP::transition_index()
{
    for (int s = 0; s < N_state; s++)
    {
        tuple_ind ti (states(s,0), states(s,1));
        new_T_index[ti] = s;
    }
}

// Transition probability vector Pr( . | s1, d)
vec DP::transition_prob(int s1, int d)
{
	vec Tp = zeros<vec>(N_state);

    int s1_ed = states(s1, 0) - actions(d, 0) - actions(d, 1);
    int s1_icu = states(s1, 1) + actions(d, 0);

    for (int arr = 0; arr <= N_A1; arr++)
    {
        for (int dept = 0; dept <= s1_icu; dept++)
        {
            int s2_ed = std::min(s1_ed + arr, N_ED);
            int s2_icu = s1_icu - dept;
            tuple_ind ti (s2_ed, s2_icu);
            int index_s2 = new_T_index[ti];
            Tp(index_s2) += condPoiPdf(arr) * nchoosek(states(index_s2, 1) + dept, dept) * std::pow(MU, dept) * std::pow(1 - MU, states(index_s2 , 1));
        }
    }
    
    double total_sum = 0.0;
    for (int i = 0; i < N_state; i++)
    {
        total_sum += Tp(i);
    }

    //There is some floating point error given certain transition states
    if (1.0 - total_sum > 1e-8)
    {
        std::cout << std::setprecision(30) << 1.0 - total_sum << std::endl;
        throw std::invalid_argument("Transition probability does not sum to 1.");
    }
        
    return Tp;
}

mat DP::prob_action(vec value_function)
{
	mat P = zeros<mat>(N_state, N_action);

	for (int i = 0; i < N_state; i++)
	{
		for (int j = 0; j < N_action; j++)
		{
			if (viable_action(i, j) == true)
			{
                //-c(s,d)
                double cost = -(PHI1 * actions(j, 1) + H1 * (states(i,0) - actions(j,0) - actions(j,1)));

                //p(s' | s,d)      
                vec Tp = transition_prob(i, j);

                //exp(-c(s,d) + beta*sum_s' p(s' |s, d)v(s')
                P(i, j) = std::exp(cost + BETA * dot(Tp, value_function));
            }
        }
	}

    colvec row_sum = sum(P, 1);
    for (int i = 0; i < N_state; i++)
    {
        for (int j = 0; j < N_action; j++)
        {
            P(i,j) = P(i,j) / row_sum(i);
        }
    }

    /* -- Overkill --
    //Check that all rows sum to 1
    colvec row_sum_check = sum(P, 1);
    for (int i = 0; i < N_state; i++)
    {
        if (1.0 - row_sum_check(i) > 1e-8)
        {
            std::cout << "Row: " << i << std::endl;
            throw std::invalid_argument("Choice probabilities do not all sum to 1.");
        }
    }
    */

    return P;
}

vec DP::value_iteration()
{
    double delta_v = 1000;
    int k=0;

    vec V_prev = zeros<vec>(N_state);
    vec V_bar = zeros<vec>(N_state);

    while (delta_v > iter_tol)
    {
        //get previous value_function
        V_prev = V_bar;
        k++;

        mat P = prob_action(V_prev);

        for (int i = 0; i < N_state; i++)
        {
            double tmp = 0.0;

            for (int j = 0; j < N_action; j++)
            {
                if (viable_action(i,j) == true)
                {
                    //-c(s,d)
                    double cost = -(PHI1 * actions(j, 1) + H1 * (states(i,0) - actions(j,0) - actions(j,1)));

                    //p(s' | s,d)  
                    vec Tp = transition_prob(i, j);

                    //sum_{d in π(s)} P(d|s) * (-c(s,d) + beta * sum_s' p(s' | s,d)v(s')
                    double M = cost + BETA * dot(Tp, V_prev);
                    tmp += P(i,j) * M;
                }
            }
            V_bar(i) = tmp;
        }
        delta_v = max(abs(V_bar - V_prev));
    }
    return V_bar;
}

int DP::choose_action(int s, double r, mat P)
{
    //Will throw error -- an action should always be chosen.
    int row_ind = -1;

    double upto = 0.0;
    for (int d = 0; d < N_action; d++)
    {
        if (P(s,d) + upto >= r)
        {
            row_ind = d;
            break;
        }
        else
        {
            upto += P(s,d);
        }
    }
    return row_ind;
}

void DP::write_outputs(int *s_indices, int *a_indicies, char *icup_c, char *rrp_c, char *rho, char *counter, char *beta)
{
    //calculate average stage cost
    double total_cost = 0;
    for (int i = 0; i < stage_cost.size(); i++)
    {
        total_cost += stage_cost.at(i);
    }

    // Print out the state-action pairs to file
	char stateCostFile[40];
    sprintf(stateCostFile, "stage_cost/full_tmp_6/beta_%s/%s/states_beta_%s.txt", beta, rho, counter);
    std::ofstream outsc(stateCostFile);

    if (total_cost == 0)
    {
        outsc << total_cost << endl;
    }
    else
    {
        outsc << total_cost / stage_cost.size() << endl;
    }
	outsc.close();

    // Print out the state-action pairs to file
	char stateFile[40];
    sprintf(stateFile, "states/full_tmp_6/beta_%s/%s/states_beta_%s.txt", beta, rho, counter);
    std::ofstream outs(stateFile);
	for (int t = 0; t < T; t++)
	{
		outs << states(s_indices[t], 0) << " " << states(s_indices[t], 1) << endl;
	}
	outs.close();

	char actionFile[40];
	sprintf(actionFile, "actions/full_tmp_6/beta_%s/%s/actions_beta_%s.txt", beta, rho, counter);
    std::ofstream outd(actionFile);
	for (int t = 0; t < T; t++)
	{
		outd << actions(a_indicies[t], 0) << " " << actions(a_indicies[t], 1) <<  endl;
	}
	outd.close();

	char saFile[40];
	sprintf(saFile, "sa_pairs/full_tmp_6/beta_%s/%s/sa_pairs_beta_%s.txt", beta, rho, counter);
    std::ofstream outsa(saFile);
	for (int t = 0; t < T; t++)
	{
		outsa << s_indices[t] << " " << a_indicies[t] << endl;
	}
	outsa.close();
}

