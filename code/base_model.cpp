#include "base_model.hpp"

StateInd::StateInd() { /*DO NOT PASS PARAMETERS TO CONSTRUCTOR -- sub class inherits base construc too*/ }

// Index all possible states
imat StateInd::list_states()
{
	imat s = zeros<imat>(N_state, 2);
	int k = 0;
	for (int n1 = 0; n1 <= N_ED; n1++)
	{
		for (int n2 = 0; n2 <= N_ICU; n2++)
			{
				s(k, 0) = n1;
				s(k, 1) = n2;
				k++;
			}
	}

	return s;
}

// Index all possible actions
imat StateInd::list_actions()
{
    //Get number of actions needed
    int tmp = N_ED + 1;

    N_action = 0;
    for (int i = 0; i <= N_ED; i++)
    {
        N_action += tmp;
        tmp--;
    }

	imat d = zeros<imat>(N_action, 2);

    tmp = N_ED;
	int k = 0;
	for (int adm1 = 0; adm1 <= N_ED; adm1++)
	{
		for (int rerout1 = 0; rerout1 <= tmp; rerout1++)
		{
            d(k, 0) = adm1;
			d(k, 1) = rerout1;
			k++;
		}
        tmp--;
	}
    return d;
}

mat StateInd::list_arrivals()
{
    std::random_device r1;
	std::default_random_engine generator1(r1());
	std::poisson_distribution<int> distribution1(lambda[0]);

	int T13 = T + 6 * 30 * 12;  // Simulate 2 yrs + 6 months and treat the 6 months as warm up
	mat arrs = zeros<mat>(T13, 1);

	for (int t = 0; t < T13; t++) 
	{
		int a1 = distribution1(generator1);

		arrs(t, 0) = std::min(a1, N_A1);

	}
    return arrs;
}

void StateInd::setParams(int ed, int icu, int a1, int t, double l1, double mu, double icu_p, double reroute_p, double beta)
{
    N_ED = ed;
    N_ICU = icu;
    N_A1 = a1;
    N_state = (N_ICU + 1) * (N_ED + 1);
    T = t;
    MU = mu;

    icup = icu_p;
    rrp = reroute_p;

    //add elements to lambda array
    lambda.push_back(l1);

    //state/action space ane viable actions
    states = list_states();
    actions = list_actions();

    //arrivals
    arrivals = list_arrivals();
}

bool StateInd::viable_action(int s, int d)
{
    if (states(s,0) != 0)
    {
        int max_admit = states(s, 0);
        while (states(s,1) + max_admit > N_ICU)
        {
            max_admit--;
        }

        if (actions(d, 0) == max_admit)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
    else
    {
        if (actions(d, 0) == 0)
        {
            return true;
        }
        else
        {
            return false;
        }
    }
}

int StateInd::nchoosek(int n, int k)
{
    if (k == 0)
    {
        return 1;
    }
    else
    {
        return (n * nchoosek(n - 1, k - 1)) / k;
    }
}
    
// Find the index of which state
int StateInd::find_index(int n1, int n2)
{
	int k = 0;
	while (k < N_state)
	{
		if (states(k, 0) == n1 && states(k, 1) == n2) return k;
		k++;
	}

    throw std::invalid_argument("State overflow.");
}

vec StateInd::value_iteration()
{
    vec V_bar = zeros<vec>(N_state);
    return V_bar;
}

mat StateInd::prob_action(vec value_function)
{
	mat P = zeros<mat>(N_state, N_action);

	for (int i = 0; i < N_state; i++)
	{
		for (int j = 0; j < N_action; j++)
		{
			if (viable_action(i, j) == true)
			{
                int num_admit = actions(j, 0);
                int num_reroute = actions(j, 1);

                P(i,j) = 1;
                break;
            }
        }
	}
    return P;
}

int StateInd::create_departures(int icu_s, int icu_a)
{
    std::random_device rd;
    std::default_random_engine generatorRd(rd());

    int departure;
    if (icu_s - actions(icu_a, 0) > 0)
    {
        std::binomial_distribution<int> distribution4(icu_s + actions(icu_a, 0), MU);
	    departure = distribution4(generatorRd);
    }
    else departure = 0;
   
    return departure;
}

int StateInd::update(int s, int t, int d, int dept, bool state_flag)
{
    if (state_flag == true)
    {
        //Do not allow arrivals into an already full ED
        int tmp_arrivals = arrivals(t,0);
        while (tmp_arrivals > 0)
        {
            if (s + tmp_arrivals - actions(d,0) - actions(d,1) > N_ED)
            {
                tmp_arrivals -= 1;
            }
            else
            {
                break;
            }
        }
        return s + tmp_arrivals - actions(d,0) - actions(d,1);
    }
    else
    {
        //update ICU
        return s + actions(d,0) - dept;
    }
}

int StateInd::choose_action(int s, double r, mat P )
{
    int row_idx = 0;

    double upto = 0;

	for (int d = 0; d < N_action; d++)
	{
        if (P(s,d) + upto >= r)
        {            
            row_idx = d;
            break;
        }
        else
        {
            upto += P(s,d);
        }
    }
    return row_idx;
}

void StateInd::getParams()
{
    std::cout << "Emergency Department Size: " << N_ED << std::endl;
    std::cout << "ICU Size: " << N_ICU << std::endl;
    std::cout << "Max arrivals per period: " << N_A1 << std::endl;
    std::cout << "Time (2-hour periods): " << T << std::endl;

    std::cout << "Arrival Rates: "; 
    for(auto i = lambda.begin(); i != lambda.end(); ++i)
    {
        std::cout << *i << " ";
    }
    std::cout << std::endl;
    
    std::cout << "Service Rate: " << MU << std::endl;
}

void StateInd::getStateActionSpace()
{
    std::cout << "State Space: " << std::endl;
    std::cout << states << std::endl;

    std::cout << "Actions: " << std::endl;
    std::cout << actions << std::endl;
}

void StateInd::write_outputs(int *s_indices, int *a_indicies, char *icup_c, char *rrp_c, char *rho, char *counter, char *beta)
{
    //CheckActionsTaken(a_indicies, icup_c, rrp_c);
    
    // Print out the state-action pairs to file
	char stateFile[40];
    sprintf(stateFile, "states/full/beta_%s/%s/states_beta_%s.txt", beta, rho, counter);
    std::ofstream outs(stateFile);
	for (int t = 0; t < T; t++)
	{
		outs << states(s_indices[t], 0) << " " << states(s_indices[t], 1) << endl;
	}
	outs.close();

	char actionFile[40];
	sprintf(actionFile, "actions/full/beta_%s/%s/actions_beta_%s.txt", beta, rho, counter);
    std::ofstream outd(actionFile);
	for (int t = 0; t < T; t++)
	{
		outd << actions(a_indicies[t], 0) << " " << actions(a_indicies[t], 1) <<  endl;
	}
	outd.close();

	char saFile[40];
	sprintf(saFile, "sa_pairs/full/beta_%s/%s/sa_pairs_beta_%s.txt", beta, rho, counter);
    std::ofstream outsa(saFile);
	for (int t = 0; t < T; t++)
	{
		outsa << s_indices[t] << " " << a_indicies[t] << endl;
	}
	outsa.close();
}

void StateInd::CheckProbRow(int s, mat P_action)
{
    //Floating point problems arise -- must use this method for comparison
    colvec c = sum(P_action, 1);
    mat CT;
    CT.ones(c.size(), 1);

    bool same_vec = approx_equal(c.row(s), CT.row(s), "absdiff", 0.0000001);
    
    if (same_vec != true)
    {
        std::cout << sum(P_action, 1) << std::endl;
        throw std::invalid_argument("Row probabilities do not add to 1.");
    }    
}

StateInd::~StateInd() {}
