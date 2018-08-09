#include "base_model.cpp"

int main(int argc, char *argv[])
{
    StateInd newSim;

    newSim.setParams(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]), atoi(argv[4]), atof(argv[5]), atof(argv[6]), atof(argv[7]), atof(argv[8]), atof(argv[9]));

    //newSim.getParams();
    
    //Start Simulation
    vec V_bar_tmp = newSim.value_iteration();
    mat P_action_0 = newSim.prob_action(V_bar_tmp);
    
	// Calculate the states
    int T13 = atoi(argv[4]) + 6 * 30 * 12;
	int *s_index_all = new int[T13]; // index the states for convenience
	int *d_index_all = new int[T13]; // index the actions
	
    //Uniform random variable generator
    std::default_random_engine generatorU;
    std::uniform_real_distribution<double> distribution(0.0,1.0);

	s_index_all[0] = 0; // start from state (0,0,0) with index 0
	mat trial_states = zeros<mat>(T13, 2);

	for (int t = 0; t < T13 - 1; t++)
	{
		// Simulate the actions taken on current period state
		double r = distribution(generatorU);
        d_index_all[t] = newSim.choose_action(s_index_all[t], r, P_action_0);
        
        // Simulate how many departs from the ICU in the current period (can be 0)
        int dept = newSim.create_departures(trial_states(t, 1), d_index_all[t]);

        // Update the state at the beginning of t+1, before any actions in t+1
        trial_states(t+1, 0) = newSim.update(trial_states(t, 0), t, d_index_all[t], dept, true);
        trial_states(t+1, 1) = newSim.update(trial_states(t, 1), t, d_index_all[t], dept, false);

        //std::cout << trial_states.row(t+1) << std::endl;

		int index = newSim.find_index(trial_states(t + 1, 0), trial_states(t + 1, 1));
		if (index >= 0) s_index_all[t + 1] = index;
	}

	double r = distribution(generatorU);
	d_index_all[T13 - 1] = newSim.choose_action(s_index_all[T13 - 1], r, P_action_0);

	// Discard the warm-up period
	int *s_index = s_index_all + 6 * 30 * 12;
	int *d_index = d_index_all + 6 * 30 * 12;

    newSim.write_outputs(s_index, d_index, argv[7], argv[8], argv[10], argv[11], argv[9]);
	
	return 0;
}
