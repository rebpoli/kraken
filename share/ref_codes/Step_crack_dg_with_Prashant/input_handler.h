/*
 * parameter_input.h
 *
 *  Created on: Apr 28, 2015
 *      Author: mital
 */
#ifndef INPUT_HANDLER_H_
#define INPUT_HANDLER_H_

#include <deal.II/base/parameter_handler.h>
#include <list>
#include <iostream>
#include <fstream>

using namespace dealii;

class InputHandler {
public:
	InputHandler (const int argc, char *const *argv) :
		flag_success(false),
		sform (1),
		max_time_steps (500),
        max_iters_AL (500),
        decompose_stress_matrix(false),
        decompose_stress_rhs(false),
        penalty_DG(1),
        penalty_BC(10),
        pf_kappa(0.01),
        pf_epsilon(0.02),
        critical_Gc(2.7),
        lame_mu(81),
        lame_lambda(122),
		timestep (0.1),
		penalty_AL (1)
    {
		//Declare and Populate
		declare_parameters();
		parse_command_line(argc, argv);
    };

	int get_sform () const { return sform; }
    int get_maxtimesteps () const { return max_time_steps; }
    int get_maxitersAL () const { return max_iters_AL; } // max_penal_iterations
	bool get_decomposestressmatrix () const { return decompose_stress_matrix; } 
    bool get_decomposestressrhs () const { return decompose_stress_rhs; } 
    double get_penaltyDG () const { return penalty_DG; } // penalty_u_face
	double get_penaltyBC () const { return penalty_BC; } // penalty_u_bc
    double get_pfkappa () const { return pf_kappa; } // constant_k
    double get_pfepsilon () const { return pf_epsilon; } // alpha_eps
    double get_criticalGc () const { return critical_Gc; } // Gc
    double get_lamemu () const { return lame_mu; } // lame_coefficient_mu
    double get_lamelambda () const { return lame_lambda; } // lame_coefficient_lambda
    double get_timestep () const { return timestep; } 
    double get_penaltyAL () const { return penalty_AL; } // penal_gamma
    std::string get_testcase () const { return testcase; } // mesh input filename
    std::string get_lintype () const { return lintype; } // linearization type
    std::string get_opttype () const { return opttype; } // optimization type 
	//void echo_parameters () const;

private:
	void print_usage_message ();
	void declare_parameters ();
	void parse_command_line (const int argc, char *const *argv);

	ParameterHandler prm;
	bool flag_success;
	int sform, max_time_steps, max_iters_AL;
    bool decompose_stress_matrix, decompose_stress_rhs;
	double penalty_DG, penalty_BC, pf_kappa, pf_epsilon, critical_Gc, lame_mu, lame_lambda,
           timestep, penalty_AL;
    std::string testcase, lintype, opttype;
};

#endif
