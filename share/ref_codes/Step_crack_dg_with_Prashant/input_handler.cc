#include "input_handler.h"

void InputHandler::print_usage_message() {
	std::cout <<"Please use the command line to specify "
			"an input parameter file with a -p flag\n"
			"You may use the following structure for your "
			"parameter file\n";

	prm.print_parameters (std::cout, ParameterHandler::Text);
}

void InputHandler::declare_parameters() {

	prm.enter_subsection("Elasticity Parameters"); {
		prm.declare_entry("Lame's Parameter Mu", "80.77e+3",
				Patterns::Double(0, -1),
				"Lame's Second Parameter a.k.a Shear Modulus [MPa]");
		prm.declare_entry("Lame's Parameter Lambda", "121.15e+3",
				Patterns::Double(0, -1),
				"Lame's First Parameter [MPa]"); 
		prm.declare_entry("Critical Fracture Energy", "2.7",
				Patterns::Double(0, -1),
				"The elastic energy restitution rate for crack growth [N/mm]");
	}
	prm.leave_subsection();

	prm.enter_subsection("DG Parameters"); {
		prm.declare_entry("Face Term Penalty", "1",
				Patterns::Double(0, -1),
				"The penalty parameter for face terms in DG");
		prm.declare_entry("Boundary Term Penalty", "10",
				Patterns::Double(0, -1),
				"The penalty parameter for enforcing Dirichlet boundary conditions in DG");
		prm.declare_entry("sform", "1",
				Patterns::Integer(-1,1),
				"The sform determines which class of DG method will be used"
				"(+1:SIPG, 0:IIPG, -1:NIPG");
	}
	prm.leave_subsection();

	prm.enter_subsection("Phase-Field Parameters"); {
			prm.declare_entry("Constant K", "1e-12",
					Patterns::Double(0,1),
					"Energy regularization constant Kappa");
			prm.declare_entry("Constant E", "4.4e-2",
					Patterns::Double(0, 1),
					"Space regularization constant Epsilon [mm]");
            prm.declare_entry("Penalty Gamma", "1e-16",
                    Patterns::Double(0, -1),
                    "Constant penalty term from Augmented Lagrange method");
            prm.declare_entry("Maximum Penalty Iterations", "500",
                    Patterns::Integer(0, -1),
                    "Maximum allowable augmented lagrange iterations at a given timestep");
            prm.declare_entry("Linearization Scheme", "interpolation",
                    Patterns::Selection("simple|interpolation"),
                    "Specify how the semi-linear form is linearized");
            // Add stopping criterion tolerance
	}
	prm.leave_subsection();
	
    prm.enter_subsection("Runtime Parameters"); {
            prm.declare_entry("Optimization Type", "simple penalization",
                    Patterns::Selection("simple penalization|augmented lagrange"),
                    "Type of outer-loop optimization to be employed");
            prm.declare_entry("Test Case", "miehe tension",
                    Patterns::Selection("miehe tension|miehe shear|symmetric bending"),
                    "Name of the test case that is to be run");
			prm.declare_entry("Timestep Size", "5e-5",
					Patterns::Double(0,1),
					"Timestep size [s]");
			prm.declare_entry("Maximum Timesteps", "500",
					Patterns::Integer(0, 10000),
					"Number of timesteps to take");
            prm.declare_entry("Decompose Stress Matrix", "false",
                    Patterns::Bool(),
                    "True/Non-zero results in spectral decomposition of stress tensor in the matrix");
            prm.declare_entry("Decompose Stress RHS", "false",
                    Patterns::Bool(),
                    "True/Non-zero results in spectral decomposition of stress tensor in the RHS");
            // Add Newton Iteration tolerance
	}
	prm.leave_subsection();
}

void InputHandler::parse_command_line(const int argc, char *const *argv) {
	if (argc <2) {
		print_usage_message ();
		exit (1);
	}

	std::list<std::string> args;
	for (int i=1; i<argc; ++i)
		args.push_back (argv[i]);

	while (args.size()) {
		if (args.front() == std::string("-p")) {
			if (args.size() == 1) {
				std::cerr << "Error: flag '-p' must be followed by the "
						<< "name of a parameter file."
						<< std::endl;
				print_usage_message ();
				exit (1);
			}
			args.pop_front ();
			const std::string parameter_file = args.front ();
			args.pop_front ();
			prm.read_input (parameter_file);

			prm.enter_subsection("Elasticity Parameters"); {
				lame_mu = prm.get_double("Lame's Parameter Mu");
				lame_lambda = prm.get_double("Lame's Parameter Lambda");
                critical_Gc = prm.get_double("Critical Fracture Energy");
			}
			prm.leave_subsection();

			prm.enter_subsection("DG Parameters"); {
				penalty_DG = prm.get_double("Face Term Penalty");
				sform = prm.get_integer("sform");
				penalty_BC = prm.get_double("Boundary Term Penalty");
			}
			prm.leave_subsection();

			prm.enter_subsection("Phase-Field Parameters"); {
				pf_kappa = prm.get_double("Constant K");
				pf_epsilon = prm.get_double("Constant E");
                penalty_AL = prm.get_double("Penalty Gamma");
                max_iters_AL = prm.get_integer("Maximum Penalty Iterations");
                lintype = prm.get("Linearization Scheme");
			}
			prm.leave_subsection();

			prm.enter_subsection("Runtime Parameters"); {
                opttype = prm.get("Optimization Type");
                testcase = prm.get("Test Case");
				timestep = prm.get_double("Timestep Size");
				max_time_steps = prm.get_integer("Maximum Timesteps");
                decompose_stress_matrix = prm.get_bool("Decompose Stress Matrix");
                decompose_stress_rhs = prm.get_bool("Decompose Stress RHS");
			}
			prm.leave_subsection();

			flag_success = true;
		}
		else {
			args.pop_front ();
		}
	}

	if (flag_success == false) {
		std::cerr << "Error: Improperly specified input file." << std::endl;
		print_usage_message ();
		exit (1);
	}
}
