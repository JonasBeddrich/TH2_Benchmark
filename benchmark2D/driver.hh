template<typename GV, typename PTree, typename GmshIndexMap>
void driver(const GV &gv, 				// GridView
			const PTree &ptree, 		// Input Parameters
			std::string problem_path,
			std::string output_path,
			GmshIndexMap boundary_index_map,
			GmshIndexMap element_index_map,
			Dune::MPIHelper &helper) {

	//	CHOOSE DOMAIN AND RANGE FIELD TYPE
	using Coord = typename GV::Grid::ctype;
	const int dim = GV::dimension;
	double time = 0.0;
	double dt = 0.0;

	//	MATERIAL PROPERTIES, NUMERICAL AND TEST PARAMETERS, CONSTITUTIVE RELATIONSHIPS
	typedef Properties<GV, PTree, GmshIndexMap> Properties;
	Properties property(gv, ptree, element_index_map, &time, &dt);

	std::string pathName = output_path;
	auto pathExt = ptree.get("output.path_name", (std::string) "test0");
	pathName += pathExt;
	pathName += "/";
	
	auto hh = ptree.get("benchmark.hydraulic_head", (std::string) "0.09"); 
	hh.erase(std::remove(hh.begin(), hh.end(), '.'), hh.end());
	auto fileName = ptree.get("output.file_name", (std::string) "test") + "_" + hh; 

	/*Non-dimensionalize time prams*/
	//PROBLEM time and dt
	auto Xc_t = property.characteristicValue.get_t();
	dt = ptree.get("time.dt_initial", (double) 0.001);
	dt *= 1. / Xc_t; /*ndim*/
	double t_END = ptree.get("time.time_end", (double) 86400.);
	t_END *= 1. / Xc_t; /*ndim*/
	// output time interval
	double t_OP = ptree.get("output.time_interval", (double) 1440.);
	t_OP *= 1. / Xc_t; /*ndim*/

	//adaptive time control
	bool adaptive_time_control = ptree.get("adaptive_time_control.flag",
			(bool) true);
	double dt_min = ptree.get("adaptive_time_control.dt_min", (double) 1.e-6);
	dt_min *= 1. / Xc_t; /*ndim*/
	double dt_max = ptree.get("adaptive_time_control.dt_max", (double) 10.);
	dt_max *= 1. / Xc_t; /*ndim*/
	int maxAllowableIterations = ptree.get(
			"adaptive_time_control.max_newton_steps", (int) 6);
	int minAllowableIterations = ptree.get(
			"adaptive_time_control.min_newton_steps", (int) 4);

	double dtstart = dt;
	double time_op = time;
	double clock_time_elapsed = 0.;

	/************************************************************************************************/
	// MAIN
	/************************************************************************************************/
	//	COMPOSITE GFS FOR PRIMARY VARIABLES 
#ifdef PARALLEL
		using CON0 = Dune::PDELab::P0ParallelConstraints;
#else
	using CON0 = Dune::PDELab::NoConstraints;
#endif

	using VBE0 = Dune::PDELab::ISTL::VectorBackend<>;
	using FEM0 = Dune::PDELab::QkDGLocalFiniteElementMap<Coord,double,0,dim,Dune::PDELab::QkDGBasisPolynomial::lagrange>;
	FEM0 fem0;
	using GFS0 = Dune::PDELab::GridFunctionSpace<GV, FEM0, CON0, VBE0>;
	GFS0 gfs0(gv, fem0);
	// gfs for composite system: Pa , Sw , T , por , Ci

#ifdef PARALLEL
		using VBE = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::fixed>;
		using GFS = Dune::PDELab::PowerGridFunctionSpace<GFS0,Indices::numOfPVs,VBE,Dune::PDELab::EntityBlockedOrderingTag>;
		GFS gfs(gfs0);
#else
	using VBE = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none>;
	using GFS = Dune::PDELab::PowerGridFunctionSpace<GFS0,Indices::numOfPVs,VBE,Dune::PDELab::LexicographicOrderingTag>;
	GFS gfs(gfs0);
#endif

	using CC = typename GFS::template ConstraintsContainer<double>::Type;
	CC cc;
	cc.clear();

	//	SUB-SPACES FOR ACCESSING PRIMARY VARIABLES
	using PathP = Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::PVId_P>>;
	using SUBGFS_P = Dune::PDELab::GridFunctionSubSpace<GFS,PathP>;
	SUBGFS_P subgfs_P(gfs);
	using PathT = Dune::TypeTree::HybridTreePath<Dune::index_constant<Indices::PVId_T>>;
	using SUBGFS_T = Dune::PDELab::GridFunctionSubSpace<GFS,PathT>;
	SUBGFS_T subgfs_T(gfs);

	//	MAKE VECTOR CONTAINER FOR THE SOLUTION
	using U = Dune::PDELab::Backend::Vector<GFS,double>;
	U u_old(gfs);
	U u_new(gfs);

	//EVALUTION FUNCTIONS FOR TIME DERIVATIVES IN LOCAL OPERATOR 
	Dune::PDELab::Evaluation<SUBGFS_P, U> evaluation_P_old(subgfs_P, &u_old); 
	Dune::PDELab::Evaluation<SUBGFS_T, U> evaluation_T_old(subgfs_T, &u_old); 

	//	INITIAL CONDITIONS
	using ICV_P = P_Initial<GV,Properties,GmshIndexMap>;
	ICV_P P_initial(gv, property, element_index_map);
	using ICV_T = T_Initial<GV,Properties,GmshIndexMap>;
	ICV_T T_initial(gv, property, element_index_map);
	
	using InitialValues = Dune::PDELab::CompositeGridFunction< ICV_P, ICV_T>;
	InitialValues icv(P_initial, T_initial);

	Dune::PDELab::interpolate(icv, gfs, u_old);
	u_new = u_old;

	//	BOUNDARY CONDITIONS
	using BoundaryConditions = ProblemBoundaryConditions<GV,Properties,GmshIndexMap>;
	BoundaryConditions bc(gv, property, boundary_index_map);

	//	MAKE INSTATIONARY GRID OPERATOR SPACE
	// spatial part 
	using LOP = LocalOperator< GV, Properties, BoundaryConditions, 
	Dune::PDELab::Evaluation<SUBGFS_P ,U>,
	Dune::PDELab::Evaluation<SUBGFS_T ,U>>;
	LOP lop(gv, property, bc, &evaluation_P_old, &evaluation_T_old, &time, &dt);

	// temporal part
	using TLOP = TimeOperator< GV, Properties >;
	TLOP tlop(gv, property, &time, &dt);

	using MBE = Dune::PDELab::ISTL::BCRSMatrixBackend<>;
	MBE mbe(5);

	using GOLOP = Dune::PDELab::GridOperator< GFS, GFS, LOP, MBE, double, double, double, CC, CC>;
	GOLOP golop(gfs, cc, gfs, cc, lop, mbe);

	typename GOLOP::Traits::Jacobian jac(golop);
	if (helper.rank() == 0) {
		std::cout << jac.patternStatistics() << std::endl;
	}

	using GOTLOP = Dune::PDELab::GridOperator< GFS, GFS, TLOP, MBE, double, double, double, CC, CC >;
	GOTLOP gotlop(gfs, cc, gfs, cc, tlop, mbe);

	using IGO = Dune::PDELab::OneStepGridOperator< GOLOP, GOTLOP>;
	IGO igo(golop, gotlop);

	// SELECT A LINEAR SOLVER BACKEND
#ifdef PARALLEL
		using LS = Dune::PDELab::ISTLBackend_BCGS_AMG_SSOR<IGO>; //works
		LS ls(gfs,100,1,false,true);
#else
	using LS = Dune::PDELab::ISTLBackend_SEQ_SuperLU;
	LS ls(0);
#endif

	//	SELECT SOLVER FOR NON-LINEAR PROBLEM
	using PDESOLVER = Dune::PDELab::Newton< IGO, LS, U >;
	PDESOLVER pdesolver(igo, ls);
	// 	select control parameters for non-linear PDE-solver
	// pdesolver.setLineSearchStrategy(PDESOLVER::Strategy::noLineSearch);

	pdesolver.setReassembleThreshold(0.0);
	pdesolver.setVerbosityLevel(2);
	pdesolver.setReduction(1e-6);
	pdesolver.setMinLinearReduction(1e-4);
	pdesolver.setMaxIterations(ptree.get("newton.max_iterations", (int) 15));
	pdesolver.setForceIteration(true);
	pdesolver.setAbsoluteLimit(ptree.get("newton.abs_error", (double) 1.e-3));

	//	SELECT TIME-STEPPER
	Dune::PDELab::ImplicitEulerParameter<double> method;
	Dune::PDELab::OneStepMethod<double, IGO, PDESOLVER, U, U> osm(method, igo,
			pdesolver);
	osm.setVerbosityLevel(2);

	/************************************************************************************************/
	//  POST-PROCESS:
	/************************************************************************************************/

	// using VBE_PP = Dune::PDELab::ISTL::VectorBackend<Dune::PDELab::ISTL::Blocking::none >;
	// using GFS_PP = Dune::PDELab::PowerGridFunctionSpace< GFS0,Indices::numOfSVs,VBE_PP, Dune::PDELab::LexicographicOrderingTag >;
	// GFS_PP gfs_pp(gfs0);
	// using CC_PP = typename GFS_PP::template ConstraintsContainer<double>::Type;
	// CC_PP cc_pp;
	// cc_pp.clear();

	// using U_PP = Dune::PDELab::Backend::Vector<GFS_PP,double>;
	// U_PP u_pp(gfs_pp, 0.0);

	// PostProcess<GV, Properties, Dune::PDELab::Evaluation<SUBGFS_Pa, U>,
	// 		Dune::PDELab::Evaluation<SUBGFS_Sw, U>,
	// 		Dune::PDELab::Evaluation<SUBGFS_T, U>,
	// 		Dune::PDELab::Evaluation<SUBGFS_por, U>,
	// 		Dune::PDELab::Evaluation<SUBGFS_Ci, U>, GFS_PP, U_PP> postprocess(
	// 		gv, property, &evaluation_Pa, &evaluation_Sw, &evaluation_T,
	// 		&evaluation_por, &evaluation_Ci, gfs_pp, &u_pp);
	// postprocess.evaluate();

	/************************************************************************************************/
	// VTK
	/************************************************************************************************/
	
	// Add names to the components for VTK output
	using namespace Dune::TypeTree::Indices;
	for (int i = 0; i < property.index.numOfPVs; i++) {
		gfs.child(i).name(property.index.getPVName(i));
		std::cout << property.index.getPVName(i) << std::endl; 
	}

	int subsampling = 1;
	using VTKWRITER = Dune::SubsamplingVTKWriter<GV>;
	VTKWRITER vtkwriter(gv, Dune::refinementIntervals(subsampling));
	using VTKSEQUENCEWRITER = Dune::VTKSequenceWriter<GV>;
	VTKSEQUENCEWRITER vtkSequenceWriter(std::make_shared < VTKWRITER > (vtkwriter), fileName, pathName, "");
	// TODO: check whether I should use u_new or u_old here ... 
	Dune::PDELab::addSolutionToVTKWriter(vtkSequenceWriter, gfs, u_old);
	vtkSequenceWriter.write(time, Dune::VTK::appendedraw);

	/***********************************************/
	//	BEGIN TIME LOOP
	/***********************************************/

	int opcount = 1;
	double timecount = time;
	double dtLast = dtstart;
	int dtFlag = 0;
	bool exceptionCaught = false;
	bool cutOutputFlag = false; 
	int newton_iterations = 0;

	while (time < t_END - 1e-8 / Xc_t) {

		if (exceptionCaught == false) {
			dt = std::max(dt, dt_min);
		}

		if (helper.rank() == 0) {
			std::cout << "_____________________________________________________"
					<< std::endl;
			std::cout << " current opcount = " << opcount - 1 << std::endl;
		}

		clock_t start = clock();

		try {
			if (helper.rank() == 0) {
				std::cout << "****************************" << std::endl;
				std::cout << "  CALLING osm.apply() !" << std::endl;
				std::cout << "****************************" << std::endl;
			}

			osm.apply(time, dt, u_old, u_new);

			newton_iterations = osm.getPDESolver().result().iterations;

			exceptionCaught = false;

		} catch (Dune::Exception &e) {
			exceptionCaught = true;
			if (dt * Xc_t > 1e-8) {

				if (helper.rank() == 0) {
					std::cout << "Catched Error, Dune reported error: " << e
							<< std::endl;
				}

				u_new = u_old;

				newton_iterations = 0;

				dt *= 0.5;

				// if the algorithm doesn't converge the time step gets small 
				if (dt < dt_min / 100){
					exit(0); 
				}

				continue;

			} else {
				if (helper.rank() == 0) {
					std::cout << "ABORTING, due to DUNE error: " << e
							<< std::endl;
				}
				exit(0);				
			}
		}

		clock_t end = clock();
		double clock_time_this_step = (double) (end - start) / CLOCKS_PER_SEC;
		clock_time_elapsed += clock_time_this_step;

		if (helper.rank() == 0) {
			std::cout << "DONE" << std::endl;
			std::cout << "_____________________________________________________"
					<< std::endl;
		}

		/*********************************************************************************************
		 * OUTPUT
		 *********************************************************************************************/
		/* At each time step: **Statistics**
		 * t_new,
		 * dt,
		 * fixed point iterations,
		 * newton iterations per fixed point iteration,
		 * total newton terations
		 */

		if (helper.rank() == 0) {
			std::string statistics_file = pathName;
			statistics_file += fileName;
			statistics_file += "_statistics";
			statistics_file += ".txt";
			property.ReportStatistics(statistics_file, time * Xc_t, dt * Xc_t,
					newton_iterations, clock_time_elapsed);
		}

		/* At t_OP
		 *
		 */
		if ((time + dt > t_OP * opcount - dt_min)
				and (time + dt <= t_OP * opcount)) {
			// POST PROCESS FOR NEW OUTPUT
			// postprocess.evaluate();

			// WRITE OUTPUT
			vtkSequenceWriter.write(time, Dune::VTK::appendedraw);

			if (helper.rank() == 0) {
				std::cout
						<< " ******************************************************************* "
						<< std::endl;
				std::cout << " OUTPUT WRITTEN " << opcount << std::endl;
				std::cout
						<< " ******************************************************************* "
						<< std::endl;
				std::cout << std::flush;
			}

			timecount = time;
			opcount = opcount + 1;
		}

		//	PREPARE FOR NEXT TIME INTEGRATION
		//	1. ASSIGN THE 'NEXT' VALUE TO 'OLD' VARIABLE
		u_old = u_new;
		//	2. ADVANCE TIME:
		time += dt;

		if (helper.rank() == 0) {
			std::cout << " " << std::endl;
			std::cout << " time = " << time * Xc_t;
			std::cout << std::flush;
		}

		if (adaptive_time_control) {
			if (newton_iterations > maxAllowableIterations) {
				dt = std::max(dt * 0.9, dt_min);
			} else if (newton_iterations <= minAllowableIterations) {
				dt = std::min(dt * 1.1, dt_max);
			}
		} else {
			dt = dtstart;
		}

		if (helper.rank() == 0) {
			std::cout << " , time+dt = " << (time + dt) * Xc_t << " , opTime = "
					<< t_OP * opcount * Xc_t;
			std::cout << std::flush;
		}

		if(cutOutputFlag){
			dt = dtLast; 
		 	cutOutputFlag = false; 
		}

		if (time + dt > t_OP * opcount) {
			dtLast = dt;
			cutOutputFlag = true; 
			dt = t_OP * opcount - time;

			if (helper.rank() == 0) {
				std::cout << " , because timeNext > opNext , dt set to : "
						<< dt * Xc_t << std::endl;
				std::cout << std::flush;
			}
			dtFlag = 0;
		}
		dtFlag += 1;

		if (helper.rank() == 0) {
			std::cout << " , dt  : " << dt * Xc_t << std::endl;
			std::cout << " " << std::endl;
			std::cout << " READY FOR NEXT ITERATION. " << std::endl;
			std::cout << std::flush;
		}
	}
}