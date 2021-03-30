// -*- tab-width: 4; indent-tabs-mode: nil; c-basic-offset: 2 -*-
// vi: set et ts=4 sw=2 sts=2:

#ifdef HAVE_CONFIG_H
# include "config.h"
#endif
#define _USE_MATH_DEFINES
#include<math.h>
#include<iostream>
#include<fstream>
#include<vector>
#include<map>
#include<string>
#include<stdlib.h>
#include<time.h>
#include<exception>
#include<chrono>

/*********************************************************************/
/* SETTINGS
 *
 * options:
 * grig backend -> USE_YSP, USE_UG
 * mpi -> PARALLEL, SEQUENTIAL
 * problem ids:
 * 		- PROBLEM_SAMPLE1D
 * 		- PROBLEM_BENCHMARK2D
 * 		- PROBLEM_KAUNER_VALLEY_SYNTHETIC
 *
 * default: PROBLEM_SAMPLE1D, USE_YASP, SEQUENTIAL
 */
#define DIMENSION 2


#define PARALLEL

#define PROBLEM_BENCHMARK2D
#define USE_YASP

/*********************************************************************/

#include"csv_reader.hh"
#include"../include_dune.hh"
#include"benchmark2D/include_problem.hh"

/*********************************************************************/

int main(int argc, char **argv) {
	try {
		// Maybe initialize MPI
		Dune::MPIHelper &helper = Dune::MPIHelper::instance(argc, argv);
		if (helper.rank() == 0) {
			std::cout << "This is SoilFreezing/permafrost project."
					<< std::endl;
		}

		if (Dune::MPIHelper::isFake) {
			std::cout << "This is a sequential program." << std::endl;
		} else {
			std::cout << "I am rank " << helper.rank() << " of "
					<< helper.size() << " processes!" << std::endl;
		}
		
		/**************************************************************************************************/
		// INPUTS
		if (argc != 4) {
			if (helper.rank() == 0) {
				// Expected input 
				std::cout
						<< "usage: ./permafrost <problem_name> <inputs_file.ini> <user_name>"
						<< '\n' << "user names: jonas, shubhangi, server, m2 ... or the path to the dune folder"
						<< '\n' << "benchmark2D"
						<< std::endl;
			}
		}

		//PROBLEM-NAME
		char PROBLEM[40];
		sscanf(argv[1], "%39s", PROBLEM);

		//USER-NAME
		char USER[40]; 
		sscanf(argv[3], "%39s", USER); 
		
		std::string DUNE_PATH(USER); 
		if (DUNE_PATH == "jonas"){
			DUNE_PATH =  "/home/jonas/dune_2_8/"; 
		}else if(DUNE_PATH == "shubhangi"){
			DUNE_PATH = "/home/sgupta/dune_2_8/";
		}else if (DUNE_PATH == "server"){
			DUNE_PATH = "/mnt/data01/beddrich/dune_2_8/"; 
		}else if (DUNE_PATH == "m2"){
			DUNE_PATH = "/home/jonas/DUNE/"; 
		}

		// PROBLEM PATH NAME
		std::string PROBLEM_PATH = DUNE_PATH + "SoilFreezing/src/TH2_Benchmark/";
		PROBLEM_PATH += PROBLEM;
		PROBLEM_PATH += "/";
		
		// INPUT PATH NAME
		std::string INPUT_PATH = PROBLEM_PATH + "inputs/";

		// OUTPUT PATH NAME
		std::string OUTPUT_PATH = DUNE_PATH + "SoilFreezing/src/outputs/TH2_Benchmark/";
		OUTPUT_PATH += PROBLEM;
		OUTPUT_PATH += "/";

		// INI-FILE FOR USER-DEFINED INPUTS
		char INI_FILE[40];
		sscanf(argv[2], "%39s", INI_FILE);
		std::string input_file = INPUT_PATH;
		input_file += INI_FILE;
		
		if (helper.rank() == 0) {
			 std::cout << "input file: " << input_file << std::endl;
			 std::cout << "problem: " << PROBLEM << std::endl;
		}
		
		// PARAMETER TREE
		Dune::ParameterTree ptree;
		Dune::ParameterTreeParser ptreeparser;
		ptreeparser.readINITree(input_file, ptree);
		ptreeparser.readOptions(argc, argv, ptree);
		
		/**************************************************************************************************/
		// !!! Fixing factors for units 
		CharacteristicValues < Dune::ParameterTree > characteristicValue(ptree);
		// MESH
		using GmshIndexMap = std::vector<int>;
		GmshIndexMap boundary_index_map;
		GmshIndexMap element_index_map;
		
#if DIMENSION==2
		const int dim = 2;
#elif DIMENSION==3
		const int dim = 3;
#else
		std::cout<< "Invalid DIMENSION " << DIMENSION << std::endl;
		exit(0);
#endif
		
#ifdef USE_YASP
		/*************
		 *  YASP
		 *************/
		double xc = characteristicValue.get_x();

		double Xmax = ptree.get("grid.yasp.LX", (double) 10.) / xc; /*ndim*/
		int nX = ptree.get("grid.yasp.NX", (int) 500);
		double Zmax = ptree.get("grid.yasp.LZ", (double) 10.) / xc; /*ndim*/
		int nZ = ptree.get("grid.yasp.NZ", (int) 500);
#if DIMENSION==3
		double Ymax = ptree.get("grid.yasp.LY", (double) 10.) / xc; /*ndim*/
		int nY = ptree.get("grid.yasp.NY", (int) 500);
#endif

		Dune::FieldVector<double, dim> L;
		std::array<int, dim> N;
		L[0] = Xmax;
		N[0] = nX;
#if DIMENSION==2
		L[1] = Zmax;
		N[1] = nZ;
#elif DIMENSION==3
		L[1] = Ymax;
		N[1] = nY;
		L[2] = Zmax;
		N[2] = nZ;
#endif
		
		std::bitset < dim > periodic(false);
		int overlap = 1;
		using Grid = Dune::YaspGrid<dim>;
		std::shared_ptr<Grid> grid = std::shared_ptr < Grid
				> (new Grid(L, N, periodic, overlap, helper.getCommunicator()));
		using GV = Grid::LeafGridView;
		GV gv = grid->leafGridView();
		grid->loadBalance();

#elif defined(USE_UG)

		/*************
		 *  UG
		 *************/

		using GridType = Dune::UGGrid<dim>;
		GridType grid_type;
		const std::string grid_file_name = ptree.get("grid.ug.name",(std::string)"sample");
		auto grid_file = PROBLEM_PATH + "grids/";
		grid_file += grid_file_name;
		Dune::GmshReader<GridType> gmshreader;
		std::shared_ptr<GridType> grid(gmshreader.read(grid_file,boundary_index_map, element_index_map,true,false));

		using GV = GridType::LeafGridView;
		GV gv = grid->leafGridView();
        grid->loadBalance();

#elif defined(USE_MESH_GENERATOR)
		/******************
		 * MESH GENERATOR 
		 *****************/
        using GridType = Dune::UGGrid<dim>;
        CharacteristicValues<Dune::ParameterTree> xc(ptree);		
		using MeshParams = MeshParameters<GmshIndexMap,Dune::ParameterTree,CharacteristicValues<Dune::ParameterTree>>;
		MeshParams meshparam(element_index_map,ptree,xc,PROBLEM_PATH);
		
#if DIMENSION==2
		ImplicitSurfaceMesh2DQuad<GridType, MeshParams> reactor(meshparam);
#elif DIMENSION==3
		ImplicitSurfaceMesh3DQuad<GridType, MeshParams> reactor(meshparam);
#endif
		reactor.grid().globalRefine(0);
		using GV = GridType::LeafGridView;
		GV gv = reactor.grid().leafGridView();
		reactor.grid().loadBalance();

#else
        std::cout<< "Incorrect grid manager. Please choose YASP or UG." << std::endl;
        exit(0);
#endif
       
        
        /**************************************************************************************************/
		// DRIVER
		driver(gv, ptree, PROBLEM_PATH, OUTPUT_PATH, element_index_map, boundary_index_map, helper);

		/**************************************************************************************************/

	} catch (Dune::Exception &e) {
		std::cerr << "Dune reported error: " << e << std::endl;
	} catch (...) {
		std::cerr << "Unknown exception thrown!" << std::endl;
	}
}