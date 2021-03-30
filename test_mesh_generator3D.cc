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

#include"../include_dune.hh"

#include"grids/implicit_surface_mesh_generator3D.hh"

class MeshParameters {
private:
	double eps = 1.0e-6;

	// dimensions of the grid 
	double Zmax;
	double Xmax;
	double Ymax;
	// cellsize 
	double dz;
	double dx;
	double dy;

public:

	MeshParameters() {
		Zmax = 10;
		Xmax = 10;
		Ymax = 10;
		dz = 0.5;
		dx = 0.5;
		dy = 0.5;
	}

	/*
	 * 3D -> X, Z, Y 
	 */

	const static int dimension = 3;

	double get_Xmax() const {
		return Xmax;
	}

	double get_Ymax() const {
		return Ymax;
	}

	double get_Zmax() const {
		return Zmax;
	}

	double get_dx() const {
		return dx;
	}

	double get_dz() const {
		return dz;
	}

	double get_dy() const {
		return dy;
	}

	bool isInsideImplicitSurface(
			Dune::FieldVector<double, dimension> vertexpos) const {
		// Failing code example 
//		if ((vertexpos[0] > 6  && vertexpos[2] > 6 && (vertexpos[1] < 4 || vertexpos[1] > 6))) {
//			return false;	

		// Successfull code example 
		if ((vertexpos[0] > 5 && vertexpos[2] > 5)) {
			return false;
		} else {
			return true;
		}
	}
};

int main(int argc, char **argv) {
	try {
		// Maybe initialize MPI
		Dune::MPIHelper &helper = Dune::MPIHelper::instance(argc, argv);
		std::cout << "Hello World! This is Mesh Generator Testing. Whuzzzup?"
				<< std::endl;
		if (Dune::MPIHelper::isFake)
			std::cout << "This is a sequential program." << std::endl;
		else
			std::cout << "I am rank " << helper.rank() << " of"
					<< helper.size() << " processes!" << std::endl;

		// MESH
		const int dim = 3;
		using GridType = Dune::UGGrid<dim>;
		MeshParameters meshparam;

		ImplicitSurfaceMesh3DQuad<GridType, MeshParameters> reactor(meshparam);

		reactor.grid().globalRefine(1);

		using GV = GridType::LeafGridView;
		GV gv = reactor.grid().leafGridView();
		reactor.grid().loadBalance();

		std::string mesh_name = "/home/sgupta/dune_2_8/SoilFreezing/src/permafrost/grids/";
		mesh_name += "reactor3D";
		Dune::VTKWriter<GV> vtkwriter(gv);
		vtkwriter.write(mesh_name, Dune::VTK::appendedraw);

		return 0;
	} catch (Dune::Exception &e) {
		std::cerr << "Dune reported error: " << e << std::endl;
	} catch (...) {
		std::cerr << "Unknown exception thrown!" << std::endl;
	}
}
