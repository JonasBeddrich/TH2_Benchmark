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
#include<stdio.h>

#include"csv_reader.hh"
#include"../include_dune.hh"
#include"grids/implicit_surface_mesh_generator.hh"

class MeshParametersKAUNER{
private:
	double eps = 1.0e-6;

	// dimensions of the grid
	double z_max;
	double z_min;
	double Zmax;
	double Xmax;

	// cell-size
	double dz;
	double dx;

	double point_distance;
	double n_points;
	double minimal_depth;

	double buffer_height;
	double n;
	double alpha;

	std::vector<double> x_Data;
	std::vector<double> z_Data;
	std::string filename;
	csv_reader csv_r;

	// Just testing ... and yes I know that's not the right place to do it ... whatever
	std::vector<std::vector<double>> my_matrix;

public:

	MeshParametersKAUNER(){

		// parameters ...
		point_distance = 25.;
		minimal_depth = 30;

		// these are not chosen randomly but to recreate a buffer_height below the uniform mesh of 70
		buffer_height = 50;
		n = 30;
		alpha = 1.05;

		filename = "/home/sgupta/dune_2_8/SoilFreezing/permafrost/kauner_valley_synthetic/data/Elevation_cross_profile.csv";

		// create x and z vector for the interpolation
		z_Data = csv_r.read_1D_axis(filename);
		n_points = z_Data.size();
		x_Data = std::vector<double>(n_points);

		for (int i = 0; i < x_Data.size(); i++) {
			x_Data[i] = i * point_distance;
		}

		// compute grid size
		z_max = *std::max_element(z_Data.begin(), z_Data.end());
		z_min = *std::min_element(z_Data.begin(), z_Data.end());
		Zmax = (z_max - z_min) + minimal_depth;
		Xmax = n_points * point_distance;
		dz = 1;
		dx = 3;
	}

	/* thanks :D
	 * 2D -> X and Z
	 */

	double get_n() const {
		return n;
	}

	double get_alpha() const {
		return alpha;
	}

	double get_Xmax() const {
		return Xmax;
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

	bool isInsideImplicitSurface(Dune::FieldVector<double,2> vertexpos,
			double buffer_height) const {
		linear_interpolator l_int(x_Data, z_Data);
		if (vertexpos[1] - buffer_height
				< l_int.interpolate(vertexpos[0]) - z_min + minimal_depth) {
			return true;
		} else
			return false;
	}
};

int main(int argc, char **argv) {
	try {
		// Maybe initialize MPI
		Dune::MPIHelper &helper = Dune::MPIHelper::instance(argc, argv);
		std::cout << "Hello World! This is Mesh Generator Testing?" << std::endl;

		if (Dune::MPIHelper::isFake)
			std::cout << "This is a sequential program." << std::endl;
		else
			std::cout << "I am rank " << helper.rank() << " of " << helper.size() << " processes!" << std::endl;

		// MESH
		const int dim = 2;
		using GridType = Dune::UGGrid<dim>;
		using MeshParameters= MeshParametersKAUNER;
		MeshParameters meshparam;

		ImplicitSurfaceMesh2DQuad<GridType, MeshParameters> reactor(meshparam);
		reactor.grid().globalRefine(0);

		using GV = GridType::LeafGridView;
		GV gv = reactor.grid().leafGridView();
		reactor.grid().loadBalance();

		std::string mesh_name = "/home/sgupta/dune_2_8/SoilFreezing/src/permafrost/grids/";
		mesh_name += "Kauner_valley_test";
		Dune::VTKWriter<GV> vtkwriter(gv);
		vtkwriter.write(mesh_name, Dune::VTK::appendedraw);

		return 0;

	} catch (Dune::Exception &e) {
		std::cerr << "Dune reported error: " << e << std::endl;
	} catch (...) {
		std::cerr << "Unknown exception thrown!" << std::endl;
	}
}
