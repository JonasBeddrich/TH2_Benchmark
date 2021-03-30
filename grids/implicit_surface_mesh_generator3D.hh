/*
 * implicit_surface_mesh_generator.hh
 *
 *  Created on: Apr 27, 2020
 *      Author: sgupta
 */

#ifndef GRIDS_IMPLICIT_SURFACE_MESH_GENERATOR_HH_
#define GRIDS_IMPLICIT_SURFACE_MESH_GENERATOR_HH_

template<typename T, typename MeshParam>
class ImplicitSurfaceMesh3DQuad {
private:
	typedef T GridType;
	const MeshParam &param;
	static const int dim = GridType::dimension;

	// the constructed grid object
	std::shared_ptr<T> grid_;

public:
	// constructor throwing exception
	ImplicitSurfaceMesh3DQuad(const MeshParam &param_) :
			param(param_) {
		double Xmax = param.get_Xmax();
		double dx = param.get_dx();
		int num_col = (int) (Xmax / dx + 1);

		double Zmax = param.get_Zmax();
		double dz = param.get_dz();
		int num_row = (int) (Zmax / dz + 1);

		double Ymax = param.get_Ymax();
		double dy = param.get_dy();
		int num_lay = (int) (Ymax / dy + 1);

		std::cout << Xmax << '\t' << Zmax << '\t' << Ymax << '\t' << dx << '\t'
				<< dz << '\t' << dy << '\t' << num_col << '\t' << num_row
				<< '\t' << num_lay << '\t' << std::endl;

		// start grid creation
		// elements will be inserted by an iteration through both dimensions followed up by create grid 
		Dune::GridFactory<GridType> factory;

		// Create position triple to iterate through everything 
		// Create vnum matrix denoting how many entries are relevant in each column 
		// Create vid tensor enumerating all the elements 

		// position pair 
		// to iterate through the grid 
		Dune::FieldVector<double, dim> v(0.0);

		// matrix
		// number of inside elements per column
		std::vector<std::vector<int>> vnum(num_col);
		for (int nc = 0; nc < vnum.size(); nc++)
			vnum[nc] = std::vector<int>(num_lay, 0);

		// tensor 
		// 0 for outside 
		// 1,2,3 ... inside elements get enumerated 
		std::vector < std::vector<std::vector<double>> > vid(num_col);
		for (int nc = 0; nc < vid.size(); nc++) {
			vid[nc] = std::vector<std::vector<double>>(num_lay);
			for (int nl = 0; nl < vid[nc].size(); nl++) {
				vid[nc][nl] = std::vector<double>(num_row, 0.);
			}
		}

		// fill these two 
		int tmp = 0;
		for (int nc = 0; nc < num_col; nc++) {
			v[0] = nc * dx;
			for (int nl = 0; nl < num_lay; nl++) {
				v[1] = nl * dy;
				vnum[nc][nl] = 0;
				for (int nr = 0; nr < num_row; nr++) {
					v[2] = nr * dz;
					if (param.isInsideImplicitSurface(v)) {
						vnum[nc][nl] += 1;
						factory.insertVertex(v);
						vid[nc][nl][nr] = tmp;
						tmp += 1;
					}
				}
			}
		}

		// Insert elements 
		// Iteration runs over 
		// ... all columns and 
		// ... v_num[current column] rows i.e. number of rows underneath the surface 
		std::vector<unsigned int> cornerIDs(8);
		std::vector<unsigned int> boundaryIDs(4);
		for (int nc = 0; nc < num_col - 1; nc++) {
			for (int nl = 0; nl < num_lay - 1; nl++) {
				for (int nr = 0; nr < vnum[nc][nl] - 1; nr++) {
					// here I'm not entirely sure about what the indices have to be ... 
					if ((nr + 1 < vnum[nc + 1][nl + 1])
							&& (nr + 1 < vnum[nc][nl + 1])
							&& (nr + 1 < vnum[nc + 1][nl])) {
						// Insert element
						// this enumeration creates are "normal" grid of cubes 
						cornerIDs[0] = vid[nc][nl][nr];
						cornerIDs[1] = vid[nc + 1][nl][nr];
						cornerIDs[2] = vid[nc][nl][nr + 1];
						cornerIDs[3] = vid[nc + 1][nl][nr + 1];
						cornerIDs[4] = vid[nc][nl + 1][nr];
						cornerIDs[5] = vid[nc + 1][nl + 1][nr];
						cornerIDs[6] = vid[nc][nl + 1][nr + 1];
						cornerIDs[7] = vid[nc + 1][nl + 1][nr + 1];
						const Dune::GeometryType type(
								Dune::GeometryTypes::hexahedron);
						factory.insertElement(type, cornerIDs);

						// Addition boundaries inside the domain 
						// top boundary - corners 2,3,6,7 
						if ((nr + 2 == vnum[nc + 1][nl + 1])
								|| (nr + 2 == vnum[nc][nl + 1])
								|| (nr + 2 == vnum[nc + 1][nl])) {
							// corners 2,3,6,7
							boundaryIDs[0] = vid[nc][nl][nr + 1];
							boundaryIDs[1] = vid[nc + 1][nl][nr + 1];
							boundaryIDs[2] = vid[nc][nl + 1][nr + 1];
							boundaryIDs[3] = vid[nc + 1][nl + 1][nr + 1];
							factory.insertBoundarySegment(boundaryIDs);
						}
						// front boundary - corners 4,5,6,7 
						if (nl < num_lay - 2
								&& (nr + 1 > vnum[nc][nl + 2]
										|| nr + 1 > vnum[nc + 1][nl + 2]
										|| nr + 1 > vnum[nc + 1][nl + 1])) {
							boundaryIDs[0] = vid[nc][nl + 1][nr];
							boundaryIDs[1] = vid[nc + 1][nl + 1][nr];
							boundaryIDs[2] = vid[nc][nl + 1][nr + 1];
							boundaryIDs[3] = vid[nc + 1][nl + 1][nr + 1];
							factory.insertBoundarySegment(boundaryIDs);

						}
						// back boundary - corners 0,1,2,3
						if (nl > 0
								&& (nr + 1 > vnum[nc][nl]
										|| nr + 1 > vnum[nc + 1][nl]
										|| nr + 1 > vnum[nc + 1][nl - 1])) {
							boundaryIDs[0] = vid[nc][nl][nr];
							boundaryIDs[1] = vid[nc + 1][nl][nr];
							boundaryIDs[2] = vid[nc][nl][nr + 1];
							boundaryIDs[3] = vid[nc + 1][nl][nr + 1];
							factory.insertBoundarySegment(boundaryIDs);
						}
						// right boundary - corners 1,3,5,7 
						if (nc < num_col - 2
								&& (nr + 1 > vnum[nc + 1][nl + 1]
										|| nr + 1 > vnum[nc + 2][nl]
										|| nr + 1 > vnum[nc + 2][nl + 1])) {
							boundaryIDs[0] = vid[nc + 1][nl][nr];
							boundaryIDs[1] = vid[nc + 1][nl][nr + 1];
							boundaryIDs[2] = vid[nc + 1][nl + 1][nr];
							boundaryIDs[3] = vid[nc + 1][nl + 1][nr + 1];
							factory.insertBoundarySegment(boundaryIDs);
						}
						// left boundary - corners 0,2,4,6 
						if (nc > 0
								&& (nr + 1 > vnum[nc][nl]
										|| nr + 1 > vnum[nc][nl + 1]
										|| nr + 1 > vnum[nc - 1][nl + 1])) {
							boundaryIDs[0] = vid[nc][nl][nr];
							boundaryIDs[1] = vid[nc][nl][nr + 1];
							boundaryIDs[2] = vid[nc][nl + 1][nr];
							boundaryIDs[3] = vid[nc][nl + 1][nr + 1];
							factory.insertBoundarySegment(boundaryIDs);
						}
					}
				}
			}
		}
		std::cout << "I also get outside of the loop" << std::endl;
		// Finish
		grid_ = factory.createGrid();
	}
	T& grid() {
		return *grid_;
	}
};

#endif /* GRIDS_IMPLICIT_SURFACE_MESH_GENERATOR_HH_ */
