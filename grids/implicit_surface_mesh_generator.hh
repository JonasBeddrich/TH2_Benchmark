/*
 * implicit_surface_mesh_generator.hh
 *
 *  Created on: Apr 27, 2020
 *      Author: sgupta
 */

#ifndef GRIDS_IMPLICIT_SURFACE_MESH_GENERATOR_HH_
#define GRIDS_IMPLICIT_SURFACE_MESH_GENERATOR_HH_

template<typename T, typename MeshParam>
class ImplicitSurfaceMesh2DQuad {
private:
	typedef T GridType;
	const MeshParam &param;
	static const int dim = GridType::dimension;

	// the constructed grid object
	std::shared_ptr<T> grid_;

public:
	// constructor throwing exception
	ImplicitSurfaceMesh2DQuad(const MeshParam &param_) :
			param(param_) {
		double Xmax = param.get_Xmax();
		double dx = param.get_dx();
		int num_col = (int) (Xmax / dx + 1);

		double Zmax = param.get_Zmax();
		double dz = param.get_dz();
		int num_row = (int) (Zmax / dz + 1);

		double alpha = param.get_alpha();
		int n = param.get_n();
		double buffer_height = (pow(alpha, n + 1) - 1) / (alpha - 1) * dz;
		std::cout << "buffer height: " << buffer_height << '\t'
				  << "n: " << n << '\t'
				  << "aspect_ratio: " << alpha << '\t'
				  << std::endl;

		std::cout << "Xmax: " << Xmax << '\t'
				  << "Zmax: " << Zmax << '\t'
				  << "dx: " << dx << '\t'
				  << "dz: " << dz << '\t'
				  << "num_col: " << num_col << '\t'
				  << "num_row: " << num_row << std::endl;

		Dune::GridFactory<GridType> factory;
		Dune::FieldVector<double, dim> v(0.0);
		std::vector<int> vnum(num_col, 0);

		std::vector<std::vector<double>> vid(num_col);
		for (int nc = 0; nc < vid.size(); nc++) {
			vid[nc] = std::vector<double>(num_row + n + 1, 0.);
		}

		// fill vnum and vid with values 
		int tmp = 0;
		for (int nc = 0; nc < num_col; nc++) {
			v[0] = nc * dx;
			vnum[nc] = 0;

			double height = 0;
			// i = 0, 1, 2 for n = 2 
			for(int i =0; i < n+1; i++ ){
				if (i == 0){
					v[1] = 0; 
				}else{
					// n+1-i = n for i = 1 and ... = n-1 for i = 2
					height += pow(alpha, n+1-i); 
					v[1] = height * dz;	
				}
				vnum[nc] += 1; 
				factory.insertVertex(v); 
				vid[nc][i] = tmp; 
				tmp += 1; 
			}
			
			// nr = 3, 4, 5, ... for n = 2 
			for (int nr = n + 1; nr < num_row + n + 1; nr++) {
				v[1] = (nr - n + height) * dz;
				if (param.isInsideImplicitSurface(v, height*dz)) {
					vnum[nc] += 1;
					factory.insertVertex(v);
					vid[nc][nr] = tmp;
					tmp += 1;
				}
			}
		}

		std::vector<unsigned int> cornerIDs(4);
		std::vector<unsigned int> boundaryIDs(2);

		for (int nc = 0; nc < num_col - 1; nc++) {
			for (int nr = 0; nr < vnum[nc] - 1; nr++) {
				if (nr + 1 < vnum[nc + 1]) {
					// 0: bottom left
					// 1: bottom right 
					// 2: upper left 
					// 3: upper right 
					cornerIDs[0] = vid[nc][nr];
					cornerIDs[1] = vid[nc + 1][nr];
					cornerIDs[2] = vid[nc][nr + 1];
					cornerIDs[3] = vid[nc + 1][nr + 1];
					const Dune::GeometryType type(
							Dune::GeometryTypes::quadrilateral);
					factory.insertElement(type, cornerIDs);
					// add top boundary 
					if (nr + 2 == vnum[nc + 1]) {
						//corners 2,3 
						boundaryIDs[0] = vid[nc][nr + 1];
						boundaryIDs[1] = vid[nc + 1][nr + 1];
						factory.insertBoundarySegment(boundaryIDs);
					}
					// add right boundaries inside the domain 
					if (nc < num_col - 2 && nr + 1 > vnum[nc + 2]) {
						//corners 2,3 
						boundaryIDs[0] = vid[nc + 1][nr];
						boundaryIDs[1] = vid[nc + 1][nr + 1];
						factory.insertBoundarySegment(boundaryIDs);
					}
					// add left boundaries inside the domain  
					// also check that we're not all the way to the left 
					if (nc > 0 && nr + 1 > vnum[nc - 1]) {
						//corners 2,3 
						boundaryIDs[0] = vid[nc][nr];
						boundaryIDs[1] = vid[nc][nr + 1];
						factory.insertBoundarySegment(boundaryIDs);
					}
				}
			}
		}
		// Finish
		grid_ = factory.createGrid();
	}

	T& grid() {
		return *grid_;
	}

};

#endif /* GRIDS_IMPLICIT_SURFACE_MESH_GENERATOR_HH_ */
