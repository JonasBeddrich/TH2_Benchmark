template<typename GV, typename Properties, typename PGMap>
class ProblemBoundaryConditions {
private:
	const GV &gv;
	const Properties &property;
	const PGMap &pgmap;
	ProblemInitialConditions<GV, Properties, PGMap> icvalue;
	const static int dim = GV::dimension;

public:

	// ! construct from gridview
	ProblemBoundaryConditions(const GV &gv_, const Properties &property_,
			const PGMap &pgmap_) :
			gv(gv_), property(property_), pgmap(pgmap_), icvalue(gv_, property_,
					pgmap_) {
	}

	/* boundary types */
	template<typename I> std::vector<int> type(I &intersection,
			const Dune::FieldVector<double, dim - 1> &xlocal, double time/*s*/,
			double dt/*s*/) const {

		auto xglobal = intersection.geometry().global(xlocal);
		std::vector<int> bct(Indices::numOfBCs, 0);

		if (property.mesh.isTopBoundary(xglobal)) {
			// zero conductive flux and zero flux 
			bct[Indices::BC_flux] = 1;
			bct[Indices::BC_heat] = 1;
		} else if (property.mesh.isBottomBoundary(xglobal)) {
			// zero conductive flux and zero flux 
			bct[Indices::BC_flux] = 1;
			bct[Indices::BC_heat] = 1;
		} else if (property.mesh.isLeftBoundary(xglobal)) {
			// T_in and imposed Head H0 + dH  
			bct[Indices::BC_flux] = 0;
			bct[Indices::BC_heat] = 0;
		} else if (property.mesh.isRightBoundary(xglobal)) {
			// zero conductive flux and imposed Head H0 
			bct[Indices::BC_flux] = 0;
			bct[Indices::BC_heat] = 1;
		}
		return bct;
	}

	/* boundary values */
	template<typename I> std::vector<double> value(I &intersection,
			const Dune::FieldVector<double, dim - 1> &xlocal,
			double time/*s*/,
			double dt/*s*/) const {

		auto xglobal = intersection.geometry().global(xlocal);
		std::vector<double> bcv(Indices::numOfBCs, 0.);

		if (property.mesh.isTopBoundary(xglobal)) {
			bcv[Indices::BC_flux] = 0.;
			bcv[Indices::BC_heat] = 0.;
		} else if (property.mesh.isBottomBoundary(xglobal)) {
			bcv[Indices::BC_flux] = 0.;
			bcv[Indices::BC_heat] = 0.;
		} else if (property.mesh.isLeftBoundary(xglobal)) {
			bcv[Indices::BC_flux] = property.parameter.get_P_left();
			bcv[Indices::BC_heat] = property.parameter.get_T0(); 
		} else if (property.mesh.isRightBoundary(xglobal)) {
			bcv[Indices::BC_flux] = property.parameter.get_P_right();
			bcv[Indices::BC_heat] = 0.;
		}
		return bcv;
	}

	// ! get a reference to the gridview
	inline const GV& getGridView() {
		return gv;
	}
};



