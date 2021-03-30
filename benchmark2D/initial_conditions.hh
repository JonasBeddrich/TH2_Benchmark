template<typename GV, typename Properties, typename PGMap>
class ProblemInitialConditions {
private:
	const GV &gv;
	const Properties &property;
	const PGMap &pgmap;
	const static int dim = GV::dimension;
	constexpr static double eps = 1.e-6;

public:

	//! construct from grid view
	ProblemInitialConditions(const GV &gv_, const Properties &property_,
			const PGMap &pgmap_) :
			gv(gv_), property(property_), pgmap(pgmap_) {
	}

	/* Initial Conditions */
	std::vector<double> evaluate(
			const typename GV::Traits::template Codim<0>::Entity &element,
			const Dune::FieldVector<double, dim> &xlocal) const {

		auto xglobal = element.geometry().global(xlocal);
		std::vector<double> icvalue(Indices::numOfPVs, 0.);

		// ice cube 
		if (property.mesh.isIceCube(xglobal)) {
			icvalue[Indices::PVId_T] = property.parameter.get_T0_ice(); 
			icvalue[Indices::PVId_P] = property.parameter.get_P0(xglobal[0]); 
		
		// outside of the ice cube 
		} else {
			icvalue[Indices::PVId_T] = property.parameter.get_T0(); 
			icvalue[Indices::PVId_P] = property.parameter.get_P0(xglobal[0]); 
		}

		return icvalue;
	}

	//! get a reference to the grid view
	inline const GV& getGridView() {
		return gv;
	}
};
