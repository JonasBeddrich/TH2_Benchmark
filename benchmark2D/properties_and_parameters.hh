template<typename GV, typename PTree, typename CharacteristicValues>
class Parameters {
private:
	const GV &gv;
	const PTree &ptree;
	const CharacteristicValues &characteristicValue;
	double Xc_x;
	double Xc_t;
	double Xc_P;
	double Xc_T;
	double Xc_KH;
	double Xc_KT;
	double Xc_lambdaR;
	double Xc_gravity;
	double *time;
	double *dt;

	const double eps = 1.e-6;
	const static int dim = GV::dimension;

	double LX; 

	// parameters 
	double beta; // compressibility 
	double gravitation; 

	// inital values 
	double T0_ice; 
	double T0; 

	// hydraulic head 
	double hh; 
	double h0; 

public:

	//! constructor
	Parameters(const GV &gv_, const PTree &ptree_,
			const CharacteristicValues &characteristicValue_,
			double *time_/*ndim*/, double *dt_/*ndim*/) :
			gv(gv_), ptree(ptree_), characteristicValue(characteristicValue_), time(
					time_), dt(dt_) {

		// Grid 
		LX = ptree.get("grid.yasp.LX", (double) 3.); 

		// parameters 
		beta = ptree.get("parameters.beta", (double) 1.e-8); 
		gravitation = ptree.get("parameters.gravitation", (double) 9.81); 

		// initial values 
		T0_ice = ptree.get("initial.T0_ice", (double) 268.15); 
		T0 = ptree.get("initial.T0", (double) 278.15);
		
		// hydraulic gradient 
		hh = ptree.get("benchmark.hydraulic_head", (double) 0.09); 
		h0 = ptree.get("benchmark.h0", (double) 1); 
	}

	double get_P_right() const {
		return h0; 
	}

	double get_P_left() const {
		return LX * hh + h0; 
	}

	double get_T0_ice () const {
		return T0_ice; 
	}

	double get_T0 () const {
		return T0; 
	}

	double get_P0 (double x) const {
		return h0 + (LX-x) * hh; 
	}

	double get_compressibility() const {
		return beta; 
	}

	double get_gravitation() const {
		return gravitation; 
	}
};

template<typename PGMap, typename PTree,typename CharacteristicValues>
class MeshParameters {
private:
	const PGMap &pgmap;
	const PTree &ptree;
	const CharacteristicValues &characteristicValue;
	double *time;
	double *dt;

	const double eps = 1.e-6;

	double Xc_x;
	double Xc_t;

	double Zmax;
	double Xmax;
	double nZ;
	double nX;

public:

	MeshParameters(const PGMap &pgmap_, const PTree &ptree_,
			const CharacteristicValues &characteristicValue_,
			double *time_/*ndim*/, double *dt_/*ndim*/) :
			pgmap(pgmap_), ptree(ptree_), characteristicValue(
					characteristicValue_), time(time_), dt(dt_) {
		Xc_x = characteristicValue.get_x();
		Xc_t = characteristicValue.get_t();
		Zmax = ptree.get("grid.yasp.LZ", (double) 10.) / Xc_x; /*ndim*/
		nZ = ptree.get("grid.yasp.NZ", (int) 500);
		Xmax = ptree.get("grid.yasp.LX", (double) 10.) / Xc_x; /*ndim*/
		nX = ptree.get("grid.yasp.NX", (double) 10.) / Xc_x;
	}

	bool isBottomBoundary(
			Dune::FieldVector<double, DIMENSION> globalPos/*ndim*/) const {
		if (globalPos[1] < 0. + eps) {
			return true;
		} else
			return false;
	}

	bool isTopBoundary(Dune::FieldVector<double, DIMENSION> globalPos/*ndim*/) const {
		if (globalPos[1] > Zmax/*ndim*/- eps) {
			return true;
		} else
			return false;
	}

	bool isLeftBoundary(
			Dune::FieldVector<double, DIMENSION> globalPos/*ndim*/) const {
		if (globalPos[0] < 0. + eps) {
			return true;
		} else
			return false;
	}

	bool isRightBoundary(
			Dune::FieldVector<double, DIMENSION> globalPos/*ndim*/) const {
		if (globalPos[0] > Xmax - eps) {
			return true;
		} else
			return false;
	}

	bool isIceCube(Dune::FieldVector<double, DIMENSION> globalPos /*ndim*/) const {
		if (globalPos[0] > 0.833333 && globalPos[0] < 1.166667
				&& globalPos[1] > 0.333333 && globalPos[1] < 0.666667) {
			return true;
		} else {
			return false;
		}
	}
};

/*************************************************************************************************************************
 * PHASE PROPERTIES
 * 1. AIR
 * 2. WATER
 * 3. ICE
 * 4. SOIL
 *************************************************************************************************************************/

template<typename CharacteristicValues>
class Air {
private:
	const CharacteristicValues &characteristicValue;
	double Xc_mu;
	double Xc_rho;
	double Xc_kth;
	double Xc_H;
public:

	//! constructor
	Air(const CharacteristicValues &characteristicValue_) :
			characteristicValue(characteristicValue_) {
		Xc_mu = characteristicValue.get_viscosity();
		Xc_rho = characteristicValue.get_density();
		Xc_kth = characteristicValue.get_thermalconductivity();
		Xc_H = characteristicValue.get_specificheat();
	}

	double Density() const {
		double rho = 1.284;
		return rho; 
	}

	double DynamicViscosity() const {
		double mu = 1.725 * 1.0e-5;
		return mu; 
	}

	double ThermalConductivity() const {
		double kth = 0.02428;
		return kth; 
	}

	double Hp() const {
		double H = 1.0038 * 1.e3;
		return H; 
	}

	double Hv() const {
		double H = 0.7167 * 1.e3;
		return H;
	}
};

template<typename CharacteristicValues>
class Water {
private:
	const CharacteristicValues &characteristicValue;
	double Xc_P;
	double Xc_mu;
	double Xc_rho;
	double Xc_kth;
	double Xc_H;
public:

	//! constructor
	Water(const CharacteristicValues &characteristicValue_) :
			characteristicValue(characteristicValue_) {
		Xc_P = characteristicValue.get_P();
		Xc_mu = characteristicValue.get_viscosity();
		Xc_rho = characteristicValue.get_density();
		Xc_kth = characteristicValue.get_thermalconductivity();
		Xc_H = characteristicValue.get_specificheat();
	}
	
	double Density() const {
		double rho = 1000;
		return rho; 
	}

	double DynamicViscosity() const {
		double mu = 0.001793;
		return mu; 
	}

	double ThermalConductivity() const {
		double kth = 0.6;
		return kth; 
	}

	double Hp() const {
		double H = 4182;
		return H; 
	}

	double Hv() const {
		double H = Hp(); 
		return H;
	}
};

template<typename CharacteristicValues, typename Parameters>
class Ice {
private:
	const CharacteristicValues &characteristicValue;
	const Parameters &param;
	double Xc_rho;
	double Xc_kth;
	double Xc_H;
	double Xc_T;
	double Xc_P;
public:

	//! constructor
	Ice(const CharacteristicValues &characteristicValue_,
			const Parameters &param_) :
			characteristicValue(characteristicValue_), param(param_) {
		Xc_rho = characteristicValue.get_density();
		Xc_kth = characteristicValue.get_thermalconductivity();
		Xc_H = characteristicValue.get_specificheat();
		Xc_T = characteristicValue.get_T();
		Xc_P = characteristicValue.get_P();
	}

	double Density() const {
		double rho = 920;
		return rho;
	}

	double ThermalConductivity() const {
		double kth = 2.14;
		return kth;
	}

	double Hp() const {
		double H = 2060;
		return H;
	}

	double Hv() const {
		double H = Hp();
		return H;
	}

	double LatentHeat() const {
		double L = 334000.;
		return L; 
	}

	double RegelationCoefficient() const {
		double R = 0.; 
		return R; 
	}

	double FreezingPointTemperature() const {
		double Tf = 273.15; 
		return Tf;
	}

	double FreezingCurve(double T) const {
		double Sw_res = 0.05; 
		double W = 0.5; 

		double tmp = (T - FreezingPointTemperature()) / W; 
		tmp *= tmp; 
		double SFC = Sw_res + (1 - Sw_res) * std::exp(-tmp); 

		if (T >= FreezingPointTemperature()){
			return 1; 
		} else {
			return SFC;
		}
	}

	double d_dT_FreezingCurve(double T) const {
		double Sw_res = 0.05; 
		double W = 0.5; 

		// prefactor 
		double pf = -(1-Sw_res); 
		pf *= 2 * (T - FreezingPointTemperature()) / W / W; 
		
		double tmp = (T - FreezingPointTemperature()) / W; 
		tmp *= tmp; 

		if (T >= FreezingPointTemperature()){
			return 0; 
		} else {
			return pf * std::exp(-tmp); 
		}
	}
};

template<typename GV, typename PGMap, typename Parameters,
		typename CharacteristicValues>
class Soil {
private:
	const GV &gv;
	const PGMap &pgmap;
	const Parameters &parameter;
	const CharacteristicValues &characteristicValue;
	double Xc_K;
	double Xc_mu;
	double Xc_rho;
	double Xc_kth;
	double Xc_H;
	double Xc_P;
	double Xc_T;
	const static int dim = GV::dimension;

public:

	//! construct from grid view
	Soil(const GV &gv_, const PGMap &pgmap_, const Parameters &parameter_,
			const CharacteristicValues &characteristicValue_) :
			gv(gv_), pgmap(pgmap_), parameter(parameter_), characteristicValue(
					characteristicValue_) {
		Xc_K = characteristicValue.get_permeability();
		Xc_mu = characteristicValue.get_viscosity();
		Xc_rho = characteristicValue.get_density();
		Xc_kth = characteristicValue.get_thermalconductivity();
		Xc_H = characteristicValue.get_specificheat();
		Xc_P = characteristicValue.get_P();
		Xc_T = characteristicValue.get_T();
	}

	double Porosity() const {
		double por = 0.37; 
		return por;
	}

	double HydraulicPermeability() const {
		double K = 1.3e-10; 
		return K;/*ndim*/
	}

	Dune::FieldVector<double, dim> HydraulicPermeabilityVector() const {
		Dune::FieldVector<double, dim> KH;
		KH[0] = HydraulicPermeability();
		KH[1] = HydraulicPermeability();
		return KH; /*ndim*/
	}

	double ThermalPermeability() const {
		double K = 0; 
		return K;/*ndim*/
	}

	// vector coefficient
	Dune::FieldVector<double, dim> ThermalPermeabilityVector() const {
		Dune::FieldVector<double, dim> KT;
		KT[0] = ThermalPermeability();
		KT[1] = ThermalPermeability();
		return KT; /*ndim*/
	}

	double Density() const {
		double rho = 2650.0;
		return rho; 
	}

	double ThermalConductivity() const {
		double kth = 9.0;
		return kth; 
	}

	double Hp() const {
		double H = 835;
		return H; 
	}

	double Hv() const {
		double H = Hp();
		return H; 
	}

	//! get a reference to the grid view
	inline const GV& getGridView() {
		return gv;
	}
};

template<typename GV, typename PGMap, typename Parameters,
		typename CharacteristicValues, typename Soil>
class HydraulicProperties {
private:
	const GV &gv;
	const PGMap &pgmap;
	const Parameters &parameter;
	const CharacteristicValues &characteristicValue;
	double Xc_P;
	double Xc_x;
	const Soil &soil;

	const static int dim = GV::dimension;

	const int numOfParams = 7;
	const int id_Pentry = 0;
	const int id_lambda = 1;
	const int id_Swr = 2;
	const int id_Sgr = 3;
	const int id_beta = 4;
	const int id_surface_tension_ratio = 5;

public:

	//! construct from grid view
	HydraulicProperties(const GV &gv_, const PGMap &pgmap_,
			const Parameters &parameter_,
			const CharacteristicValues &characteristicValue_, const Soil &soil_) :
			gv(gv_), pgmap(pgmap_), parameter(parameter_), characteristicValue(
					characteristicValue_), soil(soil_) {
		Xc_x = characteristicValue.get_x();
		Xc_P = characteristicValue.get_P();
	}

	//! get a reference to the grid view
	inline const GV& getGridView() {
		return gv;
	}

	double krw (double porosity, double Sw) const {
		double Omega = 50; 
		double krw_min = 1.e-6; 
		double krw = std::pow(10, -Omega * porosity * (1-Sw)); 
		return std::max (krw, krw_min); 
	}
};

template<typename GV, typename PTree, typename PGMap>
class Properties {
private:
	const GV &gv;
	const PTree &ptree;
	const PGMap &pgmap;
	double *time;
	double *dt;

	const int dim = GV::dimension;

public:

	//PARAMETERS
	Indices index;
	CharacteristicValues<PTree> characteristicValue;
	Parameters<GV, PTree, CharacteristicValues<PTree>> parameter;
	MeshParameters<PGMap, PTree, CharacteristicValues<PTree>> mesh;

	//PHASE PROPERTIES
	Air<CharacteristicValues<PTree>> air;
	Water<CharacteristicValues<PTree>> water;
	Ice<CharacteristicValues<PTree>,
		Parameters<GV, PTree, CharacteristicValues<PTree>> > ice;
	Soil<GV, PGMap,
		 Parameters<GV, PTree, CharacteristicValues<PTree>>,
		 CharacteristicValues<PTree> > soil;

	//CONSTITUTIVE MODELS
	HydraulicProperties<GV, PGMap,
						Parameters<GV, PTree, CharacteristicValues<PTree>>,
						CharacteristicValues<PTree>,
						Soil<GV, PGMap, Parameters<GV, PTree, CharacteristicValues<PTree>>,CharacteristicValues<PTree> > > hydraulicProperty;

	//! construct from grid view
	Properties(const GV &gv_, const PTree &ptree_, const PGMap &pgmap_,
			   double *time_, double *dt_)
	:gv(gv_), ptree(ptree_), pgmap(pgmap_), time(time_), dt(dt_),
	 characteristicValue(ptree_),
	 parameter(gv_, ptree_, characteristicValue, time_, dt_),
	 mesh(pgmap_, ptree_, characteristicValue, time_, dt_),
	 air(characteristicValue),
	 water(characteristicValue),
	 ice(characteristicValue, parameter),
	 soil(gv_, pgmap, parameter, characteristicValue),
	 hydraulicProperty(gv_, pgmap, parameter, characteristicValue, soil)
	{}

	/******************************************************************************/

	void ReportStatistics(std::string file_name, double time /*s*/,
			double dt /*s*/, int total_newton_iterations,
			double clock_time_elapsed /*s*/) {

		std::fstream result;

		if (time == 0.) {
			result.open(file_name, std::fstream::out | std::fstream::trunc);
			result << "time [s]" << '\t' << "dt [s]" << '\t'
					<< "total no. of newton iterations" << '\t'
					<< "clock time [s]" << std::endl;
			result.close();
		}

		result.open(file_name, std::fstream::app);
		double t_new = time + dt;

		result << time << '\t' << dt << '\t' << total_newton_iterations << '\t'
				<< clock_time_elapsed << std::endl;
		result.close();
	}

	/******************************************************************************/

	//! get a reference to the grid view
	inline const GV& getGridView() {
		return gv;
	}
};
