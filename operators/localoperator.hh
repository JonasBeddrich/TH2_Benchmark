template<class GV, class Params, class BC, class Evaluation_P_old,
		class Evaluation_T_old>
class LocalOperator: 
	public Dune::PDELab::NumericalJacobianApplyVolume<LocalOperator<GV, Params, BC, Evaluation_P_old, Evaluation_T_old>>,
	public Dune::PDELab::NumericalJacobianVolume<LocalOperator<GV, Params, BC, Evaluation_P_old, Evaluation_T_old>>,
	public Dune::PDELab::NumericalJacobianApplySkeleton<LocalOperator<GV, Params, BC, Evaluation_P_old, Evaluation_T_old>>,
	public Dune::PDELab::NumericalJacobianSkeleton<LocalOperator<GV, Params, BC, Evaluation_P_old, Evaluation_T_old>>,
	public Dune::PDELab::NumericalJacobianApplyBoundary<LocalOperator<GV, Params, BC, Evaluation_P_old, Evaluation_T_old>>,
	public Dune::PDELab::NumericalJacobianBoundary<LocalOperator<GV, Params, BC, Evaluation_P_old, Evaluation_T_old>>,

	public Dune::PDELab::FullSkeletonPattern,     // matrix entries skeleton
	public Dune::PDELab::FullVolumePattern,		
	public Dune::PDELab::LocalOperatorDefaultFlags,
	public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double> {

private:
	const GV &gv;
	const Params &param;
	const BC &bc;

	Evaluation_P_old *evaluation_P_old;
	Evaluation_T_old *evaluation_T_old;

	double *time;
	double *dt;
	double Xc_conv_m;
	double Xc_conv_h;
	double Xc_diff_h;
	double Xc_source_m;
	double Xc_source_h;
	double Xc_K;
	double Xc_mu;
	double Xc_rho;
	double Xc_kth;
	double Xc_H;
	double Xc_P;
	double Xc_T;
	double Xc_t;
	double Xc_x;

public:
	// pattern assembly flags
	enum {
		doPatternVolume = true
	};
	enum {
		doPatternSkeleton = true
	};

	// residual assembly flags
	enum {
		doAlphaVolume = true
	};
	enum {
		doAlphaSkeleton = true
	};
	enum {
		doAlphaBoundary = true
	};

	typedef typename GV::IndexSet IndexSet;

	// constructor stores parameters
	LocalOperator(const GV &gv_, const Params &param_, const BC &bc_, 
			Evaluation_P_old *evaluation_P_old_,
			Evaluation_T_old *evaluation_T_old_, double *time_, double *dt_) :
			gv(gv_), param(param_), bc(bc_), evaluation_P_old(evaluation_P_old_), 
			evaluation_T_old(evaluation_T_old_), time(time_), dt(dt_) {

		Xc_conv_m = param.characteristicValue.get_convective_mass_factor();
		Xc_conv_h = param.characteristicValue.get_convective_heat_factor();
		Xc_diff_h = param.characteristicValue.get_diffusive_heat_factor();
		Xc_source_m = param.characteristicValue.get_source_mass_factor();
		Xc_source_h = param.characteristicValue.get_source_heat_factor();

		Xc_K = param.characteristicValue.get_permeability();
		Xc_mu = param.characteristicValue.get_viscosity();
		Xc_rho = param.characteristicValue.get_density();
		Xc_kth = param.characteristicValue.get_thermalconductivity();
		Xc_H = param.characteristicValue.get_specificheat();
		Xc_P = param.characteristicValue.get_P();
		Xc_T = param.characteristicValue.get_T();
		Xc_t = param.characteristicValue.get_t();
		Xc_x = param.characteristicValue.get_x();
	}

	// volume integral depending on test and ansatz functions
	template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	void alpha_volume(const EG &eg, const LFSU &lfsu, const X &x,
			const LFSV &lfsv, R &r) const {
		
		const auto &cell = eg.entity(); 
		const IndexSet &indexSet = gv.indexSet(); 
		int cell_number = indexSet.index(cell); 
		
		auto geo = eg.geometry(); 
		const auto dim = geo.mydimension; 
		
		auto ref_el = referenceElement(geo);
		auto cell_center_local = ref_el.position(0, 0);
		auto cell_center_global = geo.center();
		auto cell_volume = geo.volume();

		double T_old = 0.; 
		double P_old = 0.; 

		// load old values 
		evaluation_T_old->evalFunction(cell, cell_center_local, &T_old); 
		evaluation_P_old->evalFunction(cell, cell_center_local, &P_old); 

		// compute PVs 
		double P = x(lfsu.child(param.index.PVId_P), 0); 
		double T = x(lfsu.child(param.index.PVId_T), 0); 
			
		// compute SVs 
		double Sw = param.ice.FreezingCurve(T); 
		double Sw_old = param.ice.FreezingCurve(T_old); 
		double Ci = 1-Sw; 

		// load parameters 
		double beta = param.parameter.get_compressibility(); 
		double g = param.parameter.get_gravitation(); 
			
		double rho_w = param.water.Density();
		double C_w = param.water.Hp(); 
						
		double rho_i = param.ice.Density();
		double C_i = param.ice.Hp(); 
		double Lf = param.ice.LatentHeat(); 
			
		double rho_s = param.soil.Density();
		double C_s = param.soil.Hp();   
		double porosity = param.soil.Porosity();

		// calculate Parameters 
		double rho_C_eq = porosity * (Sw * rho_w * C_w + Ci * rho_i * C_i) + (1 - porosity) * rho_s * C_s; 
		
		// derivative of the SFC -- it does not matter which one 
		// linear version 
		// double dSFC_dT = (Sw - Sw_old) / (T - T_old + 1.e-24); 

		// analytic version 
		double dSFC_dT = param.ice.dFreezingCurve_dT((T + T_old) / 2); 
		
		// Flow Equation 
		double tmp = 0.; 
		tmp -= Sw * porosity * rho_w * g * beta * (P - P_old) / (*dt);
		tmp += porosity * (rho_i - rho_w) / rho_w * (Sw - Sw_old) / (*dt); 
		r.accumulate(lfsu.child(param.index.Eq_flux), 0, tmp * cell_volume);

		// Heat Equation 
		tmp = 0.; 
		tmp -= rho_C_eq * (T - T_old) / (*dt); 
		tmp -= porosity * rho_i * Lf * dSFC_dT * (T - T_old) / (*dt); 
		r.accumulate(lfsu.child(param.index.Eq_heat), 0, tmp*cell_volume);
	}

	template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	void alpha_skeleton(const IG &ig, const LFSU &lfsu_s, const X &x_s,
			const LFSV &lfsv_s, const LFSU &lfsu_n, const X &x_n,
			const LFSV &lfsv_n, R &r_s, R &r_n) const {

		// References to inside and outside cells
		const auto &cell_inside = ig.inside();
		const auto &cell_outside = ig.outside();
		const IndexSet &indexSet = gv.indexSet();
		int inside_cell_number = indexSet.index(cell_inside);
		int outside_cell_number = indexSet.index(cell_outside);

		// get geometries
		auto geo = ig.geometry();
		const auto dim = geo.mydimension;
		auto geo_inside = cell_inside.geometry();
		auto geo_outside = cell_outside.geometry();

		// cell geometries
		auto ref_el_inside = referenceElement(geo_inside);
		auto ref_el_outside = referenceElement(geo_outside);
		auto inside_cell_center_local = ref_el_inside.position(0, 0);
		auto outside_cell_center_local = ref_el_outside.position(0, 0);
		auto inside_cell_center_global = geo_inside.center();
		auto outside_cell_center_global = geo_outside.center();

		// distance of cell centers
		auto d = outside_cell_center_global;
		d -= inside_cell_center_global;
		auto distance = d.two_norm();

		// face geometry
		auto ref_el = referenceElement(geo);
		auto face_center_local = ref_el.position(0, 0);
		auto face_center_global = geo.center();
		auto face_volume = geo.volume();
		auto normal = ig.unitOuterNormal(face_center_local);

		// compute PVs 
		double P_s = x_s(lfsu_s.child(param.index.PVId_P),0); 
		double P_n = x_n(lfsu_n.child(param.index.PVId_P),0);  

		double T_s = x_s(lfsu_s.child(param.index.PVId_T),0); 
		double T_n = x_n(lfsu_n.child(param.index.PVId_T),0);  

		// gradients 
		double grad_P = (P_n - P_s) / distance; 
		double grad_T = (T_n - T_s) / distance; 		

		// compute SVs 
		double Sw_s = param.ice.FreezingCurve(T_s);
		double Sw_n = param.ice.FreezingCurve(T_n);

		double Ci_s = 1 - Sw_s; 
		double Ci_n = 1 - Sw_n; 

		// load parameters
		double porosity = param.soil.Porosity();  
		double g = param.parameter.get_gravitation(); 

		if (normal[0] == -1){
			double hh = param.parameter.get_hydraulic_gradient(); 
			grad_P += hh; 
		} 
		
		// std::cout << normal << std::endl; 
		double rho_w = param.water.Density(); 
		double rho_i = param.ice.Density(); 
		double rhow_s = param.soil.Density(); 

		double C_w = param.water.Hp();

		// thermal conductivity 
		double lambda_w = param.water.ThermalConductivity(); 
		double lambda_i = param.ice.ThermalConductivity(); 
		double lambda_s = param.soil.ThermalConductivity(); 

		// permeability 
		double kint = param.soil.HydraulicPermeability(); 
		double mu = param.water.DynamicViscosity(); 
		
		double krw_s = param.hydraulicProperty.krw(porosity, Sw_s);
		double krw_n = param.hydraulicProperty.krw(porosity, Sw_n);
		
		double KH_s = krw_s * rho_w * g * kint / mu;
		double KH_n = krw_n * rho_w * g * kint / mu;
		double KH = 2 * KH_s * KH_n / (KH_s + KH_n);
		// double KH = (KH_s + KH_n) / 2;  
		// double KH = std::min(KH_s, KH_n); 

		double lambda_eq_s = porosity * (Sw_s * lambda_w + Ci_s * lambda_i) + (1-porosity) * lambda_s; 
		double lambda_eq_n = porosity * (Sw_n * lambda_w + Ci_n * lambda_i) + (1-porosity) * lambda_s; 
		double lambda_eq = 2 * lambda_eq_s * lambda_eq_n / (lambda_eq_s + lambda_eq_n); 
		// double lambda_eq = (lambda_eq_s + lambda_eq_n) / 2; 
		// double lambda_eq = std::min(lambda_eq_s, lambda_eq_n); 

		double T; 
		// upwinding with regard to pressure 
		// grad_P < 0 means flux from x_s to x_n
		// grad_P = 0 is irrelevant since then the entire term is 0 
		if(grad_P < 0){
			T = T_s;
		} else {
			T = T_n;
		}
		
		// FLOW EQUATION 
		double tmp = 0.; 
		tmp += grad_P * KH; 

		r_s.accumulate(lfsu_s.child(param.index.Eq_flux), 0, +tmp * face_volume);
		r_n.accumulate(lfsu_n.child(param.index.Eq_flux), 0, -tmp * face_volume);

		// HEAT EQUATION 
		tmp = 0.; 
		// term 1  
		tmp += lambda_eq * grad_T; 
		// term 2  
		tmp +=  rho_w * C_w * T * KH * grad_P; 

		r_s.accumulate(lfsu_s.child(param.index.Eq_heat), 0, +tmp * face_volume);
		r_n.accumulate(lfsu_n.child(param.index.Eq_heat), 0, -tmp * face_volume);
	}

	template<typename IG, typename LFSU, typename X, typename LFSV, typename R>
	void alpha_boundary(const IG &ig, const LFSU &lfsu, const X &x,
			const LFSV &lfsv, R &r) const {

		// References to inside and outside cells
		const auto &cell_inside = ig.inside();
		const auto &cell_outside = cell_inside;
		const IndexSet &indexSet = gv.indexSet();
		int inside_cell_number = indexSet.index(cell_inside);

		// get geometries
		auto geo = ig.geometry();
		const auto dim = geo.mydimension;
		auto geo_inside = cell_inside.geometry();

		// cell geometries
		auto ref_el_inside = referenceElement(geo_inside);
		auto inside_cell_center_local = ref_el_inside.position(0, 0);
		auto outside_cell_center_local = inside_cell_center_local;
		auto inside_cell_center_global = geo_inside.center();

		// face geometry
		auto ref_el = referenceElement(geo);
		auto face_center_local = ref_el.position(0, 0);
		auto face_center_global = geo.center();
		auto face_volume = geo.volume();
		auto normal = ig.unitOuterNormal(face_center_local);

		// distance of cell centers
		auto d = geo.global(face_center_local);
		d -= inside_cell_center_global;
		auto distance = d.two_norm();
		
		// load boundary conditions 
		auto bctype = bc.type(ig, face_center_local, (*time) * Xc_t, (*dt) * Xc_t);
		auto bcvalue = bc.value(ig, face_center_local, (*time) * Xc_t, (*dt) * Xc_t);

		// compute PVs 
		double P = x(lfsu.child(param.index.PVId_P),0); 
		double T = x(lfsu.child(param.index.PVId_T),0); 

		// depends on the boundary conditions (already multiplied with n)
		double grad_P; 
		double grad_T; 

		// compute SVs 
		double Sw = param.ice.FreezingCurve(T);
		double Ci = 1 - Sw; 

		// load parameters 
		double porosity = param.soil.Porosity();  
		double g = param.parameter.get_gravitation(); 

		double rho_w = param.water.Density(); 
		double rho_i = param.ice.Density(); 
		double rho_s = param.soil.Density(); 

		double C_w = param.water.Hp();

		// permeability 
		double kint = param.soil.HydraulicPermeability(); 
		double mu = param.water.DynamicViscosity(); 
		double krw = param.hydraulicProperty.krw(porosity, Sw);

		double KH = krw * kint * rho_w * g / mu;
		
		// thermal conductivity
		double lambda_w = param.water.ThermalConductivity(); 
		double lambda_i = param.ice.ThermalConductivity(); 
		double lambda_s = param.soil.ThermalConductivity(); 
		double lambda_eq = porosity * (Sw * lambda_w + Ci * lambda_i) + (1-porosity) * lambda_s; 

		// boundary conditions
		if (bctype[Indices::BC_flux] == Indices::dirichletP) {
			grad_P = (bcvalue[Indices::BC_flux] - P) / distance; 
		} else if (bctype[Indices::BC_flux] == Indices::neumannf) {
			grad_P = bcvalue[Indices::BC_flux]; 
		} 

		if (normal[0] == -1){
			double hh = param.parameter.get_hydraulic_gradient(); 
			grad_P += hh; 
		}

		if (bctype[Indices::BC_heat] == Indices::dirichletT) {
			grad_T = (bcvalue[Indices::BC_heat] - T) / distance;
			T = bcvalue[Indices::BC_heat]; 
		} else if (bctype[Indices::BC_heat] == Indices::neumannQ) {
			grad_T = bcvalue[Indices::BC_heat]; 
		}

		double tmp = 0.;  
 
		// FLOW EQUATION
		tmp += grad_P * KH; 
		r.accumulate(lfsu.child(param.index.Eq_flux), 0, +tmp * face_volume); 

		// HEAT EQUATION 
		tmp = 0.; 
		// term 1  
		tmp += lambda_eq * grad_T; 
		// term 2  
		tmp += rho_w * C_w * T * KH * grad_P; 
		r.accumulate(lfsu.child(param.index.Eq_heat), 0, +tmp * face_volume); 
	}
};