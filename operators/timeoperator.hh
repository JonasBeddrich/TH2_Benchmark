template < class GV, class Params >
class TimeOperator
  : public Dune::PDELab::NumericalJacobianApplyVolume< TimeOperator<GV,Params> >,
    public Dune::PDELab::NumericalJacobianVolume	 < TimeOperator<GV,Params> >,
    public Dune::PDELab::FullVolumePattern,
    public Dune::PDELab::LocalOperatorDefaultFlags,
    public Dune::PDELab::InstationaryLocalOperatorDefaultMethods<double>
{
private:
	const GV& gv;
	const Params& param;
	double *time;
	double *dt;
	double Xc_P;
	double Xc_T;

public:
	  // pattern assembly flags
	  enum { doPatternVolume = true };

	  // residual assembly flags
	  enum { doAlphaVolume = true };

	  typedef typename GV::IndexSet IndexSet;

	  // constructor remembers parameters
	  TimeOperator( const GV& gv_, const Params& param_, double	*time_, double *dt_ )
	  :  gv(gv_), param(param_), time( time_ ), dt( dt_ )
	  {
		  Xc_P 		= param.characteristicValue.get_P();
		  Xc_T 		= param.characteristicValue.get_T();
	  }

	  // volume integral depending on test and ansatz functions
	  template<typename EG, typename LFSU, typename X, typename LFSV, typename R>
	  void alpha_volume (const EG& eg, const LFSU& lfsu, const X& x, const LFSV& lfsv, R& r) const {
			
			const auto& cell = eg.entity();
			const IndexSet &indexSet = gv.indexSet();
			int cell_number = indexSet.index(cell);
			
			auto geo = eg.geometry();
			const auto dim = geo.mydimension;

	        // cell geometry
	        auto ref_el = referenceElement(geo);
	        auto cell_center_local = ref_el.position(0,0);
	        auto cell_volume = geo.volume();

			// compute PVs 
			double P = x(lfsu.child(param.index.PVId_P), 0); 
			double T = x(lfsu.child(param.index.PVId_T), 0); 
			
			// compute SVs 
			double Sw = param.ice.FreezingCurve(T); 
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
			double rhoC_eq = porosity * (Sw * rho_w * C_w + Ci * rho_i * C_i) + (1 - porosity) * rho_s * C_s; 

			// CURRENTLY NOT IN USE 			
			// r.accumulate(lfsu.child(param.index.Eq_heat), 0, tmp*cell_volume); 
	}	
};