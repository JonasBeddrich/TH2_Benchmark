/*********************************************************
 * EVALUATE OUTPUT VARIABLES
 *********************************************************/

template< class GV,
		  class Params,
		  class Evaluation_Pa,
		  class Evaluation_Sw,
		  class Evaluation_T,
		  class Evaluation_por,
		  class Evaluation_Ci,
		  class GFS_PP, typename U_pp>
class PostProcess{
private:
	const GV&  gv;
	const Params& param;
	Evaluation_Pa  *evaluation_Pa;
	Evaluation_Sw  *evaluation_Sw;
	Evaluation_T   *evaluation_T;
	Evaluation_por *evaluation_por;
	Evaluation_Ci  *evaluation_Ci;
	GFS_PP gfs_pp;
	U_pp *u_pp;

	double Xc_KT ;
	double Xc_KH ;
	double Xc_mu ;
	double Xc_rho;
	double Xc_H;
	double Xc_P;
	double Xc_T;
	double Xc_t;
	double Xc_x;

	typedef typename GV::Traits::template Codim<0>::Iterator LeafIterator;
    typedef typename GV::IndexSet IndexSet;

public:

	PostProcess(const GV& gv_,
				const Params& param_,
				Evaluation_Pa	*evaluation_Pa_,
				Evaluation_Sw	*evaluation_Sw_,
				Evaluation_T	*evaluation_T_,
				Evaluation_por	*evaluation_por_,
				Evaluation_Ci	*evaluation_Ci_,
				GFS_PP	 gfs_pp_,
				U_pp	*u_pp_ )
	: gv(gv_),
	  param(param_),
	  evaluation_Pa(evaluation_Pa_),
	  evaluation_Sw(evaluation_Sw_),
	  evaluation_T(evaluation_T_),
	  evaluation_por(evaluation_por_),
	  evaluation_Ci(evaluation_Ci_),
	  gfs_pp(gfs_pp_),
	  u_pp(u_pp_)
	{
		  Xc_KT		= param.characteristicValue.get_thermalcoefficient();
		  Xc_KH		= param.characteristicValue.get_permeability();
		  Xc_mu 	= param.characteristicValue.get_viscosity();
		  Xc_rho 	= param.characteristicValue.get_density();
		  Xc_H		= param.characteristicValue.get_specificheat();
		  Xc_P 		= param.characteristicValue.get_P();
		  Xc_T 		= param.characteristicValue.get_T();
		  Xc_t 		= param.characteristicValue.get_t();
		  Xc_x 		= param.characteristicValue.get_x();
	  }

	virtual ~PostProcess()
	{}

	void evaluate(){

		typedef Dune::PDELab::LocalFunctionSpace< GFS_PP > LFS_PP;
		LFS_PP lfs_pp(gfs_pp);
		typedef Dune::PDELab::LFSIndexCache<LFS_PP> LFSCache_PP;
		LFSCache_PP lfs_cache_pp(lfs_pp);
		typedef typename U_pp::template LocalView<LFSCache_PP> VectorView_PP;
		VectorView_PP u_pp_view( (*u_pp) );

		using RF = typename LFS_PP::template Child<0>::Type::Traits::FiniteElementType::Traits::LocalBasisType::Traits::RangeFieldType;

		// Loop over each volume - That's not how it should be done ... 

		LeafIterator beginElem = gv.template begin< 0 >();
		LeafIterator endElem = gv.template end< 0 >();

		// Iterate over each element
		for ( LeafIterator self = beginElem; self!= endElem; ++self )
		{
			// Reference to cell
	        const auto& cell = *self;
			const IndexSet &indexSet = gv.indexSet();
			int cell_number = indexSet.index(cell);
	        // get geometry
	        auto geo = cell.geometry();
			// dimension
			const auto dim = geo.mydimension;
	        // cell geometry
	        auto ref_el = referenceElement(geo);
	        auto cell_center_local = ref_el.position(0,0);
	        auto cell_volume = geo.volume();

	        RF Pa=0.;
	        evaluation_Pa->evalFunction(cell,cell_center_local,&Pa);
	        RF Sw=0.;
	        evaluation_Sw->evalFunction(cell,cell_center_local,&Sw);
	        RF T=0.;
	        evaluation_T->evalFunction(cell,cell_center_local,&T);
	        RF Ci=0.;
	        evaluation_Ci->evalFunction(cell,cell_center_local,&Ci);
	        RF por=0.;
	        evaluation_por->evalFunction(cell,cell_center_local,&por);

	        lfs_pp.bind(*self);
			lfs_cache_pp.update();
			u_pp_view.bind(lfs_cache_pp);
	        std::vector<double> ul_pp( lfs_pp.size() );
	        for(int i = 0. ; i < lfs_pp.size() ; i++){
	        	ul_pp[i] = 0.;
	        }

	        auto Cs = 1.-Ci;
	        auto Sa = 1.-Sw;
	        auto Pcaw = param.hydraulicProperty.CapillaryPressure( cell,cell_center_local,Sw,por );
			auto Pwsat = param.water.SaturationPressure( T*Xc_T );
			auto Pw = Pa-Pcaw;//std::max( Pa-Pcaw , Pwsat );
	        auto Pciw = param.ice.CryostaticPressure( Pcaw*Xc_P );
	        auto Pi = Pw + Pciw;
	        auto Cpw = param.water.Hp(T*Xc_T,Pw*Xc_P) * T + param.ice.LatentHeat(T*Xc_T);
	        auto KH = param.hydraulicProperty.HydraulicPermeability(cell,cell_center_local,por)[dim-1];
	        auto KT = param.hydraulicProperty.ThermalPermeability(cell,cell_center_local,por)[dim-1];
	        auto kra = param.hydraulicProperty.kra(cell,cell_center_local,Sw );
	        auto krw = param.hydraulicProperty.krw(cell,cell_center_local,Sw);

	    	ul_pp[lfs_pp.child(param.index.SVId_Pa	).localIndex(0)] = Pa*Xc_P ;
	    	ul_pp[lfs_pp.child(param.index.SVId_Pw	).localIndex(0)] = Pw*Xc_P ;
	    	ul_pp[lfs_pp.child(param.index.SVId_Pi	).localIndex(0)] = Pi*Xc_P ;
	    	ul_pp[lfs_pp.child(param.index.SVId_Pcaw).localIndex(0)] = Pcaw*Xc_P ;
	    	ul_pp[lfs_pp.child(param.index.SVId_Pciw).localIndex(0)] = Pciw*Xc_P ;
	    	ul_pp[lfs_pp.child(param.index.SVId_Sa	).localIndex(0)] = Sa ;
	    	ul_pp[lfs_pp.child(param.index.SVId_Sw	).localIndex(0)] = Sw ;
	    	ul_pp[lfs_pp.child(param.index.SVId_Ci	).localIndex(0)] = Ci ;
	    	ul_pp[lfs_pp.child(param.index.SVId_Cs	).localIndex(0)] = Cs ;
	    	ul_pp[lfs_pp.child(param.index.SVId_por	).localIndex(0)] = por ;
	    	ul_pp[lfs_pp.child(param.index.SVId_T	).localIndex(0)] = T*Xc_T ;
	    	ul_pp[lfs_pp.child(param.index.SVId_Cpw	).localIndex(0)] = Cpw*Xc_H*Xc_T ;
	     	ul_pp[lfs_pp.child(param.index.SVId_KH	).localIndex(0)] = KH*Xc_KH ;
	     	ul_pp[lfs_pp.child(param.index.SVId_KT	).localIndex(0)] = KT*Xc_KT ;
	     	ul_pp[lfs_pp.child(param.index.SVId_krw	).localIndex(0)] = krw ;
	     	ul_pp[lfs_pp.child(param.index.SVId_kra	).localIndex(0)] = kra ;

			u_pp_view.write( ul_pp );
			u_pp_view.commit();
			u_pp_view.unbind();

		}//END:iterate over each volume
	}
};