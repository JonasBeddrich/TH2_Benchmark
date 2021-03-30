/** \brief A function for initial values of P
 */
template<typename GV, typename Properties, typename PGMap>
class P_Initial
		: public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1> >, P_Initial<GV,Properties,PGMap> >
{
private:
	  const GV& gv;
	  const Properties& property;
	  const PGMap& pgmap;
	  ProblemInitialConditions<GV,Properties,PGMap> icvalue;
	  const static int dim = GV::dimension;

public:

  //! construct from grid view
  P_Initial (const GV& gv_, const Properties& property_, const PGMap& pgmap_)
  : gv( gv_ ),
	property(property_),
	pgmap(pgmap_),
	icvalue(gv_,property_,pgmap_)
  {}

  //! evaluate extended function on element
  inline void evaluate (const typename GV::Traits::template Codim<0>::Entity& element,
                        const Dune::FieldVector<double,dim>& xlocal,
                        double& y) const
  {
    y /*ndim*/ = icvalue.evaluate(element,xlocal)[Indices::PVId_P] ;
    return;
  }
  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};


/** \brief A function for initial values of T
 */
template<typename GV, typename Properties, typename PGMap>
class T_Initial
		: public Dune::PDELab::GridFunctionBase<Dune::PDELab::GridFunctionTraits<GV,double,1,Dune::FieldVector<double,1> >, T_Initial<GV,Properties,PGMap> >
{
private:
	  const GV& gv;
	  const Properties& property;
	  const PGMap& pgmap;
	  ProblemInitialConditions<GV,Properties,PGMap> icvalue;
	  const static int dim = GV::dimension;

public:

  //! construct from grid view
  T_Initial (const GV& gv_, const Properties& property_, const PGMap& pgmap_)
  : gv( gv_ ),
	property(property_),
	pgmap(pgmap_),
	icvalue(gv_,property_,pgmap_)
  {}

  //! evaluate extended function on element
  inline void evaluate (const typename GV::Traits::template Codim<0>::Entity& element,
                        const Dune::FieldVector<double,dim>& xlocal,
                        double& y) const
  {
    y /*ndim*/ = icvalue.evaluate(element,xlocal)[Indices::PVId_T] ;
    return;
  }
  //! get a reference to the grid view
  inline const GV& getGridView () {return gv;}
};