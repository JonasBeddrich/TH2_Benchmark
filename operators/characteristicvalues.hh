template<typename PTree>
class CharacteristicValues
{
private:
	const PTree& ptree;

	double x_c; 		//space
	double t_c; 		//time
	double P_c; 		//pressure
	double T_c; 		//temperature
	double permeability_c;
	double density_c;
	double viscosity_c;
	double specificheat_c;
	double thermalconductivity_c;

public:

	//! constructor
	CharacteristicValues (const PTree& ptree_)
	:ptree(ptree_)
	{
		x_c = ptree.get("characteristic_value.length",(double)1.); /*m*/ //space
		t_c = ptree.get("characteristic_value.time",(double)1.); /*s*/ 	//time
		P_c = ptree.get("characteristic_value.pressure",(double)1.); /*Pa*/ //pressure
		T_c = ptree.get("characteristic_value.temperature",(double)1.); /*K*/ //temperature
		permeability_c = ptree.get("characteristic_value.permeability",(double)1.); /*m^2*/ //permeability
		density_c = ptree.get("characteristic_value.density",(double)1.); /*kg/m^3*/ //density
		viscosity_c	= ptree.get("characteristic_value.dynamic_viscosity",(double)1.); /*Pa.s*/ //dyn. viscosity
		specificheat_c = ptree.get("characteristic_value.speacific_heat",(double)1.); /*J/kg/K*/ //specific heat capacity
		thermalconductivity_c = ptree.get("characteristic_value.thermal_conductivity",(double)1.); /*W/m/K*/ //thermal conductivity
	}

	double get_x()const{
		return x_c;
	}

	double get_t()const{
		return t_c;
	}

	double get_P()const{
		return P_c;
	}

	double get_T()const{
		return T_c;
	}

	double get_permeability()const{
		return permeability_c;
	}

	double get_viscosity()const{
		return viscosity_c;
	}

	double get_density()const{
		return density_c;
	}

	double get_specificheat()const{
		return specificheat_c;
	}

	double get_thermalconductivity()const{
		return thermalconductivity_c;
	}

	double get_gravity()const{
		return P_c/(density_c*x_c);
	}

	double get_thermalcoefficient()const{
		/*m^2*s^-1*K^-1*/
		return P_c*permeability_c/(T_c*viscosity_c);
	}

	double get_source_mass_factor()const{
		return t_c/density_c;
	}

	double get_source_heat_factor()const{
		return t_c/(density_c*specificheat_c*T_c);
	}

	double get_convective_mass_factor()const{
		return ( permeability_c/(x_c*x_c) ) * ( P_c*t_c/viscosity_c );
	}

	double get_convective_heat_factor()const{
		return ( permeability_c/(x_c*x_c) ) * ( P_c*t_c/viscosity_c );
	}

	double get_diffusive_heat_factor()const{
		return ( t_c*thermalconductivity_c ) / ( density_c*specificheat_c*x_c*x_c ) ;
	}

};
