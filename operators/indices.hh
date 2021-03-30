class Indices {

public:

	// Primary Variables 
	const static int numOfPVs 	= 2;
	const static int PVId_P 	= 0; 		// pressure 
	const static int PVId_T 	= 1; 		// temperature 

	// Equations 
	const static int numOfEqs 	= numOfPVs; 
	const static int Eq_flux 	= 0; 
	const static int Eq_heat 	= 1; 
	
	// Boundary Conditions 
	const static int numOfBCs	= 2;
	const static int BC_flux	= 0;
	const static int BC_heat	= 1;

	// flux 
	const static int dirichletP = 0; 
	const static int neumannf 	= 1; 

	// heat 
	const static int dirichletT = 0; 
	const static int neumannQ 	= 1; 

	std::string getPVName(int tag) const {
		std::string name; 
		switch(tag){
			case PVId_P: 
				name = "P"; 
				break;
			case PVId_T: 
				name = "T"; 
				break;
		}
		return name; 
	}
};
