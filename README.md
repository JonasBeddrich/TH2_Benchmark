Warning: Due to the time dependend boundary conditions the benchmark and the snythetic valley data do not operate on the same independent structure (localoperator, postprocess, ...)
Planned to clarfiy on the 23th of June. 
To run ensure that the correct include_problem.hh is included and that USE_Yasp is defined for the benchmark and USE_MESH_GENERATOR for the synthetic valley data. 

# Soilfreezing/permafrost benchmark2D 
## How to run the code: 
0. Clone this repository in **<your_dune_folder>/SoilFreezing/src/permafrost**. 
1. Add **add_executable("permafrost" permafrost/main.cc)** to the **CMakeLists.txt** file in the **src** folder. The name of the executable is currently _permafrost_ and can be changed here. 
2. Create the following folder hierarchy in the **src** folder: **outputs/permafrost/benchmark2D/benchmark_simulations**. Currently that's the folder where the results will be stored. This path can be adjusted in **permafrost/main.cc** and **permafrost/benchmark2D/inputs/sample.ini**. 
3. Run _make_ in **<your_dune_folder>/SoilFreezing/release-build/src**. 
4. Run __./permafrost benchmark2D sample.ini__.  

## Benchmark specifics: 
The different relevant variations of the imposed head can be changed in "permafrost/benchmark2D/inputs/sample.ini"