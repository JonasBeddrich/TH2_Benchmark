/*
 * include_problem.hh
 *
 *  Created on: Apr 6, 2020
 *      Author: sgupta
 */

#ifndef PERMAFROST_SAMPLE1D_PROBLEM_INCLUDE_PROBLEM_HH_
#define PERMAFROST_SAMPLE1D_PROBLEM_INCLUDE_PROBLEM_HH_

#include "../../extras/ParameterTraits.hh"
#include "../../extras/Evaluation.hh"

#include "../operators/indices.hh"
#include "../operators/characteristicvalues.hh"

/*******************************************/
// PROBLEM SPECIFICATION
#include "properties_and_parameters.hh"
#include "initial_conditions.hh"
#include "../operators/initialvaluefunction.hh"
#include "boundary_conditions.hh"
#include "../operators/localoperator.hh"
#include "../operators/timeoperator.hh"
#include "../operators/postprocess.hh"
#include "driver.hh"


#endif /* PERMAFROST_SAMPLE1D_PROBLEM_INCLUDE_PROBLEM_HH_ */
