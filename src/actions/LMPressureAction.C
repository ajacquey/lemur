/******************************************************************************/
/*                            This file is part of                            */
/*                       LEMUR, a MOOSE-based application                     */
/*          muLtiphysics of gEomaterials using MUltiscale Rheologies          */
/*                                                                            */
/*                  Copyright (C) 2019 by Antoine B. Jacquey                  */
/*             GFZ Potsdam, German Research Centre for Geosciences            */
/*                                                                            */
/*            Licensed under GNU Lesser General Public License v2.1           */
/*                       please see LICENSE for details                       */
/*                 or http://www.gnu.org/licenses/lgpl.html                   */
/******************************************************************************/

#include "LMPressureAction.h"
#include "Factory.h"
#include "FEProblem.h"
#include "Conversion.h"

registerMooseAction("LemurApp", LMPressureAction, "add_bc");

template <>
InputParameters
validParams<LMPressureAction>()
{
  InputParameters params = validParams<Action>();
  params.addClassDescription("Set up pressure boundary conditions.");
  params.addRequiredParam<std::vector<BoundaryName>>(
      "boundary", "The list of boundary IDs from the mesh where the pressure will be applied.");
  params.addRequiredParam<std::vector<NonlinearVariableName>>(
      "displacements",
      "The displacements appropriate for the simulation geometry and coordinate system");
  params.addParam<Real>("value", 1.0, "Value of the pressure applied.");
  params.addParam<FunctionName>("function", "Function giving the velocity applied.");
  return params;
}

LMPressureAction::LMPressureAction(const InputParameters & params) : Action(params) {}

void
LMPressureAction::act()
{
  const std::string kernel_name = "LMPressureBC";

  std::vector<NonlinearVariableName> displacements =
      getParam<std::vector<NonlinearVariableName>>("displacements");

  // Create pressure BCs
  for (unsigned int i = 0; i < displacements.size(); ++i)
  {
    // Create unique kernel name for each of the components
    std::string unique_kernel_name = kernel_name + "_" + _name + "_" + Moose::stringify(i);

    InputParameters params = _factory.getValidParams(kernel_name);
    params.applyParameters(parameters());
    params.set<bool>("use_displaced_mesh") = true;
    params.set<unsigned int>("component") = i;
    params.set<NonlinearVariableName>("variable") = displacements[i];

    _problem->addBoundaryCondition(kernel_name, unique_kernel_name, params);
  }
}
