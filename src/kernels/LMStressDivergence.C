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

#include "LMStressDivergence.h"
#include "libmesh/quadrature.h"

registerMooseObject("LemurApp", LMStressDivergence);

InputParameters
LMStressDivergence::validParams()
{
  InputParameters params = ADKernel::validParams();
  params.addClassDescription("Solid momentum kernel.");
  params.set<bool>("use_displaced_mesh") = false;
  params.addRequiredParam<unsigned int>("component",
                                        "An integer corresponding to the direction "
                                        "the variable this kernel acts in (0 for x, "
                                        "1 for y, 2 for z).");
  params.addRangeCheckedParam<Real>(
      "density", 0.0, "density >= 0.0", "The density of the material.");
  params.addParam<RealVectorValue>("gravity", RealVectorValue(), "The gravity vector.");
  return params;
}

LMStressDivergence::LMStressDivergence(const InputParameters & parameters)
  : ADKernel(parameters),
    _component(getParam<unsigned int>("component")),
    _rho(getParam<Real>("density")),
    _gravity(getParam<RealVectorValue>("gravity")),
    _stress(getADMaterialProperty<RankTwoTensor>("stress"))
{
}

ADReal
LMStressDivergence::computeQpResidual()
{
  RealVectorValue grav_term = -_rho * _gravity;

  return _stress[_qp].row(_component) * _grad_test[_i][_qp] +
         grav_term(_component) * _test[_i][_qp];
}
