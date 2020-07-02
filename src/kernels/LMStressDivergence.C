/******************************************************************************/
/*                            This file is part of                            */
/*                       LEMUR, a MOOSE-based application                     */
/*          muLtiphysics of gEomaterials using MUltiscale Rheologies          */
/*                                                                            */
/*                  Copyright (C) 2020 by Antoine B. Jacquey                  */
/*                    Massachusetts Institute of Technology                   */
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
  params.addCoupledVar("fluid_pressure", 0, "The fluid pressure variable.");
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
    _pf(adCoupledValue("fluid_pressure")),
    _component(getParam<unsigned int>("component")),
    _rho(getParam<Real>("density")),
    _gravity(getParam<RealVectorValue>("gravity")),
    _coupled_pf(isCoupled("fluid_pressure")),
    _stress(getADMaterialProperty<RankTwoTensor>("stress")),
    _biot(_coupled_pf ? &getADMaterialProperty<Real>("biot_coefficient") : nullptr)
{
}

ADReal
LMStressDivergence::computeQpResidual()
{
  RealVectorValue grav_term = -_rho * _gravity;

  ADRealVectorValue stress_row = _stress[_qp].row(_component);
  if (_coupled_pf)
    stress_row(_component) -= (*_biot)[_qp] * _pf[_qp];

  return stress_row * _grad_test[_i][_qp] + grav_term(_component) * _test[_i][_qp];
}
