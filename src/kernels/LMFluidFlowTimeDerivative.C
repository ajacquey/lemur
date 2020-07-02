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

#include "LMFluidFlowTimeDerivative.h"

registerMooseObject("LemurApp", LMFluidFlowTimeDerivative);

InputParameters
LMFluidFlowTimeDerivative::validParams()
{
  InputParameters params = ADTimeKernel::validParams();
  params.addClassDescription("Time derivative of fluid pressure for poro-mechanics.");
  return params;
}

LMFluidFlowTimeDerivative::LMFluidFlowTimeDerivative(const InputParameters & parameters)
  : ADTimeKernel(parameters),
    _C_biot(getADMaterialProperty<Real>("biot_compressibility")),
    _poro_mech(getADMaterialProperty<Real>("poro_mech"))
{
}

ADReal
LMFluidFlowTimeDerivative::computeQpResidual()
{
  return (_C_biot[_qp] * _u_dot[_qp] + _poro_mech[_qp]) * _test[_i][_qp];
}