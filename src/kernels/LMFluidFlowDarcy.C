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

#include "LMFluidFlowDarcy.h"

registerMooseObject("LemurApp", LMFluidFlowDarcy);

InputParameters
LMFluidFlowDarcy::validParams()
{
  InputParameters params = ADKernelGrad::validParams();
  params.addClassDescription("Divergence of Darcy's velocity for the fluid pressure equation.");
  return params;
}

LMFluidFlowDarcy::LMFluidFlowDarcy(const InputParameters & parameters)
  : ADKernelGrad(parameters), _fluid_mob(getMaterialProperty<Real>("fluid_mobility"))
{
}

ADRealVectorValue
LMFluidFlowDarcy::precomputeQpResidual()
{
  return _fluid_mob[_qp] * _grad_u[_qp];
}