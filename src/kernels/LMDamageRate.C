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

#include "LMDamageRate.h"

registerMooseObject("LemurApp", LMDamageRate);

InputParameters
LMDamageRate::validParams()
{
  InputParameters params = ADKernelValue::validParams();
  params.addClassDescription("Damage rate kernel.");
  return params;
}

LMDamageRate::LMDamageRate(const InputParameters & parameters)
  : ADKernelValue(parameters), _damage_rate(getADMaterialProperty<Real>("damage_rate"))
{
}

ADReal
LMDamageRate::precomputeQpResidual()
{
  return -_damage_rate[_qp];
}