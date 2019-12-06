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

registerADMooseObject("LemurApp", LMDamageRate);

defineADValidParams(
    LMDamageRate, ADKernel, params.addClassDescription("Damage rate kernel."););

template <ComputeStage compute_stage>
LMDamageRate<compute_stage>::LMDamageRate(const InputParameters & parameters)
  : ADKernel<compute_stage>(parameters),
    _damage_rate(getADMaterialProperty<Real>("damage_rate"))
{
}

template <ComputeStage compute_stage>
ADReal
LMDamageRate<compute_stage>::computeQpResidual()
{
  return -_damage_rate[_qp] * _test[_i][_qp];
}