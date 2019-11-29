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

#include "LMDamageMechMaterial.h"

registerADMooseObject("LemurApp", LMDamageMechMaterial);

defineADValidParams(
    LMDamageMechMaterial,
    LMMechMaterial,
    params.addClassDescription("Base class calculating the strain and stress of a damaged material.");
    // Coupled variables
    params.addRequiredCoupledVar(
        "damage",
        "The damage variable."););

template <ComputeStage compute_stage>
LMDamageMechMaterial<compute_stage>::LMDamageMechMaterial(const InputParameters & parameters)
  : LMMechMaterial<compute_stage>(parameters),
    // Coupled variables
    _damage_dot(adCoupledDot("damage")),
    _damage_old(coupledValueOld("damage"))
{
}

template <ComputeStage compute_stage>
void
LMDamageMechMaterial<compute_stage>::computeQpElasticGuess()
{
  _elastic_strain_incr[_qp] = _strain_increment[_qp];
  ADReal damage_corr = 1.0 - _damage_dot[_qp] * _dt / (1.0 - _damage_old[_qp]);
  _stress[_qp] = damage_corr * spinRotation(_stress_old[_qp]) +
                 (1.0 - _damage_old[_qp]) * _Cijkl * _strain_increment[_qp];
}