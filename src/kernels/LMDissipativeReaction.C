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

#include "LMDissipativeReaction.h"

registerMooseObject("LemurApp", LMDissipativeReaction);

InputParameters
LMDissipativeReaction::validParams()
{
  InputParameters params = ADKernelValue::validParams();
  params.addClassDescription(
      "Kernel for mechanical dissipative reaction in a diffusion-reaction equation.");
  params.addParam<Real>("coefficient", 1.0, "The coefficient in front of the dissipative term.");
  return params;
}

LMDissipativeReaction::LMDissipativeReaction(const InputParameters & parameters)
  : ADKernelValue(parameters),
    _alpha(getParam<Real>("coefficient")),
    _plastic_strain_incr(getADMaterialProperty<RankTwoTensor>("plastic_strain_increment"))
{
}

ADReal
LMDissipativeReaction::precomputeQpResidual()
{
  ADReal gamma_vp =
      (_dt != 0.0) ? std::sqrt(0.5) * _plastic_strain_incr[_qp].deviatoric().L2norm() / _dt : 0.0;
  return -_alpha * gamma_vp;
}