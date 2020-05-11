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

#include "LMDiffusion.h"

registerMooseObject("LemurApp", LMDiffusion);

InputParameters
LMDiffusion::validParams()
{
  InputParameters params = ADKernelGrad::validParams();
  params.addClassDescription("Diffusion kernel with diffusivity coefficient.");
  params.addRangeCheckedParam<Real>(
      "diffusivity", "diffusivity>0.0", "The diffusivity coefficient.");
  return params;
}

LMDiffusion::LMDiffusion(const InputParameters & parameters)
  : ADKernelGrad(parameters), _diffusivity(getParam<Real>("diffusivity"))
{
}

ADRealVectorValue
LMDiffusion::precomputeQpResidual()
{
  return _diffusivity * _grad_u[_qp];
}