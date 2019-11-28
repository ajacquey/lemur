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

#pragma once

#include "ADKernel.h"

template <ComputeStage>
class LMStressDivergence;
template <typename>
class RankTwoTensorTempl;
typedef RankTwoTensorTempl<Real> RankTwoTensor;
typedef RankTwoTensorTempl<DualReal> DualRankTwoTensor;

declareADValidParams(LMStressDivergence);

template <ComputeStage compute_stage>
class LMStressDivergence : public ADKernel<compute_stage>
{
public:
  LMStressDivergence(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const unsigned int _component;
  const Real _rho;
  const RealVectorValue _gravity;
  const ADMaterialProperty(RankTwoTensor) & _stress;

  usingKernelMembers;
};