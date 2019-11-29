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

#include "LMMechMaterial.h"

template <ComputeStage>
class LMDamageMechMaterial;
template <typename>
class RankTwoTensorTempl;
typedef RankTwoTensorTempl<Real> RankTwoTensor;
typedef RankTwoTensorTempl<DualReal> DualRankTwoTensor;
template <ComputeStage>
class LMViscoPlasticUpdate;

declareADValidParams(LMDamageMechMaterial);

template <ComputeStage compute_stage>
class LMDamageMechMaterial : public LMMechMaterial<compute_stage>
{
public:
  LMDamageMechMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpElasticGuess() override;

  // Coupled variables
  const ADVariableValue & _damage_dot;
  const VariableValue & _damage_old;

  usingMechMaterialMembers;
};