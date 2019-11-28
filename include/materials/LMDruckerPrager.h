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

#include "LMSingleVarUpdate.h"

template <ComputeStage>
class LMDruckerPrager;

declareADValidParams(LMDruckerPrager);

template <ComputeStage compute_stage>
class LMDruckerPrager : public LMSingleVarUpdate<compute_stage>
{
public:
  LMDruckerPrager(const InputParameters & parameters);

protected:
  virtual ADReal yieldFunction(const ADReal & gamma_vp) override;
  virtual ADReal yieldFunctionDeriv(const ADReal & gamma_vp) override;
  virtual void preReturnMap() override;
  virtual void postReturnMap() override;
  virtual ADRankTwoTensor reformPlasticStrainTensor(const ADReal & gamma_vp) override;

  const Real _phi;
  const Real _psi;
  const Real _C;
  Real _alpha;
  Real _beta;
  Real _k;

  ADReal _pressure_tr;
  ADReal _eqv_stress_tr;

  usingSingleVarUpdateMembers;
};