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

#include "LMAlphaGammaYield.h"

template <ComputeStage>
class LMDamageAlphaGammaYield;

declareADValidParams(LMDamageAlphaGammaYield);

template <ComputeStage compute_stage>
class LMDamageAlphaGammaYield : public LMAlphaGammaYield<compute_stage>
{
public:
  LMDamageAlphaGammaYield(const InputParameters & parameters);

protected:
  virtual void preReturnMap() override;
  virtual void updateYieldParameters(const ADReal & gamma_v) override;
  virtual void updateYieldParametersDerivV(ADReal & dA, ADReal & dB) override;
  virtual void postReturnMap(const ADReal & gamma_v, const ADReal & gamma_d) override;

  // Coupled variables
  const ADVariableValue & _damage;
  // Damage viscosity
  const Real _eta_a;

  // Properties
  ADMaterialProperty(Real) & _damage_rate;

  usingAlphaGammaYieldMembers;
};