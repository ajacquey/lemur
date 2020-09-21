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

#pragma once

#include "LMAlphaGammaYield.h"

class LMDamageAlphaGammaYield : public LMAlphaGammaYield
{
public:
  static InputParameters validParams();
  LMDamageAlphaGammaYield(const InputParameters & parameters);

protected:
  virtual void preReturnMap() override;
  virtual void overStress(const ADReal & gamma_v,
                          const ADReal & gamma_d,
                          ADReal & over_v,
                          ADReal & over_d) override;
  virtual void overStressDerivV(const ADReal & gamma_v,
                                const ADReal & gamma_d,
                                ADReal & over_v_v,
                                ADReal & over_d_v) override;
  virtual void overStressDerivD(const ADReal & gamma_v,
                                const ADReal & gamma_d,
                                ADReal & over_v_d,
                                ADReal & over_d_d) override;
  virtual void updateYieldParameters(const ADReal & gamma_v) override;
  virtual void updateYieldParametersDerivV(ADReal & dA, ADReal & dB) override;
  virtual void postReturnMap(const ADReal & gamma_v, const ADReal & gamma_d) override;

  // Coupled variables
  const ADVariableValue & _damage;
  // // Damage viscosity
  // const Real _eta_a;
  // Dissipation contributions
  const Real _rv;
  const Real _rs;

  // Properties
  ADMaterialProperty<Real> & _damage_rate;
};