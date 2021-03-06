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

#include "LMTwoVarUpdate.h"

class LMAlphaGammaYield : public LMTwoVarUpdate
{
public:
  static InputParameters validParams();
  LMAlphaGammaYield(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  // virtual void
  // residual(const ADReal & gamma_v, const ADReal & gamma_d, ADReal & resv, ADReal & resd)
  // override; virtual void jacobian(const ADReal & gamma_v,
  //                       const ADReal & gamma_d,
  //                       ADReal & jacvv,
  //                       ADReal & jacdd,
  //                       ADReal & jacvd,
  //                       ADReal & jacdv) override;
  virtual ADReal yieldFunction(const ADReal & gamma_v, const ADReal & gamma_d) override;
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
  virtual void preReturnMap() override;
  virtual void postReturnMap(const ADReal & gamma_v, const ADReal & gamma_d) override;
  virtual ADRankTwoTensor reformPlasticStrainTensor(const ADReal & gamma_v,
                                                    const ADReal & gamma_d) override;
  virtual void updateYieldParametersDerivV(ADReal & dA, ADReal & dB);
  virtual void updateDissipativeStress(const ADReal & gamma_v,
                                       const ADReal & gamma_d,
                                       ADReal & chi_v,
                                       ADReal & chi_d) override;
  virtual void updateYieldParameters(const ADReal & gamma_v);
  virtual ADReal dyieldFunctiondVol(const ADReal & chi_v, const ADReal & chi_d);
  virtual ADReal dyieldFunctiondDev(const ADReal & chi_v, const ADReal & chi_d);
  virtual ADReal d2yieldFunctiondVol2(const ADReal & chi_v, const ADReal & chi_d);
  virtual ADReal d2yieldFunctiondVolDev(const ADReal & chi_v, const ADReal & chi_d);
  virtual ADReal d2yieldFunctiondDevVol(const ADReal & chi_v, const ADReal & chi_d);
  virtual ADReal d2yieldFunctiondDev2(const ADReal & chi_v, const ADReal & chi_d);
  virtual ADReal dyieldFunctiondA(const ADReal & chi_v, const ADReal & chi_d);
  virtual ADReal dyieldFunctiondB(const ADReal & chi_v, const ADReal & chi_d);
  virtual ADReal d2yieldFunctiondVolA(const ADReal & chi_v, const ADReal & chi_d);
  virtual ADReal d2yieldFunctiondVolB(const ADReal & chi_v, const ADReal & chi_d);
  virtual ADReal d2yieldFunctiondDevA(const ADReal & chi_v, const ADReal & chi_d);
  virtual ADReal d2yieldFunctiondDevB(const ADReal & chi_v, const ADReal & chi_d);

  const Real _phi;
  const Real _pcr0;
  const Real _alpha;
  const Real _gamma;
  const Real _L;
  const bool _has_hardening;
  ADMaterialProperty<Real> * _intnl;
  const MaterialProperty<Real> * _intnl_old;
  Real _M;

  ADReal _pressure_tr;
  ADReal _eqv_stress_tr;

  ADReal _pcr_tr;
  ADReal _one_on_A;
  ADReal _one_on_B;

  ADReal _chi_v_tr;
  ADReal _chi_d_tr;

  ADReal _pcr;
};