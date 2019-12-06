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

#include "LMTwoVarUpdate.h"

#define usingAlphaGammaYieldMembers                                                                \
  usingTwoVarUpdateMembers;                                                                        \
  using LMAlphaGammaYield<compute_stage>::_pressure_tr;                                            \
  using LMAlphaGammaYield<compute_stage>::_eqv_stress_tr;                                          \
  using LMAlphaGammaYield<compute_stage>::_pcr_tr;                                                 \
  using LMAlphaGammaYield<compute_stage>::_pcr;                                                    \
  using LMAlphaGammaYield<compute_stage>::_M;                                                      \
  using LMAlphaGammaYield<compute_stage>::_L;                                                      \
  using LMAlphaGammaYield<compute_stage>::_alpha;                                                  \
  using LMAlphaGammaYield<compute_stage>::_gamma;                                                  \
  using LMAlphaGammaYield<compute_stage>::_one_on_A;                                               \
  using LMAlphaGammaYield<compute_stage>::_one_on_B;

template <ComputeStage>
class LMAlphaGammaYield;

declareADValidParams(LMAlphaGammaYield);

template <ComputeStage compute_stage>
class LMAlphaGammaYield : public LMTwoVarUpdate<compute_stage>
{
public:
  LMAlphaGammaYield(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
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
  virtual ADReal
  calculateProjection(const ADReal & chi_v, const ADReal & chi_d, ADReal & chi_v0, ADReal & chi_d0);
  virtual void calculateProjectionDerivV(const ADReal & chi_v,
                                         const ADReal & chi_d,
                                         ADReal & dchi_v0,
                                         ADReal & dchi_d0);
  virtual void calculateProjectionDerivD(const ADReal & chi_v,
                                         const ADReal & chi_d,
                                         ADReal & dchi_v0,
                                         ADReal & dchi_d0);
  virtual void updateYieldParametersDerivV(ADReal & dA, ADReal & dB);
  virtual ADReal
  calculateDirection(const ADReal & chi_v, const ADReal & chi_d, ADReal & ev, ADReal & ed);
  virtual void updateDissipativeStress(const ADReal & gamma_v,
                                       const ADReal & gamma_d,
                                       ADReal & chi_v,
                                       ADReal & chi_d);
  virtual void updateYieldParameters(const ADReal & gamma_v);

  const Real _phi;
  const Real _pcr0;
  const Real _alpha;
  const Real _gamma;
  const Real _L;
  const bool _has_hardening;
  ADMaterialProperty(Real) * _intnl;
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

  usingTwoVarUpdateMembers;
};