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

#include "LMAlphaGammaYield.h"

registerADMooseObject("LemurApp", LMAlphaGammaYield);

defineADValidParams(
    LMAlphaGammaYield,
    LMTwoVarUpdate,
    params.addClassDescription("Viscoplastic update based on the alpha-gamma yield functions.");
    params.addRequiredRangeCheckedParam<Real>("friction_angle",
                                              "friction_angle > 0.0",
                                              "The friction angle for the critical state line.");
    params.addRequiredRangeCheckedParam<Real>(
        "critical_pressure",
        "critical_pressure > 0.0",
        "The critical pressure of the capped yield.");
    params.addRangeCheckedParam<Real>("alpha", 1.0, "alpha >= 0.0 & alpha <= 1.0", "The alpha parameter for the yield/");
    params.addRangeCheckedParam<Real>("gamma", 1.0, "gamma >= 0.0 & gamma <= 1.0", "The gamma parameter for the yield/");
    params.addRangeCheckedParam<Real>("critical_pressure_hardening",
                                      0.0,
                                      "critical_pressure_hardening >= 0.0",
                                      "The hardening parameter in the exponential for the critical "
                                      "pressure of the capped yield."););

template <ComputeStage compute_stage>
LMAlphaGammaYield<compute_stage>::LMAlphaGammaYield(const InputParameters & parameters)
  : LMTwoVarUpdate<compute_stage>(parameters),
    _phi(getParam<Real>("friction_angle")),
    _pcr0(getParam<Real>("critical_pressure")),
    _alpha(getParam<Real>("alpha")),
    _gamma(getParam<Real>("gamma")),
    _L(getParam<Real>("critical_pressure_hardening")),
    _has_hardening(_L > 0.0),
    _intnl(_has_hardening ? &declareADProperty<Real>("volumetric_plastic_strain") : nullptr),
    _intnl_old(_has_hardening ? &getMaterialPropertyOld<Real>("volumetric_plastic_strain") : nullptr)
{
  _M = std::sqrt(3.0) * std::sin(_phi * libMesh::pi / 180.0);
}

template <ComputeStage compute_stage>
void
LMAlphaGammaYield<compute_stage>::initQpStatefulProperties()
{
  if (_has_hardening)
    (*_intnl)[_qp] = 0.0;
}

template <ComputeStage compute_stage>
ADReal
LMAlphaGammaYield<compute_stage>::yieldFunction(const ADReal & gamma_v, const ADReal & gamma_d)
{

  ADReal chi_v = 0.0, chi_d = 0.0;
  updateDissipativeStress(gamma_v, gamma_d, chi_v, chi_d);

  return Utility::pow<2>(chi_v / _A) + Utility::pow<2>(chi_d / _B) - 1.0;
}

template <ComputeStage compute_stage>
void
LMAlphaGammaYield<compute_stage>::overStress(const ADReal & gamma_v,
                                      const ADReal & gamma_d,
                                      ADReal & over_v,
                                      ADReal & over_d)
{
  ADReal chi_v = 0.0, chi_d = 0.0;
  updateDissipativeStress(gamma_v, gamma_d, chi_v, chi_d);

  ADReal chi_v0 = 0.0, chi_d0 = 0.0;
  calculateProjection(chi_v, chi_d, chi_v0, chi_d0);

  over_v = chi_v - chi_v0;
  over_d = chi_d - chi_d0;
}

template <ComputeStage compute_stage>
void
LMAlphaGammaYield<compute_stage>::overStressDerivV(const ADReal & gamma_v,
                                            const ADReal & gamma_d,
                                            ADReal & over_v_v,
                                            ADReal & over_d_v)
{
  // Dissipative stresses
  ADReal chi_v = 0.0, chi_d = 0.0;
  updateDissipativeStress(gamma_v, gamma_d, chi_v, chi_d);

  // Dissipative stress derivatives
  ADReal dchi_v = -_K * _dt - 0.5 * _gamma * _pcr * _L * _dt;
  ADReal dchi_d = 0.0;

  // Derivative of the projection
  ADReal dchi_v0 = 0.0, dchi_d0 = 0.0;
  calculateProjectionDerivV(chi_v, chi_d, dchi_v0, dchi_d0);

  over_v_v = dchi_v - dchi_v0;
  over_d_v = dchi_d - dchi_d0;
}

template <ComputeStage compute_stage>
void
LMAlphaGammaYield<compute_stage>::overStressDerivD(const ADReal & gamma_v,
                                            const ADReal & gamma_d,
                                            ADReal & over_v_d,
                                            ADReal & over_d_d)
{
  // Dissipative stresses
  ADReal chi_v = 0.0, chi_d = 0.0;
  updateDissipativeStress(gamma_v, gamma_d, chi_v, chi_d);

  // Dissipative stress derivatives
  ADReal dchi_v = 0.0;
  ADReal dchi_d = -3.0 * _G * _dt;

  // Derivative of the projection
  ADReal dchi_v0 = 0.0, dchi_d0 = 0.0;
  calculateProjectionDerivD(chi_v, chi_d, dchi_v0, dchi_d0);

  over_v_d = dchi_v - dchi_v0;
  over_d_d = dchi_d - dchi_d0;
}

template <ComputeStage compute_stage>
void
LMAlphaGammaYield<compute_stage>::preReturnMap()
{
  _pressure_tr = -_stress_tr.trace() / 3.0;
  _eqv_stress_tr = std::sqrt(1.5) * _stress_tr.deviatoric().L2norm();

  _pcr_tr = _pcr0;
  if (_has_hardening)
  { 
    (*_intnl)[_qp] = (*_intnl_old)[_qp];
    _pcr_tr = _pcr0 * std::exp(-_L * (*_intnl_old)[_qp]);
  }

  _A = (1.0 - _gamma) * _pressure_tr + 0.5 * _gamma * _pcr;
  _B = _M * ((1.0 - _alpha) * _pressure_tr + 0.5 * _alpha * _gamma * _pcr);

  _chi_v_tr = _pressure_tr - 0.5 * _gamma * _pcr_tr;
  _chi_d_tr = _eqv_stress_tr;
}

template <ComputeStage compute_stage>
void
LMAlphaGammaYield<compute_stage>::postReturnMap(const ADReal & gamma_v, const ADReal & /*gamma_d*/)
{
  if (_has_hardening)
    (*_intnl)[_qp] = (*_intnl_old)[_qp] + gamma_v * _dt;
}

template <ComputeStage compute_stage>
ADRankTwoTensor
LMAlphaGammaYield<compute_stage>::reformPlasticStrainTensor(const ADReal & gamma_v, const ADReal & gamma_d)
{
  ADRankTwoTensor flow_dir =
      (_eqv_stress_tr != 0.0) ? _stress_tr.deviatoric() / _eqv_stress_tr : ADRankTwoTensor();

  ADRankTwoTensor delta_gamma = 1.5 * gamma_d * _dt * flow_dir;
  delta_gamma.addIa(-gamma_v * _dt / 3.0);

  return delta_gamma;
}

template <ComputeStage compute_stage>
ADReal
LMAlphaGammaYield<compute_stage>::calculateProjection(const ADReal & chi_v,
                                               const ADReal & chi_d,
                                               ADReal & chi_v0,
                                               ADReal & chi_d0)
{
  // Directions
  ADReal ev = 0.0, ed = 0.0;
  calculateDirection(chi_v, chi_d, ev, ed);

  ADReal rho_0 = std::sqrt(1.0 / (Utility::pow<2>(ev / _A) + Utility::pow<2>(ed / _B)));

  chi_v0 = rho_0 * ev;
  chi_d0 = rho_0 * ed;

  return rho_0;
}

template <ComputeStage compute_stage>
void
LMAlphaGammaYield<compute_stage>::calculateProjectionDerivV(const ADReal & chi_v,
                                                     const ADReal & chi_d,
                                                     ADReal & dchi_v0,
                                                     ADReal & dchi_d0)
{
  // Directions
  ADReal ev = 0.0, ed = 0.0;
  ADReal rho_tr = calculateDirection(chi_v, chi_d, ev, ed);
  // Dissipative stress derivative
  ADReal dchi_v = -_K * _dt - 0.5 * _gamma * _pcr * _L * _dt;

  ADReal rho_0 = std::sqrt(1.0 / (Utility::pow<2>(ev / _A) + Utility::pow<2>(ed / _B)));

  dchi_v0 =
      rho_0 / rho_tr *
          (Utility::pow<2>(rho_0) * (1.0 / Utility::pow<2>(_B) - 1.0 / Utility::pow<2>(_A)) *
               Utility::pow<2>(ev) +
           1.0) *
          Utility::pow<2>(ed) * dchi_v +
      Utility::pow<3>(rho_0) *
          (Utility::pow<2>(ev) / Utility::pow<3>(_A) + Utility::pow<2>(ed) / Utility::pow<3>(_B)) *
          ev;
  dchi_d0 =
      rho_0 / rho_tr *
          (Utility::pow<2>(rho_0) * (1.0 / Utility::pow<2>(_B) - 1.0 / Utility::pow<2>(_A)) *
               Utility::pow<2>(ed) -
           1.0) *
          ev * ed * dchi_v +
      Utility::pow<3>(rho_0) *
          (Utility::pow<2>(ev) / Utility::pow<3>(_A) + Utility::pow<2>(ed) / Utility::pow<3>(_B)) *
          ed;
}

template <ComputeStage compute_stage>
void
LMAlphaGammaYield<compute_stage>::calculateProjectionDerivD(const ADReal & chi_v,
                                                     const ADReal & chi_d,
                                                     ADReal & dchi_v0,
                                                     ADReal & dchi_d0)
{
  // Directions
  ADReal ev = 0.0, ed = 0.0;
  ADReal rho_tr = calculateDirection(chi_v, chi_d, ev, ed);
  // Dissipative stress derivative
  ADReal dchi_d = -3.0 * _G * _dt;

  ADReal rho_0 = std::sqrt(1.0 / (Utility::pow<2>(ev / _A) + Utility::pow<2>(ed / _B)));

  dchi_v0 = rho_0 / rho_tr *
            (Utility::pow<2>(rho_0) * (1.0 / Utility::pow<2>(_A) - 1.0 / Utility::pow<2>(_B)) *
                 Utility::pow<2>(ev) -
             1.0) *
            ev * ed * dchi_d;
  dchi_d0 = rho_0 / rho_tr *
            (Utility::pow<2>(rho_0) * (1.0 / Utility::pow<2>(_A) - 1.0 / Utility::pow<2>(_B)) *
                 Utility::pow<2>(ed) +
             1.0) *
            Utility::pow<2>(ev) * dchi_d;
}

template <ComputeStage compute_stage>
ADReal
LMAlphaGammaYield<compute_stage>::calculateDirection(const ADReal & chi_v,
                                              const ADReal & chi_d,
                                              ADReal & ev,
                                              ADReal & ed)
{
  ADReal rho_tr = std::sqrt(Utility::pow<2>(chi_v) + Utility::pow<2>(chi_d));
  ev = (rho_tr != 0.0) ? chi_v / rho_tr : 0.0;
  ed = (rho_tr != 0.0) ? chi_d / rho_tr : 0.0;

  return rho_tr;
}

template <ComputeStage compute_stage>
void
LMAlphaGammaYield<compute_stage>::updateDissipativeStress(const ADReal & gamma_v,
                                                   const ADReal & gamma_d,
                                                   ADReal & chi_v,
                                                   ADReal & chi_d)
{
  // Here we calculate the yield function in the dissipative stress space:
  // chi_v = pressure - 0.5 * gamma * pc
  // chi_d = eqv_stress
  chi_v = _chi_v_tr - _K * gamma_v * _dt + 0.5 * _gamma * _pcr_tr * (1.0 - std::exp(_L * gamma_v * _dt));
  chi_d = _chi_d_tr - 3.0 * _G * gamma_d * _dt;

  // Update yield parameters
  ADReal pressure = _pressure_tr - _K * gamma_v * _dt;
  _pcr = _pcr_tr * std::exp(_L * gamma_v * _dt);
  _A = (1.0 - _gamma) * pressure + 0.5 * _gamma * _pcr;
  _B = _M * ((1.0 - _alpha) * pressure + 0.5 * _alpha * _gamma * _pcr);
}