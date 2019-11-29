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

#include "LMDamageAlphaGammaYield.h"

registerADMooseObject("LemurApp", LMDamageAlphaGammaYield);

defineADValidParams(
    LMDamageAlphaGammaYield,
    LMAlphaGammaYield,
    params.addClassDescription("Viscoplastic update based on the damaged alpha-gamma yield functions.");
    // Coupled variables
    params.addRequiredCoupledVar(
        "damage",
        "The damage variable."););

template <ComputeStage compute_stage>
LMDamageAlphaGammaYield<compute_stage>::LMDamageAlphaGammaYield(const InputParameters & parameters)
  : LMAlphaGammaYield<compute_stage>(parameters),
    // Coupled variables
    _damage(adCoupledValue("damage"))
{
}

template <ComputeStage compute_stage>
void
LMDamageAlphaGammaYield<compute_stage>::preReturnMap()
{
  // Damage correction
  _K *= (1.0 - _damage[_qp]);
  _G *= (1.0 - _damage[_qp]);

  LMAlphaGammaYield<compute_stage>::preReturnMap();
}

// template <ComputeStage compute_stage>
// void
// LMDamageAlphaGammaYield<compute_stage>::postReturnMap(const ADReal & gamma_v, const ADReal & /*gamma_d*/)
// {
//   if (_has_hardening)
//     (*_intnl)[_qp] = (*_intnl_old)[_qp] + gamma_v * _dt;
// }


// template <ComputeStage compute_stage>
// void
// LMDamageAlphaGammaYield<compute_stage>::calculateProjectionDerivV(const ADReal & chi_v,
//                                                      const ADReal & chi_d,
//                                                      ADReal & dchi_v0,
//                                                      ADReal & dchi_d0)
// {
//   // Directions
//   ADReal ev = 0.0, ed = 0.0;
//   ADReal rho_tr = calculateDirection(chi_v, chi_d, ev, ed);
//   // Dissipative stress derivative
//   ADReal dchi_v = -_K * _dt - 0.5 * _gamma * _pcr * _L * _dt;

//   ADReal rho_0 = std::sqrt(1.0 / (Utility::pow<2>(ev / _A) + Utility::pow<2>(ed / _B)));

//   dchi_v0 =
//       rho_0 / rho_tr *
//           (Utility::pow<2>(rho_0) * (1.0 / Utility::pow<2>(_B) - 1.0 / Utility::pow<2>(_A)) *
//                Utility::pow<2>(ev) +
//            1.0) *
//           Utility::pow<2>(ed) * dchi_v +
//       Utility::pow<3>(rho_0) *
//           (Utility::pow<2>(ev) / Utility::pow<3>(_A) + Utility::pow<2>(ed) / Utility::pow<3>(_B)) *
//           ev;
//   dchi_d0 =
//       rho_0 / rho_tr *
//           (Utility::pow<2>(rho_0) * (1.0 / Utility::pow<2>(_B) - 1.0 / Utility::pow<2>(_A)) *
//                Utility::pow<2>(ed) -
//            1.0) *
//           ev * ed * dchi_v +
//       Utility::pow<3>(rho_0) *
//           (Utility::pow<2>(ev) / Utility::pow<3>(_A) + Utility::pow<2>(ed) / Utility::pow<3>(_B)) *
//           ed;
// }

// template <ComputeStage compute_stage>
// void
// LMDamageAlphaGammaYield<compute_stage>::calculateProjectionDerivD(const ADReal & chi_v,
//                                                      const ADReal & chi_d,
//                                                      ADReal & dchi_v0,
//                                                      ADReal & dchi_d0)
// {
//   // Directions
//   ADReal ev = 0.0, ed = 0.0;
//   ADReal rho_tr = calculateDirection(chi_v, chi_d, ev, ed);
//   // Dissipative stress derivative
//   ADReal dchi_d = -3.0 * _G * _dt;

//   ADReal rho_0 = std::sqrt(1.0 / (Utility::pow<2>(ev / _A) + Utility::pow<2>(ed / _B)));

//   dchi_v0 = rho_0 / rho_tr *
//             (Utility::pow<2>(rho_0) * (1.0 / Utility::pow<2>(_A) - 1.0 / Utility::pow<2>(_B)) *
//                  Utility::pow<2>(ev) -
//              1.0) *
//             ev * ed * dchi_d;
//   dchi_d0 = rho_0 / rho_tr *
//             (Utility::pow<2>(rho_0) * (1.0 / Utility::pow<2>(_A) - 1.0 / Utility::pow<2>(_B)) *
//                  Utility::pow<2>(ed) +
//              1.0) *
//             Utility::pow<2>(ev) * dchi_d;
// }

template <ComputeStage compute_stage>
void
LMDamageAlphaGammaYield<compute_stage>::updateYieldParameters(const ADReal & gamma_v)
{
  ADReal pressure = _pressure_tr - _K * gamma_v * _dt;
  _pcr = _pcr_tr * std::exp(_L * gamma_v * _dt);
  _A = (1.0 - _gamma) / (1.0 - _damage[_qp]) * pressure + 0.5 * _gamma * _pcr;
  _B = _M * (pressure - _alpha * std::sqrt(1.0 - _damage[_qp]) * (pressure - 0.5 * _gamma * _pcr));
}

template <ComputeStage compute_stage>
void
LMDamageAlphaGammaYield<compute_stage>::updateYieldParametersDerivV(ADReal & dA, ADReal & dB)
{
  dA = -(1.0 - _gamma) * _K * _dt + 0.5 * _gamma * _L * _dt * _pcr;
  dB = _M * (-_K * _dt + _alpha * std::sqrt(1.0 - _damage[_qp]) * (_K * _dt + 0.5 * _gamma * _L * _dt));
}