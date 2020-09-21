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

#include "LMDamageAlphaGammaYield.h"

registerMooseObject("LemurApp", LMDamageAlphaGammaYield);

InputParameters
LMDamageAlphaGammaYield::validParams()
{
  InputParameters params = LMAlphaGammaYield::validParams();
  params.addClassDescription(
      "Viscoplastic update based on the damaged alpha-gamma yield functions.");
  // Coupled variables
  params.addRequiredCoupledVar("damage", "The damage variable.");
  // // Damage viscosity
  // params.addRequiredRangeCheckedParam<Real>(
  //     "damage_viscosity", "damage_viscosity > 0.0", "The damage viscosity.");
  // Dissipation contributions
  params.addRequiredRangeCheckedParam<Real>(
      "rv", "rv > 0.0", "The volumetric dissipation contribution.");
  params.addRequiredRangeCheckedParam<Real>(
      "rs", "rs > 0.0", "The shear dissipation contribution.");
  return params;
}

LMDamageAlphaGammaYield::LMDamageAlphaGammaYield(const InputParameters & parameters)
  : LMAlphaGammaYield(parameters),
    // Coupled variables
    _damage(adCoupledValue("damage")),
    // // Damage viscosity
    // _eta_a(getParam<Real>("damage_viscosity")),
    // Dissipation contributions
    _rv(getParam<Real>("rv")),
    _rs(getParam<Real>("rs")),
    // Properties
    _damage_rate(declareADProperty<Real>("damage_rate"))
{
}

void
LMDamageAlphaGammaYield::preReturnMap()
{
  // Damage correction
  _K *= (1.0 - _damage[_qp]);
  _G *= (1.0 - _damage[_qp]);

  // Damage driving
  _damage_rate[_qp] = 0.0;

  LMAlphaGammaYield::preReturnMap();
}

void
LMDamageAlphaGammaYield::overStress(const ADReal & gamma_v,
                                    const ADReal & gamma_d,
                                    ADReal & over_v,
                                    ADReal & over_d)
{
  // Dissipative stresses
  ADReal chi_v = 0.0, chi_d = 0.0;
  updateDissipativeStress(gamma_v, gamma_d, chi_v, chi_d);

  ADReal f = yieldFunction(chi_v, chi_d);
  ADReal df_dchi_v = dyieldFunctiondVol(chi_v, chi_d);
  ADReal df_dchi_d = dyieldFunctiondDev(chi_v, chi_d);

  over_v = Utility::pow<2>(_rv) * std::pow(f, _n) * df_dchi_v;
  over_d = Utility::pow<2>(_rs) * std::pow(f, _n) * df_dchi_d;
}

void
LMDamageAlphaGammaYield::overStressDerivV(const ADReal & gamma_v,
                                          const ADReal & gamma_d,
                                          ADReal & over_v_v,
                                          ADReal & over_d_v)
{
  // Dissipative stresses
  ADReal chi_v = 0.0, chi_d = 0.0;
  updateDissipativeStress(gamma_v, gamma_d, chi_v, chi_d);

  // Dissipative stress derivatives
  ADReal dchi_v = -_K * _dt - 0.5 * _gamma * _pcr * _L * _dt;

  // Yield and derivatives
  ADReal f = yieldFunction(chi_v, chi_d);
  ADReal df_dchi_v = dyieldFunctiondVol(chi_v, chi_d);
  ADReal df_dchi_d = dyieldFunctiondDev(chi_v, chi_d);
  ADReal d2f_dchi_v2 = d2yieldFunctiondVol2(chi_v, chi_d);
  ADReal d2f_dchi_d_dchi_v = d2yieldFunctiondDevVol(chi_v, chi_d);
  ADReal df_dA = dyieldFunctiondA(chi_v, chi_d);
  ADReal df_dB = dyieldFunctiondA(chi_v, chi_d);
  ADReal d2f_dchi_v_dA = d2yieldFunctiondVolA(chi_v, chi_d);
  ADReal d2f_dchi_v_dB = d2yieldFunctiondVolB(chi_v, chi_d);
  ADReal d2f_dchi_d_dA = d2yieldFunctiondDevA(chi_v, chi_d);
  ADReal d2f_dchi_d_dB = d2yieldFunctiondDevA(chi_v, chi_d);

  // Yield parameters derivatives
  ADReal dA = 0.0, dB = 0.0;
  updateYieldParametersDerivV(dA, dB);

  // Over stress derivatives wrt dissipative stress
  ADReal over_v_dchi_v = Utility::pow<2>(_rv) * std::pow(f, _n - 1.0) *
                         (_n * Utility::pow<2>(df_dchi_v) + f * d2f_dchi_v2);
  ADReal over_d_dchi_v = Utility::pow<2>(_rs) * std::pow(f, _n - 1.0) *
                         (_n * df_dchi_v * df_dchi_d + f * d2f_dchi_d_dchi_v);

  // Over stress derivatives wrt yield parameters
  ADReal over_v_dA =
      Utility::pow<2>(_rv) * std::pow(f, _n - 1.0) * (_n * df_dA * df_dchi_v + f * d2f_dchi_v_dA);
  ADReal over_v_dB =
      Utility::pow<2>(_rv) * std::pow(f, _n - 1.0) * (_n * df_dB * df_dchi_v + f * d2f_dchi_v_dB);
  ADReal over_d_dA =
      Utility::pow<2>(_rs) * std::pow(f, _n - 1.0) * (_n * df_dA * df_dchi_d + f * d2f_dchi_d_dA);
  ADReal over_d_dB =
      Utility::pow<2>(_rs) * std::pow(f, _n - 1.0) * (_n * df_dB * df_dchi_d + f * d2f_dchi_d_dB);

  over_v_v = over_v_dchi_v * dchi_v + over_v_dA * dA + over_v_dB * dB;
  over_d_v = over_d_dchi_v * dchi_v + over_d_dA * dA + over_d_dB * dB;
}

void
LMDamageAlphaGammaYield::overStressDerivD(const ADReal & gamma_v,
                                          const ADReal & gamma_d,
                                          ADReal & over_v_d,
                                          ADReal & over_d_d)
{
  // Dissipative stresses
  ADReal chi_v = 0.0, chi_d = 0.0;
  updateDissipativeStress(gamma_v, gamma_d, chi_v, chi_d);

  // Dissipative stress derivatives
  ADReal dchi_d = -3.0 * _G * _dt;

  // Yield and derivatives
  ADReal f = yieldFunction(chi_v, chi_d);
  ADReal df_dchi_v = dyieldFunctiondVol(chi_v, chi_d);
  ADReal df_dchi_d = dyieldFunctiondDev(chi_v, chi_d);
  ADReal d2f_dchi_v_dchi_d = d2yieldFunctiondVolDev(chi_v, chi_d);
  ADReal d2f_dchi_d2 = d2yieldFunctiondDev2(chi_v, chi_d);

  // Over stress derivatives wrt dissipative stress
  ADReal over_v_dchi_d = Utility::pow<2>(_rv) * std::pow(f, _n - 1.0) *
                         (_n * df_dchi_v * df_dchi_d + f * d2f_dchi_v_dchi_d);
  ADReal over_d_dchi_d = Utility::pow<2>(_rs) * std::pow(f, _n - 1.0) *
                         (_n * Utility::pow<2>(df_dchi_d) + f * d2f_dchi_d2);

  over_v_d = over_v_dchi_d * dchi_d;
  over_d_d = over_d_dchi_d * dchi_d;
}

void
LMDamageAlphaGammaYield::updateYieldParameters(const ADReal & gamma_v)
{
  ADReal pressure = _pressure_tr - _K * gamma_v * _dt;
  _pcr = _pcr_tr * std::exp(_L * gamma_v * _dt);
  _one_on_A = (1.0 - _damage[_qp]) /
              ((1.0 - _gamma) * pressure + 0.5 * (1.0 - _damage[_qp]) * _gamma * _pcr);
  _one_on_B =
      1.0 /
      (_M * (pressure - _alpha * std::sqrt(1.0 - _damage[_qp]) * (pressure - 0.5 * _gamma * _pcr)));
}

void
LMDamageAlphaGammaYield::updateYieldParametersDerivV(ADReal & dA, ADReal & dB)
{
  dA = (_damage[_qp] != 1.0)
           ? Utility::pow<2>(_one_on_A) *
                 ((1.0 - _gamma) / (1.0 - _damage[_qp]) * _K * _dt - 0.5 * _gamma * _L * _dt * _pcr)
           : 0.0;
  dB = Utility::pow<2>(_one_on_B) * _M *
       ((1.0 - _alpha * std::sqrt(1.0 - _damage[_qp])) * _K * _dt -
        0.5 * std::sqrt(1.0 - _damage[_qp]) * _alpha * _gamma * _L * _dt * _pcr);
}

void
LMDamageAlphaGammaYield::postReturnMap(const ADReal & gamma_v, const ADReal & gamma_d)
{
  LMAlphaGammaYield::postReturnMap(gamma_v, gamma_d);

  ADReal pressure = _pressure_tr - _K * gamma_v * _dt;
  ADReal eqv_stress = _eqv_stress_tr - 3.0 * _G * gamma_d * _dt;
  ADReal chi_v = 0.0, chi_d = 0.0;
  updateDissipativeStress(gamma_v, gamma_d, chi_v, chi_d);
  // Damage driving force
  ADReal Ya = 0.5 / (1.0 - _damage[_qp]) *
              (Utility::pow<2>(pressure) / _K + Utility::pow<2>(eqv_stress) / (3.0 * _G));
  // _damage_rate[_qp] = _yield_function[_qp] / (_eta_a * Ya);
  _damage_rate[_qp] = chi_v / Ya * (1.0 - Utility::pow<2>(_rv)) / Utility::pow<2>(_rv) * gamma_v +
                      chi_d / Ya * (1.0 - Utility::pow<2>(_rs)) / Utility::pow<2>(_rs) * gamma_d;
}