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

#include "LMAlphaGammaYield.h"

registerMooseObject("LemurApp", LMAlphaGammaYield);

InputParameters
LMAlphaGammaYield::validParams()
{
  InputParameters params = LMTwoVarUpdate::validParams();
  params.addClassDescription("Viscoplastic update based on the alpha-gamma yield functions.");
  params.addRequiredRangeCheckedParam<Real>(
      "friction_angle", "friction_angle > 0.0", "The friction angle for the critical state line.");
  params.addRequiredRangeCheckedParam<Real>(
      "critical_pressure", "critical_pressure > 0.0", "The critical pressure of the capped yield.");
  params.addRangeCheckedParam<Real>(
      "alpha", 1.0, "alpha >= 0.0 & alpha <= 1.0", "The alpha parameter for the yield/");
  params.addRangeCheckedParam<Real>(
      "gamma", 1.0, "gamma >= 0.0 & gamma <= 1.0", "The gamma parameter for the yield/");
  params.addRangeCheckedParam<Real>("critical_pressure_hardening",
                                    0.0,
                                    "critical_pressure_hardening >= 0.0",
                                    "The hardening parameter in the exponential for the critical "
                                    "pressure of the capped yield.");
  params.addRangeCheckedParam<Real>("porosity_hardening",
                                    0.0,
                                    "porosity_hardening>=0.0",
                                    "The coefficient for porosity hardening of the viscosity.");
  return params;
}

LMAlphaGammaYield::LMAlphaGammaYield(const InputParameters & parameters)
  : LMTwoVarUpdate(parameters),
    _phi(getParam<Real>("friction_angle")),
    _pcr0(getParam<Real>("critical_pressure")),
    _alpha(getParam<Real>("alpha")),
    _gamma(getParam<Real>("gamma")),
    _L(getParam<Real>("critical_pressure_hardening")),
    _has_hardening(_L != 0.0),
    _intnl(_has_hardening ? &declareADProperty<Real>("volumetric_plastic_strain") : nullptr),
    _intnl_old(_has_hardening ? &getMaterialPropertyOld<Real>("volumetric_plastic_strain")
                              : nullptr)
{
  _M = std::sqrt(3.0) * std::sin(_phi * libMesh::pi / 180.0);
}

void
LMAlphaGammaYield::initQpStatefulProperties()
{
  if (_has_hardening)
    (*_intnl)[_qp] = 0.0;
}

ADReal
LMAlphaGammaYield::yieldFunction(const ADReal & chi_v, const ADReal & chi_d)
{
  return std::sqrt(Utility::pow<2>(chi_v * _one_on_A) + Utility::pow<2>(chi_d * _one_on_B)) - 1.0;
}

void
LMAlphaGammaYield::overStress(const ADReal & gamma_v,
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

  over_v = std::pow(f, _n) * df_dchi_v;
  over_d = std::pow(f, _n) * df_dchi_d;
}

void
LMAlphaGammaYield::overStressDerivV(const ADReal & gamma_v,
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
  ADReal over_v_dchi_v =
      std::pow(f, _n - 1.0) * (_n * Utility::pow<2>(df_dchi_v) + f * d2f_dchi_v2);
  ADReal over_d_dchi_v =
      std::pow(f, _n - 1.0) * (_n * df_dchi_v * df_dchi_d + f * d2f_dchi_d_dchi_v);

  // Over stress derivatives wrt yield parameters
  ADReal over_v_dA = std::pow(f, _n - 1.0) * (_n * df_dA * df_dchi_v + f * d2f_dchi_v_dA);
  ADReal over_v_dB = std::pow(f, _n - 1.0) * (_n * df_dB * df_dchi_v + f * d2f_dchi_v_dB);
  ADReal over_d_dA = std::pow(f, _n - 1.0) * (_n * df_dA * df_dchi_d + f * d2f_dchi_d_dA);
  ADReal over_d_dB = std::pow(f, _n - 1.0) * (_n * df_dB * df_dchi_d + f * d2f_dchi_d_dB);

  over_v_v = over_v_dchi_v * dchi_v + over_v_dA * dA + over_v_dB * dB;
  over_d_v = over_d_dchi_v * dchi_v + over_d_dA * dA + over_d_dB * dB;
}

void
LMAlphaGammaYield::overStressDerivD(const ADReal & gamma_v,
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
  ADReal over_v_dchi_d =
      std::pow(f, _n - 1.0) * (_n * df_dchi_v * df_dchi_d + f * d2f_dchi_v_dchi_d);
  ADReal over_d_dchi_d =
      std::pow(f, _n - 1.0) * (_n * Utility::pow<2>(df_dchi_d) + f * d2f_dchi_d2);

  over_v_d = over_v_dchi_d * dchi_d;
  over_d_d = over_d_dchi_d * dchi_d;
}

void
LMAlphaGammaYield::preReturnMap()
{
  _pressure_tr = -_stress_tr.trace() / 3.0;
  _eqv_stress_tr = std::sqrt(1.5) * _stress_tr.deviatoric().L2norm();

  _pcr_tr = _pcr0;
  if (_has_hardening)
  {
    (*_intnl)[_qp] = (*_intnl_old)[_qp];
    _pcr_tr = _pcr0 * std::exp(_L * (*_intnl_old)[_qp]);
  }

  _chi_v_tr = _pressure_tr - 0.5 * _gamma * _pcr_tr;
  _chi_d_tr = _eqv_stress_tr;

  _eta_p = getParam<Real>("plastic_viscosity");
}

void
LMAlphaGammaYield::postReturnMap(const ADReal & gamma_v, const ADReal & /*gamma_d*/)
{
  if (_has_hardening)
    (*_intnl)[_qp] = (*_intnl_old)[_qp] + gamma_v * _dt;
}

ADRankTwoTensor
LMAlphaGammaYield::reformPlasticStrainTensor(const ADReal & gamma_v, const ADReal & gamma_d)
{
  ADRankTwoTensor flow_dir =
      (_eqv_stress_tr != 0.0) ? _stress_tr.deviatoric() / _eqv_stress_tr : ADRankTwoTensor();

  ADRankTwoTensor delta_gamma = 1.5 * gamma_d * _dt * flow_dir;
  delta_gamma.addIa(-gamma_v * _dt / 3.0);

  return delta_gamma;
}

void
LMAlphaGammaYield::updateYieldParametersDerivV(ADReal & dA, ADReal & dB)
{
  dA = Utility::pow<2>(_one_on_A) * ((1.0 - _gamma) * _K * _dt - 0.5 * _gamma * _L * _dt * _pcr);
  dB = Utility::pow<2>(_one_on_B) * _M *
       ((1.0 - _alpha) * _K * _dt - 0.5 * _alpha * _gamma * _L * _dt * _pcr);
}

void
LMAlphaGammaYield::updateDissipativeStress(const ADReal & gamma_v,
                                           const ADReal & gamma_d,
                                           ADReal & chi_v,
                                           ADReal & chi_d)
{
  // Here we calculate the yield function in the dissipative stress space:
  // chi_v = pressure - 0.5 * gamma * pc
  // chi_d = eqv_stress
  chi_v = _chi_v_tr - _K * gamma_v * _dt +
          0.5 * _gamma * _pcr_tr * (1.0 - std::exp(_L * gamma_v * _dt));
  chi_d = _chi_d_tr - 3.0 * _G * gamma_d * _dt;

  // Update yield parameters
  updateYieldParameters(gamma_v);
}

void
LMAlphaGammaYield::updateYieldParameters(const ADReal & gamma_v)
{
  ADReal pressure = _pressure_tr - _K * gamma_v * _dt;
  _pcr = _pcr_tr * std::exp(_L * gamma_v * _dt);
  _one_on_A = 1.0 / ((1.0 - _gamma) * pressure + 0.5 * _gamma * _pcr);
  _one_on_B = 1.0 / (_M * ((1.0 - _alpha) * pressure + 0.5 * _alpha * _gamma * _pcr));
}

ADReal
LMAlphaGammaYield::dyieldFunctiondVol(const ADReal & chi_v, const ADReal & chi_d)
{
  ADReal f = yieldFunction(chi_v, chi_d);
  return Utility::pow<2>(_one_on_A) * chi_v / (1.0 + f);
}

ADReal
LMAlphaGammaYield::dyieldFunctiondDev(const ADReal & chi_v, const ADReal & chi_d)
{
  ADReal f = yieldFunction(chi_v, chi_d);
  return Utility::pow<2>(_one_on_B) * chi_d / (1.0 + f);
}

ADReal
LMAlphaGammaYield::d2yieldFunctiondVol2(const ADReal & chi_v, const ADReal & chi_d)
{
  ADReal f = yieldFunction(chi_v, chi_d);
  ADReal df_dchi_v = dyieldFunctiondVol(chi_v, chi_d);
  return (Utility::pow<2>(_one_on_A) - Utility::pow<2>(df_dchi_v)) / (1.0 + f);
}

ADReal
LMAlphaGammaYield::d2yieldFunctiondVolDev(const ADReal & chi_v, const ADReal & chi_d)
{
  ADReal f = yieldFunction(chi_v, chi_d);
  ADReal df_dchi_v = dyieldFunctiondVol(chi_v, chi_d);
  ADReal df_dchi_d = dyieldFunctiondDev(chi_v, chi_d);
  return -df_dchi_v * df_dchi_d / (1.0 + f);
}

ADReal
LMAlphaGammaYield::d2yieldFunctiondDevVol(const ADReal & chi_v, const ADReal & chi_d)
{
  ADReal f = yieldFunction(chi_v, chi_d);
  ADReal df_dchi_v = dyieldFunctiondVol(chi_v, chi_d);
  ADReal df_dchi_d = dyieldFunctiondDev(chi_v, chi_d);
  return -df_dchi_v * df_dchi_d / (1.0 + f);
}

ADReal
LMAlphaGammaYield::d2yieldFunctiondDev2(const ADReal & chi_v, const ADReal & chi_d)
{
  ADReal f = yieldFunction(chi_v, chi_d);
  ADReal df_dchi_d = dyieldFunctiondDev(chi_v, chi_d);
  return (Utility::pow<2>(_one_on_B) - Utility::pow<2>(df_dchi_d)) / (1.0 + f);
}

ADReal
LMAlphaGammaYield::dyieldFunctiondA(const ADReal & chi_v, const ADReal & chi_d)
{
  ADReal f = yieldFunction(chi_v, chi_d);
  return Utility::pow<2>(chi_v) * _one_on_A / (1.0 + f);
}

ADReal
LMAlphaGammaYield::dyieldFunctiondB(const ADReal & chi_v, const ADReal & chi_d)
{
  ADReal f = yieldFunction(chi_v, chi_d);
  return Utility::pow<2>(chi_d) * _one_on_B / (1.0 + f);
}

ADReal
LMAlphaGammaYield::d2yieldFunctiondVolA(const ADReal & chi_v, const ADReal & chi_d)
{
  ADReal f = yieldFunction(chi_v, chi_d);
  ADReal df_dA = dyieldFunctiondA(chi_v, chi_d);
  return _one_on_A * chi_v / (1.0 + f) * (2.0 - _one_on_A / (1.0 + f) * df_dA);
}

ADReal
LMAlphaGammaYield::d2yieldFunctiondVolB(const ADReal & chi_v, const ADReal & chi_d)
{
  ADReal f = yieldFunction(chi_v, chi_d);
  ADReal df_dchi_v = dyieldFunctiondVol(chi_v, chi_d);
  ADReal df_dB = dyieldFunctiondB(chi_v, chi_d);
  return -df_dchi_v * df_dB / (1.0 + f);
}

ADReal
LMAlphaGammaYield::d2yieldFunctiondDevA(const ADReal & chi_v, const ADReal & chi_d)
{
  ADReal f = yieldFunction(chi_v, chi_d);
  ADReal df_dchi_d = dyieldFunctiondDev(chi_v, chi_d);
  ADReal df_dA = dyieldFunctiondA(chi_v, chi_d);
  return -df_dchi_d * df_dA / (1.0 + f);
}

ADReal
LMAlphaGammaYield::d2yieldFunctiondDevB(const ADReal & chi_v, const ADReal & chi_d)
{
  ADReal f = yieldFunction(chi_v, chi_d);
  ADReal df_dB = dyieldFunctiondB(chi_v, chi_d);
  return _one_on_B * chi_d / (1.0 + f) * (2.0 - _one_on_B / (1.0 + f) * df_dB);
}