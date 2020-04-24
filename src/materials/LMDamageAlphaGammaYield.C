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

registerMooseObject("LemurApp", LMDamageAlphaGammaYield);

InputParameters
LMDamageAlphaGammaYield::validParams()
{
  InputParameters params = LMAlphaGammaYield::validParams();
  params.addClassDescription(
      "Viscoplastic update based on the damaged alpha-gamma yield functions.");
  // Coupled variables
  params.addRequiredCoupledVar("damage", "The damage variable.");
  // Damage viscosity
  params.addRequiredRangeCheckedParam<Real>(
      "damage_viscosity", "damage_viscosity > 0.0", "The damage viscosity.");
  return params;
}

LMDamageAlphaGammaYield::LMDamageAlphaGammaYield(const InputParameters & parameters)
  : LMAlphaGammaYield(parameters),
    // Coupled variables
    _damage(adCoupledValue("damage")),
    // Damage viscosity
    _eta_a(getParam<Real>("damage_viscosity")),
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
  // Damage driving force
  ADReal Ya = 0.5 / (1.0 - _damage[_qp]) *
              (Utility::pow<2>(pressure) / _K + Utility::pow<2>(eqv_stress) / (3.0 * _G));
  _damage_rate[_qp] = _yield_function[_qp] / (_eta_a * Ya);
}