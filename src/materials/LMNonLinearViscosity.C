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

#include "LMNonLinearViscosity.h"

registerMooseObject("LemurApp", LMNonLinearViscosity);

InputParameters
LMNonLinearViscosity::validParams()
{
  InputParameters params = LMViscoElasticUpdate::validParams();
  params.addClassDescription("Viscoelastic update based on a Maxwell medium.");
  params.addRequiredRangeCheckedParam<Real>(
      "viscosity", "viscosity > 0.0", "The viscosity coefficient for the non-linear medium.");
  params.addRequiredRangeCheckedParam<Real>(
      "exponent", "exponent>1.0 & exponent<=2.0", "The exponent for the non-linear viscosity.");
  return params;
}

LMNonLinearViscosity::LMNonLinearViscosity(const InputParameters & parameters)
  : LMViscoElasticUpdate(parameters),
    _eta(getParam<Real>("viscosity")),
    _n(getParam<Real>("exponent"))
{
}

ADReal
LMNonLinearViscosity::effectiveViscosity(const ADReal & gamma_v)
{
  ADReal tau = stressInvariant(gamma_v);

  return _eta * std::pow(tau, (_n - 2.0) / (_n - 1.0));
}

ADReal
LMNonLinearViscosity::creepRate(const ADReal & gamma_v)
{
  ADReal tau = stressInvariant(gamma_v);
  return std::pow(tau, 1.0 / (_n - 1.0)) / (2.0 * _eta);
}

ADReal
LMNonLinearViscosity::creepRateDeriv(const ADReal & gamma_v)
{
  ADReal tau = stressInvariant(gamma_v);
  ADReal dtau = stressInvariantDeriv(gamma_v);
  return 1.0 / (_n - 1.0) * dtau * std::pow(tau, (2.0 - _n) / (_n - 1.0)) / (2.0 * _eta);
}

void
LMNonLinearViscosity::preReturnMap()
{
  _tau_tr = std::sqrt(0.5) * _stress_tr.deviatoric().L2norm();
}

void
LMNonLinearViscosity::postReturnMap(const ADReal & /*gamma_v*/)
{
}