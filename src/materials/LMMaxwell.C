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

#include "LMMaxwell.h"

registerMooseObject("LemurApp", LMMaxwell);

InputParameters
LMMaxwell::validParams()
{
  InputParameters params = LMViscoElasticUpdate::validParams();
  params.addClassDescription("Viscoelastic update based on a Maxwell medium.");
  params.addRequiredRangeCheckedParam<Real>(
      "viscosity", "viscosity > 0.0", "The viscosity of the Maxwell medium.");
  return params;
}

LMMaxwell::LMMaxwell(const InputParameters & parameters)
  : LMViscoElasticUpdate(parameters), _eta(getParam<Real>("viscosity"))
{
}

ADReal
LMMaxwell::effectiveViscosity(const ADReal & /*gamma_v*/)
{
  return _eta;
}

ADReal
LMMaxwell::creepRate(const ADReal & gamma_v)
{
  ADReal tau = stressInvariant(gamma_v);
  return tau / (2.0 * _eta);
}

ADReal
LMMaxwell::creepRateDeriv(const ADReal & gamma_v)
{
  ADReal dtau = stressInvariantDeriv(gamma_v);
  return dtau / (2.0 * _eta);
}

void
LMMaxwell::preReturnMap()
{
}

void
LMMaxwell::postReturnMap(const ADReal & /*gamma_v*/)
{
}