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

#include "LMEqvStrainRateAux.h"
#include "metaphysicl/raw_type.h"

registerMooseObject("LemurApp", LMEqvStrainRateAux);

InputParameters
LMEqvStrainRateAux::validParams()
{
  InputParameters params = LMStrainAuxBase::validParams();
  params.addClassDescription("Calculates the equivalent strain rate of the given tensor.");
  return params;
}

LMEqvStrainRateAux::LMEqvStrainRateAux(const InputParameters & parameters)
  : LMStrainAuxBase(parameters)
{
}

Real
LMEqvStrainRateAux::computeValue()
{
  return std::sqrt(2.0 / 3.0) * MetaPhysicL::raw_value((*_strain_incr)[_qp].deviatoric().L2norm()) /
         _dt;
}
