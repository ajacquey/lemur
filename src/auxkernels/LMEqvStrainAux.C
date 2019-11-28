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

#include "LMEqvStrainAux.h"

registerMooseObject("LemurApp", LMEqvStrainAux);

template <>
InputParameters
validParams<LMEqvStrainAux>()
{
  InputParameters params = validParams<LMStrainAuxBase>();
  params.addClassDescription("Calculates the equivalent strain of the given tensor.");
  return params;
}

LMEqvStrainAux::LMEqvStrainAux(const InputParameters & parameters) : LMStrainAuxBase(parameters) {}

Real
LMEqvStrainAux::computeValue()
{
  return _u_old[_qp] + std::sqrt(2.0 / 3.0) * (*_strain_incr)[_qp].deviatoric().L2norm();
}
