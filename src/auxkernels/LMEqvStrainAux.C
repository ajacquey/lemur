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

#include "LMEqvStrainAux.h"
#include "metaphysicl/raw_type.h"

registerMooseObject("LemurApp", LMEqvStrainAux);

InputParameters
LMEqvStrainAux::validParams()
{
  InputParameters params = LMStrainAuxBase::validParams();
  params.addClassDescription("Calculates the equivalent strain of the given tensor.");
  return params;
}

LMEqvStrainAux::LMEqvStrainAux(const InputParameters & parameters) : LMStrainAuxBase(parameters) {}

Real
LMEqvStrainAux::computeValue()
{
  return _u_old[_qp] +
         std::sqrt(2.0 / 3.0) * MetaPhysicL::raw_value((*_strain_incr)[_qp].deviatoric().L2norm());
}
