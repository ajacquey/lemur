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

#include "LMVolStrainAux.h"
#include "metaphysicl/raw_type.h"

registerMooseObject("LemurApp", LMVolStrainAux);

InputParameters
LMVolStrainAux::validParams()
{
  InputParameters params = LMStrainAuxBase::validParams();
  params.addClassDescription("Calculates the volumetric strain of the given tensor.");
  return params;
}

LMVolStrainAux::LMVolStrainAux(const InputParameters & parameters) : LMStrainAuxBase(parameters) {}

Real
LMVolStrainAux::computeValue()
{
  return _u_old[_qp] + MetaPhysicL::raw_value((*_strain_incr)[_qp].trace());
}
