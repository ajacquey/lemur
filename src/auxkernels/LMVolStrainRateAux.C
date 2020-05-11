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

#include "LMVolStrainRateAux.h"
#include "metaphysicl/raw_type.h"

registerMooseObject("LemurApp", LMVolStrainRateAux);

InputParameters
LMVolStrainRateAux::validParams()
{
  InputParameters params = LMStrainAuxBase::validParams();
  params.addClassDescription("Calculates the volumetric strain rate of the given tensor.");
  return params;
}

LMVolStrainRateAux::LMVolStrainRateAux(const InputParameters & parameters)
  : LMStrainAuxBase(parameters)
{
}

Real
LMVolStrainRateAux::computeValue()
{
  return MetaPhysicL::raw_value((*_strain_incr)[_qp].trace()) / _dt;
}
