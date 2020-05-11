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

#include "LMPressureAux.h"
#include "metaphysicl/raw_type.h"

registerMooseObject("LemurApp", LMPressureAux);

InputParameters
LMPressureAux::validParams()
{
  InputParameters params = LMStressAuxBase::validParams();
  params.addClassDescription("Calculates the pressure.");
  return params;
}

LMPressureAux::LMPressureAux(const InputParameters & parameters) : LMStressAuxBase(parameters) {}

Real
LMPressureAux::computeValue()
{
  return -MetaPhysicL::raw_value(_stress[_qp].trace()) / 3.0;
}