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

#include "LMPressureAux.h"

registerMooseObject("LemurApp", LMPressureAux);

template <>
InputParameters
validParams<LMPressureAux>()
{
  InputParameters params = validParams<LMStressAuxBase>();
  params.addClassDescription("Calculates the pressure.");
  return params;
}

LMPressureAux::LMPressureAux(const InputParameters & parameters) : LMStressAuxBase(parameters) {}

Real
LMPressureAux::computeValue()
{
  return -_stress[_qp].trace() / 3.0;
}