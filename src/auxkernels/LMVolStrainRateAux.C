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

#include "LMVolStrainRateAux.h"

registerMooseObject("LemurApp", LMVolStrainRateAux);

template <>
InputParameters
validParams<LMVolStrainRateAux>()
{
  InputParameters params = validParams<LMStrainAuxBase>();
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
  return (*_strain_incr)[_qp].trace() / _dt;
}
