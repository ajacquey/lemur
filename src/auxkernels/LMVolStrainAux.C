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

#include "LMVolStrainAux.h"

registerMooseObject("LemurApp", LMVolStrainAux);

template <>
InputParameters
validParams<LMVolStrainAux>()
{
  InputParameters params = validParams<LMStrainAuxBase>();
  params.addClassDescription("Calculates the volumetric strain of the given tensor.");
  return params;
}

LMVolStrainAux::LMVolStrainAux(const InputParameters & parameters) : LMStrainAuxBase(parameters) {}

Real
LMVolStrainAux::computeValue()
{
  return _u_old[_qp] + (*_strain_incr)[_qp].trace();
}
