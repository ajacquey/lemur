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

#include "LMStrainAux.h"
#include "metaphysicl/raw_type.h"

registerMooseObject("LemurApp", LMStrainAux);

InputParameters
LMStrainAux::validParams()
{
  InputParameters params = LMStrainAuxBase::validParams();
  params.addClassDescription("Class for outputting a component of the strain tensor based on the "
                             "indexes i and j. Returns stress(i, j).");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "index_i", "index_i>=0 & index_i<=2", "The i index.");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "index_j", "index_j>=0 & index_j<=2", "The j index.");
  return params;
}

LMStrainAux::LMStrainAux(const InputParameters & parameters)
  : LMStrainAuxBase(parameters),
    _i(getParam<unsigned int>("index_i")),
    _j(getParam<unsigned int>("index_j"))
{
}

Real
LMStrainAux::computeValue()
{
  return _u_old[_qp] + MetaPhysicL::raw_value((*_strain_incr)[_qp](_i, _j));
}