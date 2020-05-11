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

#include "LMStressAux.h"
#include "metaphysicl/raw_type.h"

registerMooseObject("LemurApp", LMStressAux);

InputParameters
LMStressAux::validParams()
{
  InputParameters params = LMStressAuxBase::validParams();
  params.addClassDescription("Class for outputting a component of the stress tensor based on the "
                             "indexes i and j. Returns stress(i, j).");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "index_i", "index_i>=0 & index_i<=2", "The i index.");
  params.addRequiredRangeCheckedParam<unsigned int>(
      "index_j", "index_j>=0 & index_j<=2", "The j index.");
  return params;
}

LMStressAux::LMStressAux(const InputParameters & parameters)
  : LMStressAuxBase(parameters),
    _i(getParam<unsigned int>("index_i")),
    _j(getParam<unsigned int>("index_j"))
{
}

Real
LMStressAux::computeValue()
{
  return MetaPhysicL::raw_value(_stress[_qp](_i, _j));
}