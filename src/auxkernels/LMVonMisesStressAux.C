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

#include "LMVonMisesStressAux.h"
#include "metaphysicl/raw_type.h"

registerMooseObject("LemurApp", LMVonMisesStressAux);

InputParameters
LMVonMisesStressAux::validParams()
{
  InputParameters params = LMStressAuxBase::validParams();
  params.addClassDescription("Calculates the Von Mises stress.");
  return params;
}

LMVonMisesStressAux::LMVonMisesStressAux(const InputParameters & parameters)
  : LMStressAuxBase(parameters)
{
}

Real
LMVonMisesStressAux::computeValue()
{
  ADRankTwoTensor stress_dev = _stress[_qp].deviatoric();
  return std::sqrt(1.5) * MetaPhysicL::raw_value(stress_dev.L2norm());
}
