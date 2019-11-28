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

registerMooseObject("LemurApp", LMVonMisesStressAux);

template <>
InputParameters
validParams<LMVonMisesStressAux>()
{
  InputParameters params = validParams<LMStressAuxBase>();
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
  RankTwoTensor stress_dev = _stress[_qp].deviatoric();
  return std::sqrt(1.5) * stress_dev.L2norm();
}
