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

#include "LMStressAuxBase.h"

InputParameters
LMStressAuxBase::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Base class for outputting stress values.");
  return params;
}

LMStressAuxBase::LMStressAuxBase(const InputParameters & parameters)
  : AuxKernel(parameters), _stress(getADMaterialProperty<RankTwoTensor>("stress"))
{
}