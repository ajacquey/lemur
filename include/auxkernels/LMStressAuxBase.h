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

#pragma once

#include "AuxKernel.h"
#include "RankTwoTensor.h"

class LMStressAuxBase : public AuxKernel
{
public:
  static InputParameters validParams();
  LMStressAuxBase(const InputParameters & parameters);

protected:
  const ADMaterialProperty<RankTwoTensor> & _stress;
};