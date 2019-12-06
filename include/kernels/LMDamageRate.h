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

#include "ADKernel.h"

template <ComputeStage>
class LMDamageRate;

declareADValidParams(LMDamageRate);

template <ComputeStage compute_stage>
class LMDamageRate : public ADKernel<compute_stage>
{
public:
  LMDamageRate(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const ADMaterialProperty(Real) & _damage_rate;

  usingKernelMembers;
};