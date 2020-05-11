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

#pragma once

#include "ADKernel.h"

class LMStressDivergence : public ADKernel
{
public:
  static InputParameters validParams();
  LMStressDivergence(const InputParameters & parameters);

protected:
  virtual ADReal computeQpResidual() override;

  const unsigned int _component;
  const Real _rho;
  const RealVectorValue _gravity;
  const ADMaterialProperty<RankTwoTensor> & _stress;
};