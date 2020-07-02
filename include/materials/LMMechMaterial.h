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

#include "LMMechMaterialBase.h"

class LMMechMaterial : public LMMechMaterialBase
{
public:
  static InputParameters validParams();
  LMMechMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpElasticityTensor() override;

  // Elastic parameters
  const Real _bulk_modulus;
  const Real _shear_modulus;
};