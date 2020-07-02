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

class LMPoroMechMaterial : public LMMechMaterialBase
{
public:
  static InputParameters validParams();
  LMPoroMechMaterial(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual void computeQpElasticityTensor() override;

  // Coupled variables
  const VariableValue & _porosity;

  // Elastic parameters
  const Real _compression_idx;
  const Real _poisson_ratio;

  // Strain quantities
  const MaterialProperty<RankTwoTensor> & _elastic_strain_incr_old;
};