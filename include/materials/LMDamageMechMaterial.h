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

#include "LMMechMaterial.h"

class LMViscoPlasticUpdate;

class LMDamageMechMaterial : public LMMechMaterial
{
public:
  static InputParameters validParams();
  LMDamageMechMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpElasticGuess() override;

  // Coupled variables
  const ADVariableValue & _damage_dot;
  const VariableValue & _damage_old;
};