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

#include "AuxKernel.h"

class LMPorosityAux : public AuxKernel
{
public:
  static InputParameters validParams();
  LMPorosityAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  const VariableValue & _pf_dot;
  const MaterialProperty<Real> & _biot;
  const MaterialProperty<Real> & _K;
  const ADMaterialProperty<RankTwoTensor> & _strain_incr;
  const ADMaterialProperty<RankTwoTensor> & _plastic_strain_incr;
};