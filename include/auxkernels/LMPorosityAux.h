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

  const bool _coupled_pf;
  const VariableValue & _pf_dot;
  const bool _coupled_dam;
  const VariableValue & _damage;
  const VariableValue & _damage_dot;
  const ADMaterialProperty<Real> * _biot;
  const ADMaterialProperty<Real> & _K;
  const ADMaterialProperty<RankTwoTensor> & _strain_incr;
  const bool _has_ve;
  const ADMaterialProperty<RankTwoTensor> * _viscous_strain_incr;
  const bool _has_vp;
  const ADMaterialProperty<RankTwoTensor> * _plastic_strain_incr;
  const ADMaterialProperty<RankTwoTensor> * _stress;
};