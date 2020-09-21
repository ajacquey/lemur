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

#include "ADMaterial.h"

class LMPoroMaterial : public ADMaterial
{
public:
  static InputParameters validParams();
  LMPoroMaterial(const InputParameters & parameters);

protected:
  virtual void computeQpProperties() override;

  const VariableValue & _porosity;
  const ADVariableValue & _damage;
  const ADVariableValue & _damage_dot;
  const bool _coupled_mech;
  const Real _perm;
  const Real _fluid_visco;
  const Real _Kf;
  const Real _Ks;

  const ADMaterialProperty<Real> * _K;
  const ADMaterialProperty<RankTwoTensor> * _strain_increment;
  const bool _has_ve;
  const ADMaterialProperty<RankTwoTensor> * _viscous_strain_incr;
  const bool _has_vp;
  const ADMaterialProperty<RankTwoTensor> * _plastic_strain_incr;
  const bool _coupled_dam;
  const ADMaterialProperty<RankTwoTensor> * _stress;
  ADMaterialProperty<Real> & _C_biot;
  MaterialProperty<Real> & _fluid_mob;
  ADMaterialProperty<Real> & _biot;
  ADMaterialProperty<Real> & _poro_mech;
};