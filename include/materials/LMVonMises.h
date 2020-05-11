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

#include "LMSingleVarUpdate.h"

class LMVonMises : public LMSingleVarUpdate
{
public:
  static InputParameters validParams();
  LMVonMises(const InputParameters & parameters);

protected:
  virtual void initQpStatefulProperties() override;
  virtual ADReal yieldFunction(const ADReal & gamma_vp) override;
  virtual ADReal yieldFunctionDeriv(const ADReal & gamma_vp) override;
  virtual void preReturnMap() override;
  virtual void postReturnMap(const ADReal & gamma_vp) override;
  virtual ADRankTwoTensor reformPlasticStrainTensor(const ADReal & gamma_vp) override;

  const bool _coupled_v;
  const ADVariableValue & _v;
  const Real _yield_strength;
  const Real _hg;
  const Real _ht;
  const bool _has_hardening;
  ADMaterialProperty<Real> * _intnl;
  const MaterialProperty<Real> * _intnl_old;

  ADReal _tau_tr;
  ADReal _yield_strength_tr;
};