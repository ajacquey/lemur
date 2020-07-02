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

class LMViscoPlasticUpdate : public ADMaterial
{
public:
  static InputParameters validParams();
  LMViscoPlasticUpdate(const InputParameters & parameters);
  void setQp(unsigned int qp);
  virtual void viscoPlasticUpdate(ADRankTwoTensor & stress,
                                  const ADRankFourTensor & Cijkl,
                                  ADRankTwoTensor & elastic_strain_incr) = 0;
  void resetQpProperties() final {}
  void resetProperties() final {}

protected:
  const ADVariableValue & _pf;
  const Real _abs_tol;
  const Real _rel_tol;
  const unsigned int _max_its;
  ADReal _eta_p;
  const Real _n;
  const Real _pf0;
  const Real _Ar;

  ADMaterialProperty<Real> & _yield_function;
  ADMaterialProperty<RankTwoTensor> & _plastic_strain_incr;
};