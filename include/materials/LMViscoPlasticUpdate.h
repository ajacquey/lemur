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

#include "ADMaterial.h"

class LMViscoPlasticUpdate : public ADMaterial
{
public:
  static InputParameters validParams();
  LMViscoPlasticUpdate(const InputParameters & parameters);
  void setQp(unsigned int qp);
  virtual void viscoPlasticUpdate(ADRankTwoTensor & stress,
                                  const RankFourTensor & Cijkl,
                                  ADRankTwoTensor & elastic_strain_incr) = 0;
  void resetQpProperties() final {}
  void resetProperties() final {}

protected:
  Real _abs_tol;
  Real _rel_tol;
  unsigned int _max_its;
  Real _eta_p;
  Real _n;

  ADMaterialProperty<Real> & _yield_function;
  ADMaterialProperty<RankTwoTensor> & _plastic_strain_incr;
};