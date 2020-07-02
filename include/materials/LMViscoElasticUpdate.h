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

class LMViscoElasticUpdate : public ADMaterial
{
public:
  static InputParameters validParams();
  LMViscoElasticUpdate(const InputParameters & parameters);
  void setQp(unsigned int qp);
  virtual void viscoElasticUpdate(ADRankTwoTensor & stress,
                                  const ADRankFourTensor & Cijkl,
                                  ADRankTwoTensor & elastic_strain_incr);
  void resetQpProperties() final {}
  void resetProperties() final {}

protected:
  virtual ADReal returnMap();
  virtual ADReal residual(const ADReal & gamma_v);
  virtual ADReal jacobian(const ADReal & gamma_v);
  virtual ADReal stressInvariant(const ADReal & gamma_v);
  virtual ADReal stressInvariantDeriv(const ADReal & gamma_v);
  virtual ADRankTwoTensor reformViscousStrainTensor(const ADReal & gamma_v);
  virtual ADReal effectiveViscosity(const ADReal & gamma_v) = 0;
  virtual ADReal creepRate(const ADReal & gamma_v) = 0;
  virtual ADReal creepRateDeriv(const ADReal & gamma_v) = 0;
  virtual void preReturnMap() = 0;
  virtual void postReturnMap(const ADReal & gamma_v) = 0;

  const Real _abs_tol;
  const Real _rel_tol;
  unsigned int _max_its;

  ADRankTwoTensor _stress_tr;
  ADReal _tau_tr;
  ADReal _G;

  ADMaterialProperty<Real> & _viscosity;
  ADMaterialProperty<RankTwoTensor> & _viscous_strain_incr;
};