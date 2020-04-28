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

#include "LMViscoPlasticUpdate.h"

class LMSingleVarUpdate : public LMViscoPlasticUpdate
{
public:
  static InputParameters validParams();
  LMSingleVarUpdate(const InputParameters & parameters);
  virtual void viscoPlasticUpdate(ADRankTwoTensor & stress,
                                  const RankFourTensor & Cijkl,
                                  ADRankTwoTensor & elastic_strain_incr) override;

protected:
  virtual ADReal returnMap();
  virtual ADReal residual(const ADReal & gamma_vp);
  virtual ADReal jacobian(const ADReal & gamma_vp);
  virtual ADReal yieldFunction(const ADReal & gamma_vp) = 0;
  virtual ADReal yieldFunctionDeriv(const ADReal & gamma_vp) = 0;
  virtual ADRankTwoTensor reformPlasticStrainTensor(const ADReal & gamma_vp) = 0;
  virtual void preReturnMap() = 0;
  virtual void postReturnMap(const ADReal & gamma_vp) = 0;

  ADRankTwoTensor _stress_tr;
  ADReal _K;
  ADReal _G;
};