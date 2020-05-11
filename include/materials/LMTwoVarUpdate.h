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

#include "LMViscoPlasticUpdate.h"

class LMTwoVarUpdate : public LMViscoPlasticUpdate
{
public:
  static InputParameters validParams();
  LMTwoVarUpdate(const InputParameters & parameters);
  virtual void viscoPlasticUpdate(ADRankTwoTensor & stress,
                                  const RankFourTensor & Cijkl,
                                  ADRankTwoTensor & elastic_strain_incr) override;

protected:
  virtual void returnMap(ADReal & gamma_v, ADReal & gamma_d);
  virtual void
  residual(const ADReal & gamma_v, const ADReal & gamma_d, ADReal & resv, ADReal & resd);
  virtual void jacobian(const ADReal & gamma_v,
                        const ADReal & gamma_d,
                        ADReal & jacvv,
                        ADReal & jacdd,
                        ADReal & jacvd,
                        ADReal & jacdv);
  virtual ADReal yieldFunction(const ADReal & chi_v, const ADReal & chi_d) = 0;
  virtual void
  overStress(const ADReal & gamma_v, const ADReal & gamma_d, ADReal & over_v, ADReal & over_d) = 0;
  virtual void overStressDerivV(const ADReal & gamma_v,
                                const ADReal & gamma_d,
                                ADReal & over_v_v,
                                ADReal & over_d_v) = 0;
  virtual void overStressDerivD(const ADReal & gamma_v,
                                const ADReal & gamma_d,
                                ADReal & over_d_v,
                                ADReal & over_d_d) = 0;
  virtual ADRankTwoTensor reformPlasticStrainTensor(const ADReal & gamma_v,
                                                    const ADReal & gamma_d) = 0;
  virtual void preReturnMap() = 0;
  virtual void postReturnMap(const ADReal & gamma_v, const ADReal & gamma_d) = 0;
  virtual void updateDissipativeStress(const ADReal & gamma_v,
                                       const ADReal & gamma_d,
                                       ADReal & chi_v,
                                       ADReal & chi_d) = 0;

  ADRankTwoTensor _stress_tr;
  ADReal _K;
  ADReal _G;
};