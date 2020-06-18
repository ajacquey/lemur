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

#include "LMViscoElasticUpdate.h"

class LMNonLinearViscosity : public LMViscoElasticUpdate
{
public:
  static InputParameters validParams();
  LMNonLinearViscosity(const InputParameters & parameters);

protected:
  virtual ADReal effectiveViscosity(const ADReal & gamma_v) override;
  virtual ADReal creepRate(const ADReal & gamma_v) override;
  virtual ADReal creepRateDeriv(const ADReal & gamma_v) override;
  virtual void preReturnMap() override;
  virtual void postReturnMap(const ADReal & gamma_v) override;

  const Real _eta;
  const Real _n;
};