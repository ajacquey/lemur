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

#include "LMStrainAux.h"

class LMStrainAux : public LMStrainAuxBase
{
public:
  static InputParameters validParams();
  LMStrainAux(const InputParameters & parameters);

protected:
  virtual Real computeValue() override;

  const unsigned int _i;
  const unsigned int _j;
};