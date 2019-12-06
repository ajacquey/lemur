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

#include "IntegratedBC.h"

class LMPressureBC;
class Function;

template <>
InputParameters validParams<LMPressureBC>();

class LMPressureBC : public IntegratedBC
{
public:
  LMPressureBC(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();

  const int _component;
  const Real _value;
  const Function * _function;
};