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

#include "Action.h"

class LMPressureAction;

template <>
InputParameters validParams<LMPressureAction>();

class LMPressureAction : public Action
{
public:
  LMPressureAction(const InputParameters & params);

  virtual void act() override;

protected:
  std::vector<std::vector<AuxVariableName>> _save_in_vars;
  std::vector<bool> _has_save_in_vars;
};