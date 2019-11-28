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

#include "AuxKernel.h"
#include "RankTwoTensor.h"
#include "DerivativeMaterialInterface.h"

class LMStrainAuxBase;

template <>
InputParameters validParams<LMStrainAuxBase>();

class LMStrainAuxBase : public DerivativeMaterialInterface<AuxKernel>
{
public:
  LMStrainAuxBase(const InputParameters & parameters);
  static MooseEnum strainType();

protected:
  MooseEnum _strain_type;
  std::string _strain_name;
  const MaterialProperty<RankTwoTensor> * _strain_incr;
};