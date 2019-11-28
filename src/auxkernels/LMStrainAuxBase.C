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

#include "LMStrainAuxBase.h"

template <>
InputParameters
validParams<LMStrainAuxBase>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("Base class for outputting strain values.");
  params.addParam<MooseEnum>("strain_type",
                             LMStrainAuxBase::strainType() = "total",
                             "The type of the strain tensor to output.");
  return params;
}

LMStrainAuxBase::LMStrainAuxBase(const InputParameters & parameters)
  : DerivativeMaterialInterface<AuxKernel>(parameters),
    _strain_type(getParam<MooseEnum>("strain_type"))
{
  switch (_strain_type)
  {
    case 1:
      _strain_name = "strain_increment";
      break;
    case 2:
      _strain_name = "elastic_strain_increment";
      break;
    case 3:
      _strain_name = "plastic_strain_increment";
      break;
    default:
      mooseError("LMStrainAuxBase: unknown strain type!");
  }
  _strain_incr = &getDefaultMaterialProperty<RankTwoTensor>(_strain_name);
}

MooseEnum
LMStrainAuxBase::strainType()
{
  return MooseEnum("total=1 elastic=2 plastic=3");
}
