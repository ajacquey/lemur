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

#include "LMStrainAuxBase.h"

InputParameters
LMStrainAuxBase::validParams()
{
  InputParameters params = AuxKernel::validParams();
  params.addClassDescription("Base class for outputting strain values.");
  params.addParam<MooseEnum>("strain_type",
                             LMStrainAuxBase::strainType() = "total",
                             "The type of the strain tensor to output.");
  return params;
}

LMStrainAuxBase::LMStrainAuxBase(const InputParameters & parameters)
  : AuxKernel(parameters), _strain_type(getParam<MooseEnum>("strain_type"))
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
      _strain_name = "viscous_strain_increment";
      break;
    case 4:
      _strain_name = "plastic_strain_increment";
      break;
    default:
      mooseError("LMStrainAuxBase: unknown strain type!");
  }
  _strain_incr = &getADMaterialProperty<RankTwoTensor>(_strain_name);
}

MooseEnum
LMStrainAuxBase::strainType()
{
  return MooseEnum("total=1 elastic=2 viscous=3 plastic=4");
}
