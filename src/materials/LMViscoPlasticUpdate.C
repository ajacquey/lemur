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

#include "LMViscoPlasticUpdate.h"

InputParameters
LMViscoPlasticUpdate::validParams()
{
  InputParameters params = ADMaterial::validParams();
  params.addClassDescription("Base class for the viscoplastic correction.");
  params.set<bool>("compute") = false;
  params.suppressParameter<bool>("compute");
  params.addCoupledVar("fluid_pressure", 0, "The fluid pressure variable.");
  params.addRangeCheckedParam<Real>("abs_tolerance",
                                    1.0e-10,
                                    "abs_tolerance > 0.0",
                                    "The absolute tolerance for the iterative update.");
  params.addRangeCheckedParam<Real>("rel_tolerance",
                                    1.0e-10,
                                    "rel_tolerance > 0.0",
                                    "The relative tolerance for the iterative update.");
  params.addRangeCheckedParam<unsigned int>(
      "max_iterations",
      200,
      "max_iterations >= 1",
      "The maximum number of iterations for the iterative update");
  params.addRequiredRangeCheckedParam<Real>(
      "plastic_viscosity", "plastic_viscosity > 0.0", "The plastic viscosity.");
  params.addRangeCheckedParam<Real>(
      "exponent", 1.0, "exponent > 0.0", "The exponent for Perzyna-like flow rule.");
  params.addRangeCheckedParam<Real>("reference_fluid_pressure",
                                    1.0,
                                    "reference_fluid_pressure>0.0",
                                    "The reference fluid pressure.");
  params.addRangeCheckedParam<Real>(
      "Arrhenius_coefficient",
      0.0,
      "Arrhenius_coefficient>=0",
      "The Arrhenius coefficient for the fludi pressure activated viscosity.");
  return params;
}

LMViscoPlasticUpdate::LMViscoPlasticUpdate(const InputParameters & parameters)
  : ADMaterial(parameters),
    _pf(adCoupledValue("fluid_pressure")),
    _abs_tol(getParam<Real>("abs_tolerance")),
    _rel_tol(getParam<Real>("rel_tolerance")),
    _max_its(getParam<unsigned int>("max_iterations")),
    _eta_p(getParam<Real>("plastic_viscosity")),
    _n(getParam<Real>("exponent")),
    _pf0(getParam<Real>("reference_fluid_pressure")),
    _Ar(getParam<Real>("Arrhenius_coefficient")),
    _yield_function(declareADProperty<Real>("yield_function")),
    _plastic_strain_incr(declareADProperty<RankTwoTensor>("plastic_strain_increment"))
{
}

void
LMViscoPlasticUpdate::setQp(unsigned int qp)
{
  _qp = qp;
}