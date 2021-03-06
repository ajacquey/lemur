#include "LemurApp.h"
#include "Moose.h"
#include "AppFactory.h"
#include "ModulesApp.h"
#include "MooseSyntax.h"

InputParameters
LemurApp::validParams()
{
  InputParameters params = MooseApp::validParams();

  // Do not use legacy DirichletBC, that is, set DirichletBC default for preset = true
  params.set<bool>("use_legacy_dirichlet_bc") = false;
  
  return params;
}

LemurApp::LemurApp(InputParameters parameters) : MooseApp(parameters)
{
  LemurApp::registerAll(_factory, _action_factory, _syntax);
}

LemurApp::~LemurApp() {}

static void
associateSyntaxInner(Syntax & syntax, ActionFactory & /*action_factory*/)
{
  registerSyntax("EmptyAction", "BCs/LMPressure");
  registerSyntax("LMPressureAction", "BCs/LMPressure/*");
}

void
LemurApp::registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  ModulesApp::registerAll(f, af, s);
  Registry::registerObjectsTo(f, {"LemurApp"});
  Registry::registerActionsTo(af, {"LemurApp"});

  /* register custom execute flags, action syntax, etc. here */
  associateSyntaxInner(s, af);
}

void
LemurApp::registerApps()
{
  registerApp(LemurApp);
}

/***************************************************************************************************
 *********************** Dynamic Library Entry Points - DO NOT MODIFY ******************************
 **************************************************************************************************/
extern "C" void
LemurApp__registerAll(Factory & f, ActionFactory & af, Syntax & s)
{
  LemurApp::registerAll(f, af, s);
}
extern "C" void
LemurApp__registerApps()
{
  LemurApp::registerApps();
}
