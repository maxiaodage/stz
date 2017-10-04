#ifndef STZAPP_H
#define STZAPP_H

#include "MooseApp.h"

class StzApp;

template<>
InputParameters validParams<StzApp>();

class StzApp : public MooseApp
{
public:
  StzApp(InputParameters parameters);
  virtual ~StzApp();

  static void registerApps();
  static void registerObjects(Factory & factory);
  static void associateSyntax(Syntax & syntax, ActionFactory & action_factory);
};

#endif /* STZAPP_H */
