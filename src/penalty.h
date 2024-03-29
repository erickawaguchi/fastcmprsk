#ifndef PENALTY_H
#define PENALTY_H

// This is the content of the .h file, which is where the declarations go
double getRidge(double grad, double hess, double a, double lam);
double getLasso(double grad, double hess, double a, double lam);
double getScad(double grad, double hess, double a, double lam, double gamma);
double getMcp(double grad, double hess, double a, double lam, double gamma);
double getElasticNet(double grad, double hess, double a, double lam, double alpha);
double getGroupLasso(double z, double hess, double lam);
double getGroupScad(double z, double hess, double lam, double gamma);
double getGroupMcp(double z, double hess, double lam, double gamma);
#endif
