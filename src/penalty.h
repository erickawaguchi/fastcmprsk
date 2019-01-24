#ifndef PENALTY_H
#define PENALTY_H

// This is the content of the .h file, which is where the declarations go
double getRidgeDelta(double grad, double hess, double a, double lam);
double getLassoDelta(double grad, double hess, double a, double lam);
double getScadDelta(double grad, double hess, double a, double lam, double gamma);
double getMcpDelta(double grad, double hess, double a, double lam, double gamma);
#endif
