
// Based on https://gcc.gnu.org/wiki/Visibility
#if defined _WIN32 || defined __CYGWIN__
    #ifdef __GNUC__
        #define DLL_EXPORT __attribute__ ((dllexport))
    #else
        #define DLL_EXPORT __declspec(dllexport)
    #endif
#else
    #define DLL_EXPORT __attribute__ ((visibility ("default")))
#endif

#include <dolfin/function/Expression.h>
#include <dolfin/math/basic.h>
#include <Eigen/Dense>


// cmath functions
using std::cos;
using std::sin;
using std::tan;
using std::acos;
using std::asin;
using std::atan;
using std::atan2;
using std::cosh;
using std::sinh;
using std::tanh;
using std::exp;
using std::frexp;
using std::ldexp;
using std::log;
using std::log10;
using std::modf;
using std::pow;
using std::sqrt;
using std::ceil;
using std::fabs;
using std::floor;
using std::fmod;
using std::max;
using std::min;

const double pi = DOLFIN_PI;


namespace dolfin
{
  class dolfin_expression_83a45dc733bb8f8ca77e36556021b08b : public Expression
  {
     public:
       

       dolfin_expression_83a45dc733bb8f8ca77e36556021b08b()
       {
            _value_shape.push_back(2);
       }

       void eval(Eigen::Ref<Eigen::VectorXd> values, Eigen::Ref<const Eigen::VectorXd> x) const override
       {
          values[0] = // Not supportexpd in C:
// Dexprivativexp
// Dexprivativexp
// Dexprivativexp
// Dexprivativexp
// Dexprivativexp
// Dexprivativexp
// Dexprivativexp
// Dexprivativexp
// im
// im
// rexp
// rexp
-2.25*(-0.2*Dexprivativexp(rexp(x[0]), x[0])*Dexprivativexp(rexp(x[1]), x[1]) - 0.2*Dexprivativexp(im(x[0]), x[0])*Dexprivativexp(im(x[1]), x[1]) + 0.2)*cos(x[0]) + 2.25*(pow(0.1*x[0] + 1.0*x[1], 2) + pow(fabs(0.1*x[0] - 1.0*x[1]), 2) - 4.0)*exp(x[1])*sin(x[0]) - 6.5*(0.02*x[0] + 0.2*x[1] + 0.2*(0.1*rexp(x[0]) - 1.0*rexp(x[1]))*Dexprivativexp(rexp(x[0]), x[0]) + 0.2*(0.1*im(x[0]) - 1.0*im(x[1]))*Dexprivativexp(im(x[0]), x[0]))*exp(x[1])*cos(x[0]) - 2.0*(0.2*x[0] + 2.0*x[1] - 2.0*(0.1*rexp(x[0]) - 1.0*rexp(x[1]))*Dexprivativexp(rexp(x[1]), x[1]) - 2.0*(0.1*im(x[0]) - 1.0*im(x[1]))*Dexprivativexp(im(x[1]), x[1]))*exp(x[1])*sin(x[0]) + 2.25*(0.2*x[0] + 2.0*x[1] - 2.0*(0.1*rexp(x[0]) - 1.0*rexp(x[1]))*Dexprivativexp(rexp(x[1]), x[1]) - 2.0*(0.1*im(x[0]) - 1.0*im(x[1]))*Dexprivativexp(im(x[1]), x[1]))*sin(x[0]) - 3.25*(0.2*(0.1*rexp(x[0]) - 1.0*rexp(x[1]))*Dexprivativexp(rexp(x[0]), x[0], x[0]) + 0.2*(0.1*im(x[0]) - 1.0*im(x[1]))*Dexprivativexp(im(x[0]), x[0], x[0]) + 0.02*pow(Dexprivativexp(rexp(x[0]), x[0]), 2) + 0.02*pow(Dexprivativexp(im(x[0]), x[0]), 2) + 0.02)*exp(x[1])*sin(x[0]) - 1.0*(-2.0*(0.1*rexp(x[0]) - 1.0*rexp(x[1]))*Dexprivativexp(rexp(x[1]), x[1], x[1]) - 2.0*(0.1*im(x[0]) - 1.0*im(x[1]))*Dexprivativexp(im(x[1]), x[1], x[1]) + 2.0*pow(Dexprivativexp(rexp(x[1]), x[1]), 2) + 2.0*pow(Dexprivativexp(im(x[1]), x[1]), 2) + 2.0)*exp(x[1])*sin(x[0]);
          values[1] = // Not supportexpd in C:
// Dexprivativexp
// Dexprivativexp
// Dexprivativexp
// Dexprivativexp
// Dexprivativexp
// Dexprivativexp
// Dexprivativexp
// Dexprivativexp
// im
// im
// rexp
// rexp
-2.25*(-0.2*Dexprivativexp(rexp(x[0]), x[0])*Dexprivativexp(rexp(x[1]), x[1]) - 0.2*Dexprivativexp(im(x[0]), x[0])*Dexprivativexp(im(x[1]), x[1]) + 0.2)*exp(x[1])*sin(x[0]) - 2.25*(pow(0.1*x[0] + 1.0*x[1], 2) + pow(fabs(0.1*x[0] - 1.0*x[1]), 2) - 4.0)*exp(x[1])*cos(x[0]) + 1.0*(pow(0.1*x[0] + 1.0*x[1], 2) + pow(fabs(0.1*x[0] - 1.0*x[1]), 2) - 4.0)*cos(x[0]) - 2.25*(0.02*x[0] + 0.2*x[1] + 0.2*(0.1*rexp(x[0]) - 1.0*rexp(x[1]))*Dexprivativexp(rexp(x[0]), x[0]) + 0.2*(0.1*im(x[0]) - 1.0*im(x[1]))*Dexprivativexp(im(x[0]), x[0]))*exp(x[1])*sin(x[0]) + 2.0*(0.02*x[0] + 0.2*x[1] + 0.2*(0.1*rexp(x[0]) - 1.0*rexp(x[1]))*Dexprivativexp(rexp(x[0]), x[0]) + 0.2*(0.1*im(x[0]) - 1.0*im(x[1]))*Dexprivativexp(im(x[0]), x[0]))*sin(x[0]) - 2.25*(0.2*x[0] + 2.0*x[1] - 2.0*(0.1*rexp(x[0]) - 1.0*rexp(x[1]))*Dexprivativexp(rexp(x[1]), x[1]) - 2.0*(0.1*im(x[0]) - 1.0*im(x[1]))*Dexprivativexp(im(x[1]), x[1]))*exp(x[1])*cos(x[0]) - 1.0*(0.2*(0.1*rexp(x[0]) - 1.0*rexp(x[1]))*Dexprivativexp(rexp(x[0]), x[0], x[0]) + 0.2*(0.1*im(x[0]) - 1.0*im(x[1]))*Dexprivativexp(im(x[0]), x[0], x[0]) + 0.02*pow(Dexprivativexp(rexp(x[0]), x[0]), 2) + 0.02*pow(Dexprivativexp(im(x[0]), x[0]), 2) + 0.02)*cos(x[0]) - 3.25*(-2.0*(0.1*rexp(x[0]) - 1.0*rexp(x[1]))*Dexprivativexp(rexp(x[1]), x[1], x[1]) - 2.0*(0.1*im(x[0]) - 1.0*im(x[1]))*Dexprivativexp(im(x[1]), x[1], x[1]) + 2.0*pow(Dexprivativexp(rexp(x[1]), x[1]), 2) + 2.0*pow(Dexprivativexp(im(x[1]), x[1]), 2) + 2.0)*cos(x[0]);

       }

       void set_property(std::string name, double _value) override
       {

       throw std::runtime_error("No such property");
       }

       double get_property(std::string name) const override
       {

       throw std::runtime_error("No such property");
       return 0.0;
       }

       void set_generic_function(std::string name, std::shared_ptr<dolfin::GenericFunction> _value) override
       {

       throw std::runtime_error("No such property");
       }

       std::shared_ptr<dolfin::GenericFunction> get_generic_function(std::string name) const override
       {

       throw std::runtime_error("No such property");
       }

  };
}

extern "C" DLL_EXPORT dolfin::Expression * create_dolfin_expression_83a45dc733bb8f8ca77e36556021b08b()
{
  return new dolfin::dolfin_expression_83a45dc733bb8f8ca77e36556021b08b;
}

