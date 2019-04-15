// Marko Imbrisak
// 16/4/2019
// Code to generate midpoint, trapezoidal and Simpson integral estimates
// Usage:
// integrator <double>I(3,0.,1.,user_function);
// Inegrator takes the number of steps, integral's limits and the integrand
// Values are recovered by calling I.call() or by accessing
// I.midpoint_rule, I.trapezoidal_rule, or I.simpson_rule
// Their corresponding error estimates can be found in
//  I.midpoint_rule_error, I.trapezoidal_rule_error, I.simpson_rule_error

#include "integrators.cpp"

template <class R>
const R funk(const R& x){
    return 1./(1+x*x);
}
int main(){
    cout<< "Integral should equal pi/4 "<< M_PI_4<<endl;
    cout<< "Float precision"<< endl;
    integrator <float>J(100,0.,1.,funk);
    J.call(float(M_PI_4));
    cout<<endl;

    cout<< "double precision"<< endl;
    integrator <double>I(100,0.,1.,funk);
    I.call(M_PI_4);
    cout<<endl;


    cout<< "double precision, long"<< endl;
    integrator <double,long>L(10000,0.,1.,funk);
    L.call(M_PI_4);
    return 0.;
}
