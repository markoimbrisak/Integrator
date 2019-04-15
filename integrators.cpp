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

#include <iostream>
#include <cmath>
using namespace std;

template <class R, class n=int>
class integrator{
    n N;
    R a,b;
    const R (*f)(const R&);
    public:
        R midpoint_rule,trapezoidal_rule,simpson_rule;
        R midpoint_rule_error, trapezoidal_rule_error, simpson_rule_error; 
        integrator(const n& number_of_steps, const R& lower_limit, const R& upper_limit, const R (*integrand)(const R&)){
            N = number_of_steps;
            a = lower_limit;
            b = upper_limit;
            f = integrand;
            calc();
        }

        void call(){
            cout << "midpoint integrator   ";
            repr(midpoint_rule,midpoint_rule_error);
            cout << "trapezoidal  integrator   ";
            repr(trapezoidal_rule,trapezoidal_rule_error);
            cout << "Simpson     integrator   ";
            repr(simpson_rule,simpson_rule_error);
        }

        void call(const R& expectation){
            cout << "midpoint integrator   ";
            repr(midpoint_rule,midpoint_rule_error,relative_error(midpoint_rule,expectation));
            cout << "trapezoidal  integrator   ";
            repr(trapezoidal_rule,trapezoidal_rule_error, relative_error(trapezoidal_rule,expectation));
            cout << "Simpson     integrator   ";
            repr(simpson_rule,simpson_rule_error,relative_error(simpson_rule,expectation));
        }

        void repr(const R&A, const R& B, const R&C){
            cout << A << " +/- " << B << " (" << C << "\%)" << endl;
        }

        void repr(const R&A, const R& B){
            cout << A << " +/- " << B << endl;
        }

    private:
        void calc(){
                compute_midpoint_method();
                compute_midpoint_rule_error();
                compute_simpson_method();
                compute_simpson_rule_error();
                compute_trapezoidal_rule_error();
                compute_trapezoidal_method();
            }
        
        R relative_error(const R& a, const R& b){
                return 100.*Rabs(a-b)/Rabs(b);
            }
            
        R Rabs(const R& x){
                if(x>0) return x;
                if(x<0) return -x;
            }
        
        R compute_midpoint_method(){
            R h = (b-a)/R(N);
            midpoint_rule = 0.;
            for(n i = 0; i<N; i++){
                midpoint_rule+= h*f(a+i*h);
            }
            return midpoint_rule;
        }

        R compute_trapezoidal_method(){
            R h = (b-a)/R(N);
            trapezoidal_rule = h*(f(a)+f(b))/2.;
            for(n i = 1; i<N; i++){
                trapezoidal_rule+= h*f(a+i*h);
            }
            return trapezoidal_rule;  
        }

        R compute_simpson_method(){
            R h = (b-a)/R(2*N);
            simpson_rule = h*(f(a)+f(b))/3.;
            for(n i = 1; i<N; i++){
                simpson_rule+= h*(4*f(a+(2*i-1)*h)+2*f(a+2*i*h))/3.;
            }
            simpson_rule+= 4*h*f(a+(2*N-1)*h)/3.;
            return simpson_rule;  
        }

        R compute_midpoint_rule_error(){
            R h = (b-a)/R(N);
            R max = 0;
            R d;
            for(n i = 0; i<N; i++){
                d = Rabs((f(a+(i+1)*h)+f(a+(i-1)*h)-2*f(a+i*h))/(h*h));
                if(d>max) max = d;
            }
            midpoint_rule_error = (b-a)*h*h*max/24.;  
            return midpoint_rule_error;  
        }
        
        R compute_trapezoidal_rule_error(){
            R h = (b-a)/R(N);
            R max = 0;
            R d;
            for(n i = 0; i<N; i++){
                d = Rabs((f(a+(i+1)*h)+f(a+(i-1)*h)-2*f(a+i*h))/(h*h));
                if(d>max) max = d;
            }
            trapezoidal_rule_error = (b-a)*h*h*max/12.;  
            return trapezoidal_rule_error;  
        }

        R compute_simpson_rule_error(){
            R h = (b-a)/R(2*N);
            R max = 0.;
            R d;
            for(n i = 0; i<N; i++){
                d = Rabs((f(a+(i+2)*h)+f(a+(i-2)*h)+6*f(a+i*h)-4*f(a+(i+1)*h)-4*f(a+(i-1)*h)))/(h*h*h*h);
                
                if(d>max) max = d;
            }
            simpson_rule_error = (b-a)*h*h*h*h*max/180.;  
            return simpson_rule_error;  
        }
    
};