# Integrator
Code to generate midpoint, trapezoidal and Simpson integral estimates of a function.
As an example, consider a function
```C++
template <class R>
const R funk(const R& x){
    return 1./(1+x*x);
}
```
Estimates of its integral can be computed as

```C++
integrator <double>I(100,0.,1.,funk);
```
## Accessing computed estimates
The computed estimates can be accessed indirectly with
 ```C++
I.call();
```
or directly by name
```C++
I.midpoint_rule;
I.trapezoidal_rule;
I.simpson_rule;
```

## Accessing computed error estimates
The computed estimates can be accessed indirectly with
 ```C++
I.call();
```
or directly by name
```C++
I.midpoint_rule_error;
I.trapezoidal_rule_error;
I.simpson_rule_error;
```
