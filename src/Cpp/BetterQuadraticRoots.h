#pragma once
#ifndef BETTER_QUADRATIC_ROOTS
#define BETTER_QUADRATIC_ROOTS
//
//    MIT License
//    
//    Copyright (c) 2022 Chris Lomont
//    
//    Permission is hereby granted, free of charge, to any person obtaining a copy
//    of this software and associated documentation files (the "Software"), to deal
//    in the Software without restriction, including without limitation the rights
//    to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
//    copies of the Software, and to permit persons to whom the Software is
//    furnished to do so, subject to the following conditions:
//    
//    The above copyright notice and this permission notice shall be included in all
//    copies or substantial portions of the Software.
//    
//    THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
//    IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
//    FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
//    AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
//    LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
//    OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
//    SOFTWARE.
//

#include <tuple>

namespace Lomont{namespace Numerical{namespace QuadraticEquation{

// return type for roots
enum class RootType
{
    // most common cases, last two bits = 01
    SuccessReal      = 0xb'0'01, // got 2 real roots in r1 and r2
    SuccessComplex   = 0xb'1'01, // got 2 complex roots as r1 +/- i*r2

    // bad input, last 2 bits = 00
    InputHasNaN      = 0xb'0'00, // r1=r2=NaN
    InputHasInfinity = 0xb'1'00, // r1=r2=NaN

    // rare cases, last two bits 10
    OneRealRoot      = 0xb'0'0'10, // coeff 'a' was 0, one root r1 = -c/b, which may be infinite, r2 = NaN
    AllRealNumbers   = 0xb'0'1'10, // a=b=c=0, all real numbers are roots, r1=r2=NaN
};

// Compute roots using float32 (float) or float64 (double) for the quadratic equation ax^2 + bx + c = 0
// Returns (r1, r2, rootType) where root type is 
//   SuccessReal      : two real roots r1,r2
//   SuccessComplex   : two complex valued roots r1 +\- i*r2
//   InputHasNaN      : input has invalid values
//   InputHasInfinity : input has invalid values
//   OneRealRoot      : a was 0, so real root in r1, r2 = NaN
//   AllRealNumbers   : a=b=c=0, all numbers valid roots, r1=r2=NaN
//   
// Derivation of algorithms Chris Lomont, 2022, https://lomont.org/posts/2022/a-better-quadratic-formula-algorithm/


std::tuple<float , float , RootType > FloatRoots(float a, float b, float c);
//std::tuple<double, double, RootType > DoubleRoots(double a, double b, double c);

}}} // namespace

#endif // BETTER_QUADRATIC_ROOTS
