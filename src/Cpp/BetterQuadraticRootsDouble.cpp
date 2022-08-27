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

#include "BetterQuadraticRoots.h"
#include <cassert>
#include <cstdint>
#include <cmath>
#include <limits>

using namespace std;

namespace Lomont {
    namespace Numerical {
        namespace QuadraticEquation {

            namespace { // anonymous namespace

                //#region Per platform funcs

                double Sqrt(double v) { return sqrt(v); }
                int Sign(double v) { return (0 < v) - (v < 0); }
                double CopySign(double x, double y) { return copysign(x, y); }
                double Abs(double v) { return abs(v); }
                double Fma(double x, double y, double z) { return fma(x, y, z); }
                bool IsNaN(double v) { return isnan(v); }
                bool IsInfinity(double v) { return isinf(v); }
                bool IsFinite(double v) { return isfinite(v); }
                bool IsNormal(double v) { return isnormal(v); }

                const double NaN = std::numeric_limits<double>::quiet_NaN();
                const double MaxValue = std::numeric_limits<double>::max();

                /// <summary>
                /// Return x*2^exp
                /// </summary>
                /// <param name="x"></param>
                /// <param name="exp"></param>
                /// <returns></returns>
                double Scale2(double x, int exp) { return scalbn(x, exp); }

                // put in 64 bit type to make porting to other sizes easier
                uint64_t ToBits(double value)
                {
                    uint64_t i;
                    memcpy(&i, &value, sizeof(value));
                    return i;
                }

                //#endregion

                //#region per type values // double, double, half, others...

                const int PrecisionBits = 53; // # of bits in binary format + implicit 1 bit
                const int ExponentBits = 11;
                const int TotalBits = sizeof(double) * 8;
                const int ExpMask = (1 << ExponentBits) - 1;
                const int ExpBias = (1 << (ExponentBits - 1)) - 1;

                // #endregion

                    /// <summary>
                    /// Compute a*b-c*d
                    /// </summary>
                    /// <param name="a"></param>
                    /// <param name="b"></param>
                    /// <param name="c"></param>
                    /// <param name="d"></param>
                    /// <returns></returns>
                double Det2X2(double a, double b, double c, double d)
                {
                    auto v1 = c * d; // c * d loses precision
                    auto v2 = Fma(-c, d, v1); // c*d full precision - c*d low precision gives excess
                    auto v3 = Fma(a, b, -v1); // ab-cd with lost precision
                    return v3 + v2; // add back excess
                }


                /*
                move all values into form (sign,exp, frac) with frac in 1<= frac < 2
                input value = sign*2^exp * value
                input of 0 returns (1,0,0)
                 */
                 // return int sign, int exp, double frac
                tuple<int, int, double> Normalize(double value)
                {
                    assert(IsFinite(value)); // do not call on NaN, Inf
                    if (value == 0) return { 1, 0, 0 };

                    uint64_t i = ToBits(value);

                    int sign = (i & (1ULL << (TotalBits - 1))) == 0 ? 1 : -1;
                    int exp = (int)((i >> (PrecisionBits - 1)) & ExpMask) - ExpBias;
                    auto frac = sign * Scale2(value, -exp);

                    if (!IsNormal(value))
                    {
                        // above not enough - do more....
                        auto [s2, e2, f2] = Normalize(frac);
                        exp += e2;
                        frac = f2;
                        assert(s2 == 1);
                    }

                    assert(1 <= frac && frac < 2.0f);

                    assert(sign * Scale2(frac, exp) == value);


                    return { sign, exp, frac };
                }

                /// <summary>
                /// Compute sqrt(|x/y|), handling overflow and underflow if possible. 
                /// </summary>
                /// <param name="x"></param>
                /// <param name="y"></param>
                /// <returns></returns>
                double DivRoot(double x, double y)
                {
                    auto [xS, xE, xF] = Normalize(x);
                    auto [yS, yE, yF] = Normalize(y);
                    assert(xS * yS >= 0);

                    auto q = xF / yF;
                    auto e = xE - yE;
                    if (((xE + yE) & 1) == 1)
                    {
                        // exponent odd, scale so can easily update after root
                        q = Scale2(q, 1);
                        e--;
                    }

                    auto r = Sqrt(q);
                    return Scale2(r, e / 2);
                }



                /// <summary>
                /// check special cases 
                /// </summary>
                /// <param name="a"></param>
                /// <param name="b"></param>
                /// <param name="c"></param>
                /// <returns></returns>
                /// return bool isHandled, double r1, double r2, RootType type
                tuple<bool, double, double, RootType> HandleSpecialCasesFloat(double a, double b, double c)
                {
                    if (IsNaN(a) || IsNaN(b) || IsNaN(c))
                        return { true, NaN, NaN, RootType::InputHasNaN };
                    if (IsInfinity(a) || IsInfinity(b) || IsInfinity(c))
                        return { true, NaN, NaN, RootType::InputHasInfinity };

                    // cases:
                    if (a == 0)
                    { // want bx+c = 0 gives x = -c/b
                        if (b == 0 && c == 0)
                            return { true, NaN, NaN, RootType::AllRealNumbers };

                        auto r1 = -c / b;
                        return { true, r1, NaN, RootType::OneRealRoot };
                    }

                    if (b == 0)
                    { // a != 0, want ax^2+c = 0, so x = +/- sqrt(-c/a)

                        auto sgn = Sign(a) * Sign(c); // sign of quotient

                        if (sgn <= 0)
                        { // real answers
                            auto r1 = DivRoot(-c, a);
                            auto r2 = -r1;
                            return { true, r1, r2, RootType::SuccessReal };
                        }
                        else
                        { // complex answers, purely imaginary
                            auto r2 = DivRoot(c, a);
                            return { true, 0, r2, RootType::SuccessComplex }; // 0 +/- i*r2
                        }
                    }

                    if (c == 0)
                    { // a,b != 0, of form ax^2 + bx = 0, so roots are x=0 and x=-b/a
                        return { true, 0, -b / a, RootType::SuccessReal };
                    }

                    return { false, 0, 0, RootType::SuccessReal }; // not real, but will continue to work

                }



                /// <summary>
                /// Compute the discriminant D = b*b-4*a*c
                /// Return the (scaled) root r' = Sqrt(|D|), if d >= 0, and a scaling factor E
                /// such that the correct root is r = 2^E * r'
                /// </summary>
                /// <param name="a"></param>
                /// <param name="b"></param>
                /// <param name="c"></param>
                /// <returns></returns>
                /// return (double root, bool nonnegative, int scale)
                tuple<double, bool, int> DiscriminantInfo(double a, double b, double c)
                {
                    auto [aS, aE, aF] = Normalize(a);
                    auto [bS, bE, bF] = Normalize(b);
                    auto [cS, cE, cF] = Normalize(c);

                    double root = b;
                    int scale = 0;
                    bool nonnegative = true;

                    if (2 * bE > aE + cE + PrecisionBits + 5) // +5 works, is derived, seems to work( +4, +0, -2, ) , -10 fails (-10, -5, -4, -3)
                    {
                        root = bF;
                        scale = bE;
                        nonnegative = true;
                    }
                    else if (2 * bE < aE + cE - PrecisionBits - 1) // works: (-1,+2,+4), fails (+5, +7,+15,+40)
                    {
                        scale = aE + cE;
                        if ((scale & 1) != 0) // is odd
                        {
                            scale--;
                            aF = Scale2(aF, 1); // move factor back in
                        }
                        scale = scale / 2 + 1; //  +1 for the 4 in 4ac, then root, /2 is root

                        root = Sqrt(aF * cF);
                        nonnegative = aS * cS < 0;
                    }
                    else
                    {
                        // from above, we have:
                        assert(
                            -PrecisionBits - 1 <= 2 * bE - aE - cE &&
                            2 * bE - aE - cE <= PrecisionBits + 5
                        );

                        // now must align exponents for b*b and a*c so they can be subtracted... 
                        // idea, pull midpoint of (be + be) and (ae + ce) to zero, making values close to 1.0f, but still same scale

                        // scale a and c exponents to center of a and c, leaves a*c fixed in value, makes robust against underflow/overflow on following scaling
                        // in effect, we will scale (a,c) = (a*2^-dc, c*2^dc), brings them to near same in size
                        auto deltaE = (aE - cE) / 2;

                        auto mid = (bE + bE + aE + cE) / 4;
                        aF = Scale2(a, -mid + 2 - deltaE); // add 2 to handle 4 in 4ac
                        bF = Scale2(b, -mid);
                        cF = Scale2(c, -mid + deltaE);

                        assert(IsFinite(aF) && IsFinite(bF) && IsFinite(cF));

                        auto d = Det2X2(bF, bF, aF, cF);

                        root = Sqrt(Abs(d));
                        nonnegative = d >= 0;
                        scale = mid;
                    }
                    return { root, nonnegative, scale };
                }
            } // anonymous namespace


            std::tuple<double, double, RootType > DoubleRoots(double a, double b, double c)
            {
                auto [isHandled, r1, r2, type] = HandleSpecialCasesFloat(a, b, c);
                if (isHandled)
                    return { r1, r2, type };

                // so now can assume a,b,c nonzero
                auto [root, nonnegative, rootE] = DiscriminantInfo(a, b, c);

                // todo - can make this all slightly more accurate, see https://lomont.org/posts/2022/a-better-quadratic-formula-algorithm/
                root = Scale2(root, rootE);

                if (nonnegative)
                {

                    if (Abs(b) < MaxValue / 2)
                        r1 = (-b - CopySign(root, b)) / Scale2(a, 1);
                    else
                        r1 = -b / Scale2(a, 1) - CopySign(root, b) / Scale2(a, 1);
                    r2 = c / (r1 * a);
                    return { r1, r2, RootType::SuccessReal };
                }
                else
                {
                    r1 = -b / Scale2(a, 1);
                    r2 = root / Scale2(a, 1);
                    return { r1, r2, RootType::SuccessComplex };
                }
            }

        }
    }
}

