﻿//
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

using System.Diagnostics;

namespace Lomont.Numerical
{
    public static partial class QuadraticEquation
    {

        /// <summary>
        /// Compute roots using float32 for the quadratic equation ax^2 + bx + c = 0
        /// Returns (r1, r2, rootType) where root type is 
        ///   SuccessReal      : two real roots r1,r2
        ///   SuccessComplex   : two complex valued roots r1 +\- i*r2
        ///   InputHasNaN      : input has invalid values
        ///   InputHasInfinity : input has invalid values
        ///   OneRealRoot      : a was 0, so real root in r1, r2 = NaN
        ///   AllRealNumbers   : a=b=c=0, all numbers valid roots, r1=r2=NaN
        ///   
        /// Derivation of algorithms Chris Lomont, 2022, https://lomont.org/posts/2022/a-better-quadratic-formula-algorithm/
        /// </summary>
        /// <param name="a"></param>
        /// <param name="b"></param>
        /// <param name="c"></param>
        /// <returns></returns>
        public static (float r1, float r2, RootType type) FloatRoots(float a, float b, float c)
        {
            var (isHandled, r1, r2, type) = InternalF.HandleSpecialCasesFloat(a, b, c);
            if (isHandled)
                return (r1, r2, type);

            // so now can assume a,b,c nonzero
            var (root, nonnegative, rootE) = InternalF.DiscriminantInfo(a, b, c);

            // todo - can make this all slightly more accurate, see https://lomont.org/posts/2022/a-better-quadratic-formula-algorithm/
            root = InternalF.Scale2(root, rootE);

            if (nonnegative)
            {

                if (InternalF.Abs(b) < InternalF.MaxValue / 2)
                    r1 = (-b - InternalF.CopySign(root, b)) / InternalF.Scale2(a, 1);
                else
                    r1 = -b / InternalF.Scale2(a, 1) - InternalF.CopySign(root, b) / InternalF.Scale2(a, 1);
                r2 = c / (r1 * a);
                return (r1, r2, RootType.SuccessReal);
            }
            else
            {
                r1 = -b / InternalF.Scale2(a, 1);
                r2 = root / InternalF.Scale2(a, 1);
                return (r1, r2, RootType.SuccessComplex);
            }
        }

        static class InternalF
        {
            #region Per platform funcs

            internal static float Sqrt(float v) => float.Sqrt(v);
            internal static int Sign(float v) => float.Sign(v);
            internal static float CopySign(float x, float y) => float.CopySign(x, y);
            internal static float Abs(float v) => float.Abs(v);
            internal static float Fma(float x, float y, float z) => float.FusedMultiplyAdd(x, y, z);
            internal static bool IsNaN(float v) => float.IsNaN(v);
            internal static bool IsInfinity(float v) => float.IsInfinity(v);
            internal static bool IsFinite(float v) => float.IsFinite(v);
            internal static bool IsNormal(float v) => float.IsNormal(v);

            internal const float NaN = float.NaN;
            internal const float MaxValue = float.MaxValue;

            /// <summary>
            /// Return x*2^exp
            /// </summary>
            /// <param name="x"></param>
            /// <param name="exp"></param>
            /// <returns></returns>
            internal static float Scale2(float x, int exp) => float.ScaleB(x, exp);

            /// <summary>
            /// Compute a*b-c*d
            /// </summary>
            /// <param name="a"></param>
            /// <param name="b"></param>
            /// <param name="c"></param>
            /// <param name="d"></param>
            /// <returns></returns>
            internal static float Det2X2(float a, float b, float c, float d)
            {
                var v1 = c * d; // c * d loses precision
                var v2 = Fma(-c, d, v1); // c*d full precision - c*d low precision gives excess
                var v3 = Fma(a, b, -v1); // ab-cd with lost precision
                return v3 + v2; // add back excess
            }

            // put in 64 bit type to make porting to other sizes easier
            static ulong ToBits(float value) => BitConverter.SingleToUInt32Bits(value);

            #endregion

            #region per type values // double, float, half, others...

            const int PrecisionBits = 24; // # of bits in binary format + implicit 1 bit
            const int ExponentBits = 8;
            const int TotalBits = sizeof(float) * 8;
            const int ExpMask = (1 << ExponentBits) - 1;
            const int ExpBias = (1 << (ExponentBits - 1)) - 1;

            #endregion

            /*
            move all values into form (sign,exp, frac) with frac in 1<= frac < 2 
            input value = sign*2^exp * value
            input of 0 returns (1,0,0)
             */
            internal static (int sign, int exp, float frac) Normalize(float value)
            {
                Debug.Assert(IsFinite(value)); // do not call on NaN, Inf
                if (value == 0) return (1, 0, 0);

                ulong i = ToBits(value);

                int sign = (i & (1UL << (TotalBits-1))) == 0 ? 1 : -1;
                int exp = (int)((i >> (PrecisionBits-1)) & ExpMask) - ExpBias;
                var frac = sign * Scale2(value, -exp);

                if (!IsNormal(value))
                {
                    // above not enough - do more....
                    var (s2, e2, f2) = Normalize(frac);
                    exp += e2;
                    frac = f2;
                    Debug.Assert(s2 == 1);
                }

                Debug.Assert(1 <= frac && frac < 2.0f);

                Debug.Assert(sign * Scale2(frac, exp) == value);


                return (sign, exp, frac);
            }




            /// <summary>
            /// check special cases 
            /// </summary>
            /// <param name="a"></param>
            /// <param name="b"></param>
            /// <param name="c"></param>
            /// <returns></returns>
            internal static (bool isHandled, float r1, float r2, RootType type) HandleSpecialCasesFloat(float a, float b, float c)
            {
                if (IsNaN(a) || IsNaN(b) || IsNaN(c))
                    return (true, NaN, NaN, RootType.InputHasNaN);
                if (IsInfinity(a) || IsInfinity(b) || IsInfinity(c))
                    return (true, NaN, NaN, RootType.InputHasInfinity);

                // cases:
                if (a == 0)
                { // want bx+c = 0 gives x = -c/b
                    if (b == 0 && c == 0)
                        return (true, NaN, NaN, RootType.AllRealNumbers);

                    var r1 = -c / b;
                    return (true, r1, NaN, RootType.OneRealRoot);
                }

                if (b == 0)
                { // a != 0, want ax^2+c = 0, so x = +/- sqrt(-c/a)

                    var sgn = Sign(a) * Sign(c); // sign of quotient

                    if (sgn <= 0)
                    { // real answers
                        var r1 = DivRoot(-c, a);
                        var r2 = -r1;
                        return (true, r1, r2, RootType.SuccessReal);
                    }
                    else
                    { // complex answers, purely imaginary
                        var r2 = DivRoot(c, a);
                        return (true, 0, r2, RootType.SuccessComplex); // 0 +/- i*r2
                    }
                }

                if (c == 0)
                { // a,b != 0, of form ax^2 + bx = 0, so roots are x=0 and x=-b/a
                    return (true, 0, -b / a, RootType.SuccessReal);
                }

                return (false, 0, 0, RootType.SuccessReal); // not real, but will continue to work

            }


            /// <summary>
            /// Compute sqrt(|x/y|), handling overflow and underflow if possible. 
            /// </summary>
            /// <param name="x"></param>
            /// <param name="y"></param>
            /// <returns></returns>
            internal static float DivRoot(float x, float y)
            {
                var (xS, xE, xF) = Normalize(x);
                var (yS, yE, yF) = Normalize(y);
                Debug.Assert(xS*yS>=0);

                var q = xF / yF;
                var e = xE - yE;
                if (((xE + yE) & 1) == 1)
                {
                    // exponent odd, scale so can easily update after root
                    q = Scale2(q, 1);
                    e--;
                }

                var r = Sqrt(q);
                return Scale2(r, e / 2);
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
            internal static (float root, bool nonnegative, int scale) DiscriminantInfo(float a, float b, float c)
            {
                var (aS, aE, aF) = Normalize(a);
                var (bS, bE, bF) = Normalize(b);
                var (cS, cE, cF) = Normalize(c);

                var (root, scale, nonnegative) = (b, 0, true); // float,scaling, disc >= 0

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
                    Debug.Assert(
                        -PrecisionBits - 1 <= 2 * bE - aE - cE &&
                        2 * bE - aE - cE <= PrecisionBits + 5
                    );

                    // now must align exponents for b*b and a*c so they can be subtracted... 
                    // idea, pull midpoint of (be + be) and (ae + ce) to zero, making values close to 1.0f, but still same scale

                    // scale a and c exponents to center of a and c, leaves a*c fixed in value, makes robust against underflow/overflow on following scaling
                    // in effect, we will scale (a,c) = (a*2^-dc, c*2^dc), brings them to near same in size
                    var deltaE = (aE - cE) / 2;

                    var mid = (bE + bE + aE + cE) / 4;
                    aF = Scale2(a, -mid + 2 - deltaE); // add 2 to handle 4 in 4ac
                    bF = Scale2(b, -mid);
                    cF = Scale2(c, -mid + deltaE);

                    Debug.Assert(IsFinite(aF) && IsFinite(bF) && IsFinite(cF));

                    var d = Det2X2(bF, bF, aF, cF);

                    root = Sqrt(Abs(d));
                    nonnegative = d >= 0;
                    scale = mid;
                }
                return (root, nonnegative, scale);
            }
        }
    }
}

