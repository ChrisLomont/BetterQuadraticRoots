namespace Lomont.Numerical
{
    public static partial class QuadraticEquation
    {

        [Flags]
        public enum RootType
        {
            // most common cases, last two bits = 01
            SuccessReal       = 0xb_0_01, // got 2 real roots in r1 and r2
            SuccessComplex    = 0xb_1_01, // got 2 complex roots as r1 +/- i*r2

            // bad input, last 2 bits = 00
            InputHasNaN       = 0xb_0_00, // r1=r2=NaN
            InputHasInfinity  = 0xb_1_00, // r1=r2=NaN

            // rare cases, last two bits 10
            OneRealRoot       = 0xb_0_0_10, // coeff 'a' was 0, one root r1 = -c/b, which may be infinite, r2 = NaN
            AllRealNumbers    = 0xb_0_1_10, // a=b=c=0, all real numbers are roots, r1=r2=NaN
        }
    }
}
