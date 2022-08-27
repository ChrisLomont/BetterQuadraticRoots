# Better Quadratic Roots

Chris Lomont, lomont.org

v0.1 - Aug 22, 2022

Here is a much more robust algorithm for **computing roots of quadratic equations** than you will find almost anywhere. Here are accuracy results and fail rates for the float version over 1,000,000 randomly generated floats with exponents in various ranges:

| 1000000 passes each   | small range <br />(2^-32, 2^32) |          |        | large range<br />(2^-70,2^70) |           |        | huge range<br />(2^-126,2^126) |           |        |
| --------------------- | ------------------------------- | -------- | ------ | ----------------------------- | --------- | ------ | ------------------------------ | --------- | ------ |
| Function              | max ulp                         | avg ulp  | fail % | max ulp                       | avg ulp   | fail % | max ulp                        | avg ulp   | fail % |
| Naïve                 | 4.10E+12                        | 1.50E+07 | 0%     | 3.10E+23                      | -1.90E+17 | 4.10%  | 2.60E+37                       | -7.20E+30 | 34%    |
| Common                | 3.60E+04                        | 0.38     | 0%     | 1500                          | 0.34      | 3.90%  | 2.00E+22                       | 1.80E+17  | 25%    |
| Rare                  | 3.20                            | 0.36     | 0%     | 1500                          | 0.34      | 3.50%  | 2.00E+22                       | 1.80E+17  | 13%    |
| **Better (this one)** | **3.20**                        | **0.36** | **0%** | **3.00**                      | **0.33**  | **0%** | **3.2**                        | **0.31**  | **0%** |

Naïve is implementing a quadratic using the high school algorithm, which everyone knows is terrible. Common is using the trick to remove basic cancellation, which is the most used algorithm. Rare uses a special algorithm to compute the discriminant, and this version is really quite rare. But even it fails in many useful cases. **Better** is this one, which handles all sorts of numerical issues, including underflow, overflow, many types of cancellation, and all performs quite well (~4x more time than naive, ~70ns on my i7), and is small code (a few hundred lines). Fail is cases mischaracterized compared to truth data.

Goals were

1. **Accurate**. Keep overall error very low, achieved by keeping ulp error under 4 in teating, with an average of ~0.3 ulps over a wide range of floating point
2. **Robust**. Does not puke or produce bad results without warning. Also achieved
3. **Small code**, making it easy to jam into projects. Also achieved. It relies on no other software, making it easy to port. For example, I developed the float32 C# version first, then ported that to C# double, C# Half, C++ float32, and C++ double, all 4 ports totaling a few hours including testing.

Initial version are in the `src` folder, including a C# 7.0 version and C++17 versions for float32, float64, and float16 (for C# only).

I wrote how I developed it at https://lomont.org/posts/2022/a-better-quadratic-formula-algorithm/.

Have fun!


## TODO

1. post testing and development code
2. Demonstrate performance versus common other numerical libs