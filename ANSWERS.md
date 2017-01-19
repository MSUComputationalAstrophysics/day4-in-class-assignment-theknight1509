Answers for in-class assignment january 19th.

The difference between relative change in history between Euler, EulerCromer and
the RungeKutta4 method can be viewed in 'relative_error.png'.

For these three methods to reach a relative error of 10^-4 or better is:
dt = 6*10^-2 for RK4 -> 2/3 *10^2 iterations
dt = 2*10^-2 for EC -> 2 *10^2 iterations
dt = 7*10^-4 for Euler -> 4/7 *10^4 iterations
dt = 10^-2 for midpoint -> 4 *10^2 iterations
dt = 10^-1 for PC -> 4*10 iterations

Flops for each method:
a is number of flops in acceleration function (1flop)
n is number of iterations requiired

RK4: (22+4a)n -> 1.8 *10^3 flops
EC: (4+a)n -> 10^3 flops
Euler: (4+a)n -> 20/7 *10^4 flops
midpoint: (10+2a)n -> 1.6*10^3 flops
PC: (10+2a)n -> 4.8 *10^2 flops

RK4 was the third best in total, which is a little disappointing.