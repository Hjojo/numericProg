I choose to make the examination project number 1: "Yet another cubic sub-spline".

Since the cubic sub-spline have the same structure as the cubic spline from the interpolation exercise, the evaluation and free functions could be reused.
The spline alloc function did require some tweaks. The spline have the following form:
S_i(x) = y_i + b_i*(x-x_i) + c_i*(x-x_i)^2 + d_i*(x-x_i)^3 ,
where b_i is the derivative of the quadratic spline made around the points x_(i-1), x_i and x(i+1) at the point x_i. The quadratic spline was made with the made method from the interpolation exercise where the value b_2 correspond to the derivative at the point.
c_i and d_i could be found in the same ways as for the Akima sub-spline in the interpolation note.

Two figures are created, plotRND.svg and plotSin.svg, which compare the sub-spline, cubic spline and quadratic spline.
In the figures it can be seen that the sub-spline acts like a flatter version of the cubic spline.
However the sub-spline differ from the cubic spline by changing direction after the last point, whereas the cubic spline continue roughly in the same direction.
This is however no problem as long as it's only used for interpolation, as intended.
