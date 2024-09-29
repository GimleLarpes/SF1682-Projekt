function [x,y] = analytic_solutions(t,B,D,theta0,x0,y0)
    theta = D*t + theta0;
    x = x0 + (B/D)*(sin(theta)-sin(theta0));
    y = y0 - (B/D)*(cos(theta)-cos(theta0)); 
end