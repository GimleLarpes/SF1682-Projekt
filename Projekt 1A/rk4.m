%%RK4
function [tv, yv]=rk4(f, tspan, y, n)
    t = tspan(1);
    h = (tspan(2) - tspan(1)) / n;

    tv=[];
    yv=[];
    while t < tspan(2)
        s1 = f(t,y);
        s2 = f(t + 0.25*h, y + 0.25*s1*h);
        s3 = f(t + 0.375*h, y + (0.09375*s1 + 0.28125*s2)*h);
        s4 = f(t + 0.9231*h, y + (0.8794*s1 - 3.2772*s2 + 3.321*s3)*h);
        s5 = f(t + h, y + (2.0324*s1 -8*s2 + 7.1735*s3 - 0.2059*s4)*h);
        
        y = y + h*(0.1157*s1 + 0.5489*s3 + 0.5353*s4 - 0.2*s5);
        
        tv = cat(1, tv, t);
        yv = cat(1, yv, y');
        t = t + h;
    end
end