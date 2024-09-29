%%RK4
function [tv, yv]=rk4(f, tspan, y, n)
    t = tspan(1);
    h = (tspan(2) - tspan(1)) / n;

    tv=[];
    yv=[];
    while t < tspan(2)
        s1 = f(t,y);
        s2 = f(t + 0.5*h, y + 0.5*s1*h);
        s3 = f(t + 0.5*h, y + 0.5*s2*h);
        s4 = f(t + h, y + s3*h);
        
        y = y + h/6 * (s1 + 2*s2 + 2*s3 + s4);
        t = t + h;

        tv = cat(1, tv, t);
        yv = cat(1, yv, y');
    end
end