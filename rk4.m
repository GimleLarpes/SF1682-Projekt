%% RK4
function [tv, yv]=rk4(f, tspan, y)
    tv = tspan(1);
    yv = y;

    for i=2:length(tspan)
        t = tspan(i);
        h = tspan(i) - tspan(i-1);

        s1 = f(t,y');
        s2 = f(t + 0.5*h, y' + 0.5*s1*h);
        s3 = f(t + 0.5*h, y' + 0.5*s2*h);
        s4 = f(t + h, y' + s3*h);
        
        y = y + h/6 * (s1 + 2*s2 + 2*s3 + s4)';

        tv = cat(1, tv, t);
        yv = cat(1, yv, y);
    end
end