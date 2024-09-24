function [xbar, ybar, R, Tlap] = robotcircle(b, wL, wR, x0, y0, theta0)
    % Definierar B, D
    B = (wR + wL)/2;
    D = (wR - wL)/b;

    % Mittpunkt
    xbar = x0 - B/D * sin(theta0);
    ybar = y0 + B/D * cos(theta0);

    % Radie
    R = abs(B/D);

    % T_lap
    Tlap = 2*pi/D;
end