function ds=fvel(t, s, b, aL, aR, wL, wR)

x = s(1);
y = s(2);
theta = s(3);
vR = aR*t + wR;
vL = aL*t + wL;
v = (vR + vL)/2;

dx = v*cos(theta);
dy = v*sin(theta);

dtheta = (vR - vL)/b;

ds = [dx, dy, dtheta];

end
