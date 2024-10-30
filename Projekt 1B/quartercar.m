%% quartercar.m
function dv=quartercar(v, h, dh, k1, k2, c1, c2, m1, m2)
    H = [0, 0, 1, 0;
         0, 0, 0, 1;
         -k1/m1, k1/m1, -c1/m1, c1/m1;
         k1/m1, -(k1 + k2) / m2, c1/m2, -(c1 + c2) / m2];
    g = [0; (c2*dh + k2*h) / m2];

    dv = H*v + g;
end