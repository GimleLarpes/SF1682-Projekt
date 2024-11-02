function dv=quartercar(t, v_vec, k1, k2, c1, c2, m1, m2, H, L, v)
    if t <= L/v
      h = H/2*(1-cos((2*pi*v*t)/L));
      dh = ((2*pi*v)/L)*sin((2*pi*v*t)/L);
      else
      h = 0;
      dh = 0;
    end
    
    H_matrix = [0, 0, 1, 0;
                0, 0, 0, 1;
               -k1 / m1, k1 / m1, -c1 / m1, c1 / m1;
                k1 / m1, -((k1 + k2) / m2), c1 / m2, -((c1 + c2) / m2)];
    g = [0; 0; 0; (c2*dh + k2*h) / m2];

    dv = H_matrix*v_vec + g;
end
