function [vl,tl]=qc_inv_trap(tspan, h, v_vec, k1, k2, c1, c2, m1, m2, H, L, v)
    
    H_matrix = [0, 0, 1, 0;
                0, 0, 0, 1;
               -k1 / m1, k1 / m1, -c1 / m1, c1 / m1;
                k1 / m1, -((k1 + k2) / m2), c1 / m2, -((c1 + c2) / m2)];

    vl=[];
    tl=[];
    dt=h;
    t=tspan(1);
    while t<tspan(2)-h
        %Adjust step length to not overshoot
        if tspan(2) - t < dt
            dt = tspan(2) - t + eps;
        end
        
        %Calculate stuff
        g = calc_g(t, k2, c2, m2, H, L, v);

        dv = H_matrix*v_vec + g;

        %v_vec = v_vec + 0.5*dt * (H_matrix*v_vec + H_matrix*(v_vec + dv) + g + calc_g(t+dt, k2, c2, m2, H, L, v));
        theta=0.5;
        v_vec = (eye-dt*theta * H_matrix)\(v_vec + dt*((1-theta) * (H_matrix * v_vec+g)) + dt * theta * calc_g(t+dt, k2, c2, m2, H, L, v));

        vl = cat(2, vl, v_vec);
        tl = cat(2, tl, t);
        t=t+dt;
    end
end

function g=calc_g(t, k2, c2, m2, H, L, v)
    if t <= L/v
        h = H/2*(1-cos((2*pi*v*t)/L));
        dh = ((2*pi*v)/L)*sin((2*pi*v*t)/L);
    else
        h = 0;
        dh = 0;
    end
    g = [0; 0; 0; (c2*dh + k2*h) / m2];
end