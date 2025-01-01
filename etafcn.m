function etav = etafcn(yv,delta,N)
a = 1-2*delta;
hfunc = @(s) exp(2*exp(-1/s)/(s-1));

etav = zeros(size(yv));
for i = 1:length(yv)
    y = yv(i);
    if 0 <= y && round(y*N) <= delta*N
        etav(i) = 1;
    elseif delta*N < round(y*N) && round(y*N) < (1-delta)*N
        etav(i) = hfunc((y-delta)/a)/(hfunc((y-delta)/a)+hfunc((1-(y-delta)/a)));
    elseif (1-delta)*N <= round(y*N) && round(y*N) <= 1*N
        etav(i) = 0;
    end
end

return