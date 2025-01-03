% Referenslösning
N_ref = 2048;                      % Antal punkter
t = 0.1;                           % Tidsvärde
x_ref = linspace(0, 1, N_ref);     % Diskretisering av x-intervallet [0, 1]

u0_ref = etafcn(x_ref,0.05);
u0_hat_ref = fft(u0_ref);

k = (0:N_ref-1);  % Vågrökningsnummer
k = k - floor(N_ref/2);  % Justera så att k går från -N/2 till N/2-1
u_hat_t_ref = u0_hat_ref .* exp(-4 * t * pi^2 * k.^2);  % FFT av lösningen vid tid t

u_fft_ref = ifft(u_hat_t_ref);

% Lösningen för N = 2^m, m = 4, 5, 6, 7, 8, 9.
m = 4:9;
error = zeros(1,length(m));

for i=1:length(m)
    N = 2^m(i);
    x = linspace(0, 1, N);

    u0 = etafcn(x,0.05);
    u0_hat = fft(u0);
    
    k = (0:N-1);
    k = k - floor(N/2);
    u_hat_t = u0_hat .* exp(-4 * t * pi^2 * k.^2);
    
    u_fft = ifft(u_hat_t);

    u_fft_ref_interp = interp1(x_ref, u_fft_ref, x, 'linear');
    error(i) = sqrt(sum(abs(u_fft_ref_interp-u_fft))/N);
end

loglog(2.^m,error)
title("Felet e som funktion av antalet diskretiseringar N")
ylabel("e(N)"); xlabel("N")