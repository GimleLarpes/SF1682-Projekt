% Parametrar
N_ref = 2048;                   
t = 0.1;               
x_ref = linspace(0, 1, N_ref) ;  
D = 1;

% initialvillkor
u0_ref = etafcn(x_ref,0.05);
u0_hat_ref = fftshift(fft(u0_ref)) / N_ref;

k = (-N_ref/2:(N_ref/2-1));
u_hat_t_ref = u0_hat_ref .* exp(-D * (2*pi*k).^2 * t);  % FFT av lösningen vid tid t

u_fft_ref = real(ifft(ifftshift(u_hat_t_ref)*N_ref));

% Lösningen för N = 2^m, m = 4, 5, 6, 7, 8, 9.
m = 4:9;
convergence = zeros(1,length(m));
error = zeros(1,length(m));

for i=1:length(m)
    N = 2^m(i);
    x = linspace(0, 1, N);  

    u0 = etafcn(x,0.05);
    u0_hat = fftshift(fft(u0)) / N;

    k = (-N/2:(N/2-1));
    u_hat_t = u0_hat .* exp(-D * (2*pi*k).^2 * t);  % FFT av lösningen vid tid t
    
    u_fft = real(ifft(ifftshift(u_hat_t)*N));

    u_fft_ref_interp = interp1(x_ref, u_fft_ref, x, 'linear'); % Estimerar värdet av u_fft_ref i
    % x-diskretiseringarna istället för x_ref. 
    error(i) = sqrt(sum(abs(u_fft_ref_interp-u_fft))/N);
    if i > 2
        convergence(i) = log(error(i-1)/error(i))/log(error(i-2)/error(i-1));
    end
end
max(error)
figure
plot(u_fft,'--'); hold on
plot(u_fft_ref_interp)
figure
loglog(2.^m,error)
title("Felet e som funktion av antalet diskretiseringar N")
ylabel("e(N)"); xlabel("N")
