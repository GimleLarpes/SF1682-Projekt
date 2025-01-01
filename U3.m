clear all

% U3
% DEL A
D=1;
% Time
T=0.1; dt=1e-3;

% Space (-L/2,L/2)
L=1; N=16;
x=L*[0:N-1]/N;
dx=L/N;
% Frequency
xf = (2*pi/L) * [-N/2:N/2-1]; xf = fftshift(xf);


% Initial condition: u_0(x) = f(x)
beta = 6;
u_0 = cos(beta*pi * x);
uf_0 = fft(u_0);


% Solve in fourier space
t=0:dt:T;
[t,uf] = rk4(@(t, uf) HeatEQf(t, uf, xf, D) + gf(t,x,beta), t, uf_0);


% Inverse fourier
for k=1:length(t)
    u(k,:) = real(ifft(uf(k,:)));
end

% Plot solution
figure;
p=surf(u); xlabel('x'); ylabel('t'); zlabel('u')
title("Del A");


% ERROR
u_analytical = cos(beta * pi*x) * exp(-T);
error_vec = abs(u(end,:) - u_analytical);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% MedelFel då dt=1e-3:   5.68e-4, Max: 9.04e-4
%             dt=0.5e-3: 2.84e-4, Max: 4.52e-4
%             Felet avtar alltså som förväntat




% DEL B - Stabilitetsgräns
% Stabilitetsområdet för RK4 längs reel linje: [−2.8, 0]
dt_max = 0.28/N^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Max tidssteg för N=16:
%             dt_max = 0.28/N^2
%             Detta stämmer med testade värden på dt




% DEL C
VISRES=5;
% Time
T=0.1; dt=0.5*dt_max;

% Initial condition: u_0(x) = f(x)
delta = 0.05;
u_0 = x*0;
uf_0 = fft(u_0);


% Solve in fourier space
t=0:dt:T;
[t,ucf] = rk4(@(t, uf) HeatEQf(t, uf, xf, D) + gcf(t,x,delta), t, uf_0);


% Inverse fourier
for k=1:length(t)
    uc(k,:) = real(ifft(ucf(k,:)));
end

% Plot, minst 25 punkter
figure;
p=waterfall(x,t(1:VISRES:end),(uc(1:VISRES:end,:))); set(p,'LineWidth',5,'FaceAlpha',0.5);
xlabel('x'); ylabel('t'); zlabel('u')
title("Del C");




% FUNCTIONS
% Fourier'd heat equation
function dudt = HeatEQf(t, u_h, xf, D)
  dudt = -D^2 * (xf.^2)'.*u_h;
end

% Fourier'd g
function gf_vec = gf(t,x,beta)
    g = @(t,x) (beta^2*pi^2 - 1) * cos(beta*pi * x') .* exp(-t); 
    gf_vec = fft(g(t,x));
end

% Fourier'd g in C
function gf_vec = gcf(t,x,delta)
    N = 400; % This is different from N in the main function
    g = @(t,x) 100 * etafcn(2*abs(x' - 1/2), delta, N) .* exp(-100*t); 
    
    gf_vec = fft(g(t,x));
end