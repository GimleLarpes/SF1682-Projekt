clear all
N = 100; % Spacial resolution
M = 1000; % Temporal resolution
T = 0.01; % End time
D = 1;
L=0;
R=0;

dx = 1/N;
dt = T/M;
x = linspace(0,1,N);

% Initial value
vec_twos = linspace(-2,-2,N-2)';
vec_ones = linspace(1,1,N-3)';
A = diag(vec_twos)+diag(vec_ones,1)+diag(vec_ones,-1);
g = zeros(N-2,1); g(1)=L; g(end)=R; % Randvärden
f = zeros(N-2,1); % Constant vector
w_vec = sin(5*pi * x)';

u_vec = w_vec(2:end-1);


% Step
w = [w_vec];
for t=18*dt:dt:T 

   u_vec = u_vec + dt * (D/(dx^2) * (A*u_vec + g) + f);
    
    w_vec = [L;u_vec;R];
    w = cat(2,w,w_vec);
end

%Visualize
surf(w)
%mesh(1:1001,x,w)

% DEL C
% Plot difference between analytical solution and numerical solution
w_numerical = w(:,end);
w_analytical = sin(5*pi * x') * exp(-25*pi^2*D * T);
error = abs(w_numerical - w_analytical);
max_error = max(error);
disp("Max error: "+max_error)
figure;
plot(x,error)




% DEL D
% Stabilitet för euler
max_dt = 0.5 * dx^2 / D;
disp("Max dt:"+max_dt)
% Stämmer överens, är stabil vid dt=1e-5 men inte stabil vid 1e-4.