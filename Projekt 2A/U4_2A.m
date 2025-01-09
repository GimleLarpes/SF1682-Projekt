%Definitionsmängder
xl = 0; xr = 1;
yb = 0; yt = 0.1;

N = (yt-yb) / 10^-3;     % antal steg i t
M = 100;                 % antal steg i x
h = (xr-xl) / M;         % steglängd i x
k = (yt-yb) / N;         % steglängd i t

%Dirichlet randvillkor
x0 = @(t) 0 * t;
x1 = @(t) sin(exp(t));

f = @(x, t) exp(t).*x.*cos(exp(t).*x) + exp(2*t).*sin(exp(t).*x);
D = 1;
sigma = D * k / (h^2);

%Diskretiseringar
x = xl + (0:M) * h;
t = yb + (0:N) * k;

A = diag(2+2*sigma*ones(M-1,1)) + diag(-sigma*ones(M-2,1),1) + diag(-sigma*ones(M-2,1),-1);
B = diag(2-2*sigma*ones(M-1,1)) + diag(sigma*ones(M-2,1),1) + diag(sigma*ones(M-2,1),-1);

lside = x0(yb + (0:N) * k);
rside = x1(yb + (0:N) * k);

w = zeros(M + 1, N + 1);      %Skapar lösningsmatris w med insatta R.V och B.V
w(:,1)= sin(x);
w(1, :) = lside;
w(end, :) = rside;

for j = 1:N
    sides = [lside(j) + lside(j+1); zeros(M-3, 1); rside(j) + rside(j+1)];
    f_eval = f(x(2:end-1),t(j))+f(x(2:end-1),t(1+j));

    w(2:end-1, j+1) = A \ (B * w(2:end-1, j) + sigma*sides + k*f_eval');
end

% Lösningar vid T=0.1
cn_solution = w(:, end);
analytic_solution = sin(exp(0.1) * x)';

error = abs(analytic_solution - cn_solution);
maxerror = max(error)

plot(x, error);
xlabel('x'); ylabel('e(x)');
title('Felet e i C-N:s metod vid tiden t=0,1');
grid on;