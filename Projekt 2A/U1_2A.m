% b)
x = 0.75;
u = @(x) sin(exp(x));
u_prime_analytic = cos(exp(x)) * exp(x);

h_values = zeros(1,8);
u_prime_2 = zeros(1,8); u_prime_6 = zeros(1,8); u_prime_7 = zeros(1,8);
error_2 = zeros(1,8); error_6 = zeros(1,8); error_7 = zeros(1,8);

for k = 1:8
    h = 2^(-k);
    h_values(k) = h;

    u_prime_2(k) = (u(x+h)-u(x-h))/(2*h);
    u_prime_6(k) = (-u(x+2*h)+4*u(x+h)-3*u(x))/(2*h);
    u_prime_7(k) = (3*u(x)-4*u(x-h)+u(x-2*h))/(2*h);

    error_2(k) = abs(u_prime_2(k) - u_prime_analytic);
    error_6(k) = abs(u_prime_6(k) - u_prime_analytic);
    error_7(k) = abs(u_prime_7(k) - u_prime_analytic);
end

x_ref = logspace(0, 2.5, 100); %Skapar en linje som kommer ha lutning -2 i loglog-grafen
y_ref = 1.5*x_ref.^(-2);

loglog(1./h_values,error_2)
xlabel('1/h'); ylabel('e(1/h)'); title('felet e som funk. av 1/h'); hold on
loglog(1./h_values,error_6)
loglog(1./h_values,error_7)
loglog(x_ref, y_ref,'--')
legend({"formel 2", "formel 6","formel 7","referenslinje"})
