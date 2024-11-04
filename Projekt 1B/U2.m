%% a) plotta lösningen för z1, z2 med ode45 med reltol = 10^-6 vad är z1max, z2max och vad innebär det för passagerarna?
z1=0; z2=0; zprick1=0; zprick2=0;
m1=465; m2=55; k1=5350; k2=136100; c1=310; c2=1250; v=63/3.6; H=0.27; L=1.1;
v_vec = [z1; z2; zprick1; zprick2];
tspan = [0 0.5];

options = odeset('RelTol',10^(-6),'Refine',1);

[t, zode45] = ode45(@(t, z) quartercar(t, z, k1, k2, c1, c2, m1, m2, H, L, v), tspan, v_vec, options);

figure(1);
plot(zode45(:,1))
hold on
plot(zode45(:,2))
legend({"z1", "z2"})
title('Förflyttningarna z1 och z2 som funktion av tiden t')
xlabel('t')
ylabel('z(t)')

% b) bifoga en plot som innehåller tidsstegen för ode45. hur stora är tidsstegen och hur förändras de med tiden?
% förklara förändringarna i tidsstegen över tid genom att studera plottarna
% för den numeriska lösningen och tidsstegen. (anpassa refine)
time_steps = length(t)
figure(2);
plot(t)
title('Steglängden h som funktion av tiden t')
xlabel('t')
ylabel('h(t)')

%% c) jämför euler med ode45 för deltat = 5*10^-3 och deltat = 5*10^-4 
%(plotta lösningarna för z2 i samma figur) och kommentera.

%i) Euler framåt
z1=0; z2=0; zprick1=0; zprick2=0;
m1=465; m2=55; k1=5350; k2=136100; c1=310; c2=1250; v=63/3.6; H=0.27; L=1.1;
T = 0.5;
z2_vec_1=[];
z2_vec_2=[];

% Euler med steglängd 5*10^-3
delta_t = 5*10^-3;
v_vec = [z1; z2; zprick1; zprick2];
i = 0;
t = 0;

while t<T
    i = i+1;

    dv = quartercar(t, v_vec, k1, k2, c1, c2, m1, m2, H, L, v);
    v_vec = v_vec + dv*delta_t;
    t = t + delta_t;

    z2_vec_1(i) = v_vec(1); 
end

plot(z2_vec_1)
hold on

% Euler med steglängd 5*10^-4
delta_t = 5*10^-4;
v_vec = [z1; z2; zprick1; zprick2];
i = 0;
t = 0;

while t<T
    i = i+1;

    dv = quartercar(t, v_vec, k1, k2, c1, c2, m1, m2, H, L, v);
    v_vec = v_vec + dv*delta_t;
    t = t + delta_t;

    z2_vec_2(i) = v_vec(1); 
end

plot(z2_vec_2)

% ode45
v_vec = [z1; z2; zprick1; zprick2];
tspan = [0 0.5];

options = odeset('RelTol',10^(-6),'Refine',1);

[t, zode45] = ode45(@(t, z) quartercar(t, z, k1, k2, c1, c2, m1, m2, H, L, v), tspan, v_vec, options);
time_steps = length(t)

plot(zode45(:,2))
legend({"Euler 5*10^-3", "Euler 5*10^-4", "ode45"})
title('Förflyttningen z2 som funktion av tiden t')
xlabel('t')
ylabel('z2(t)')
