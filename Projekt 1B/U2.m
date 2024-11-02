% skriv en funktion quartercar.m med systemet från U1
% skriv en kod som löser z1,z2,zprick1,zprick2 numeriskt med:


% initialvillkor z1 = z2 = zprick1 = zprick2 = 0, k1, k2 deras tabellvärden

%% a) plotta lösningen för z1, z2 med ode45 med reltol = 10^-6 vad är z1max, z2max och vad innebär det för passagerarna?
z1=0; z2=0; zprick1=0; zprick2=0;
m1=465; m2=55; k1=5350; k2=136100; c1=310; c2=1250; v=63/3.6; H=0.27; L=1.1;
v_vec = [z1; z2; zprick1; zprick2];
tspan = [0 10];

options = odeset('RelTol',10^(-6),'Refine',1);

[t, zode45] = ode45(@(t, z) quartercar(t, z, k1, k2, c1, c2, m1, m2, H, L, v), tspan, v_vec, options);
time_steps = length(t)
zode45
%plot(zode45(:,1))


%% b) bifoga en plot som innehåller tidsstegen för ode45. hur stora är tidsstegen och hur förändras de med tiden?
% förklara förändringarna i tidsstegen över tid genom att studera plottarna
% för den numeriska lösningen och tidsstegen. (anpassa refine)

hold on
plot(zode45(1), zode45(2))
plot(time_steps)

%% c) jämför euler med ode45 för deltat = 5*10^-3 och deltat = 5*10^-4 
%(plotta lösningarna för z2 i samma figur) och kommentera.

%i) Euler framåt
z1=0; z2=0; zprick1=0; zprick2=0;
m1=465; m2=55; k1=5350; k2=136100; c1=310; c2=1250; v=63/3.6; H=0.27; L=1.1;
v_vec = [z1; z2; zprick1; zprick2];
t = 0;
T = 10;
delta_t = 5*10^-3;
i = 0;
z1_vec=[];
z2_vec=[];

while t<T
    i = i+1;

    dv = quartercar(t, v_vec, k1, k2, c1, c2, m1, m2, H, L, v);
    v_vec = v_vec + dv*delta_t;
    t = t + delta_t;

    z1_vec(i) = v_vec(1);
    z2_vec(i) = v_vec(2);
end

v_vec

plot(z1_vec)
hold on
plot(z2_vec)
legend({"z1", "z2"})
