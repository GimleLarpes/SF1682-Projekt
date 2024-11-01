% skriv en funktion quartercar.m med systemet från U1
% skriv en kod som löser z1,z2,zprick1,zprick2 numeriskt med:


% initialvillkor z1 = z2 = zprick1 = zprick2 = 0, k1, k2 deras tabellvärden

%% i) Euler framåt
% startvärden
z1=0; z2=0; zprick1=0; zprick2=0;
m1=465; m2=55; k1=5350; k2=136100; c1=310; c2=1250; v=63/3.6; H=0.27; L=1.1;
zEuler = [z1 z2 zprick1 zprick2];
T = 10;
n = 100;
h=T/n;

while t<T
    %Måste beräkna h,dh
    zvec = quartercar(v, h, dh, k1, k2, c1, c2, m1, m2);
    zEuler = zEuler + zvec*h;
    t = t + h;
end
zEuler


%% ii) ode45
z1=0; z2=0; zprick1=0; zprick2=0;
m1=465; m2=55; k1=5350; k2=136100; c1=310; c2=1250; v=63/3.6; H=0.27; L=1.1;
zode450 = [z1 z2 zprick1 zprick2];
tspan = [0 10];
options = odeset('RelTol',10^(-6), 'AbsTol',10^(-6),'Refine',1);

[t, zode45] = ode45(@(t,z) quartercar(m1, m2, k1, k2, c1, c2, v, H, L),tspan,zode450,options)
time_steps = length(t);



%% a) plotta lösningen för z1, z2 med ode45 med reltol = 10^-6 vad är z1max, z2max och vad innebär det för passagerarna?
hold on
plot(zode45(1))
plot(zode45(2))
% eller plotta i samma?

%% b) bifoga en plot som innehåller tidsstegen för ode45. hur stora är tidsstegen och hur förändras de med tiden?
% förklara förändringarna i tidsstegen över tid genom att studera plottarna
% för den numeriska lösningen och tidsstegen. (anpassa refine)

hold on
plot(zode45(1), zode45(2))
plot(time_steps)

%% c) jämför euler med ode45 för deltat = 5*10^-3 och deltat = 5*10^-4 
%(plotta lösningarna för z2 i samma figur) och kommentera.

z1=0; z2=0; zprick1=0; zprick2=0;
m1=465; m2=55; k1=5350; k2=136100; c1=310; c2=1250; v=63/3.6; H=0.27; L=1.1;
zEuler = [z1 z2 zprick1 zprick2];

h=[];

for h=[5*10^(-3) 5*10^(-4)]
while t<T
    zvec = quartercar(m1, m2, k1, k2, c1, c2, v, H, L);
    zEuler = zEuler + zvec*h;
    t = t + h;
end

zEuler

end
