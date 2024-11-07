%% gör en implementation av den implicita trapetsmetoden (theta = 0.5)



% a) använd värdena från U4d, kör koden med tidssteget alpha*deltatmax med
% deltatmax som beräknats för euler framåt, för alpha = 1, 10, 100. vad
% händer? stämmer det med teorin?
z1 = 0; z2 = 0; zprick1 = 0; zprick2 = 0;
m1 = 465; m2 = 55; k1 = 5350; k2 = 13610000; c1 = 310; c2 = 1250; v = 63/3.6; H = 0.27; L = 1.1;
v_vec0 = [z1; z2; zprick1; zprick2];
T = 0.5;

h_max = 1.1458e-4;%0.0395; %deltatmax beräknat i U4
tspan = [0, T];

nZ=[];
nZT=[];
for n=0:2
    h = 10^n * h_max;

    t = 0;
    %Inv trapets
    [Z,ZT] = qc_inv_trap(tspan, h, v_vec0, k1, k2, c1, c2, m1, m2, H, L, v);

    switch n
        case 0
            n0Z=Z;
            n0ZT=ZT;
        case 1
            n1Z=Z;
            n1ZT=ZT;
        otherwise
            n2Z=Z;
            n2ZT=ZT;
    end
end

%Plot
hold on
a=2;
plot(n0ZT,n0Z(a,:))
plot(n1ZT,n1Z(a,:))
plot(n2ZT,n2Z(a,:))
options = odeset('RelTol',1e-9,'Refine',1);
[Rt,Rv] = ode45(@(t,y) quartercar(t, y, k1, k2, c1, c2, m1, m2, H, L, v),tspan,v_vec0,options);
plot(Rt, Rv(:,a))
ylim([0,0.3])
legend({"1*h_{max}", "10*h_{max}", "100*h_{max}", "ode45"})


%% b) konvergensstudie fram till t = 0.05 s, k2 = 100*k2,ref
% normen av felet = max absolutbelopp av felet i z2 i intervallet

%Beräknar exakt lösning
tspan=[0,0.05];
dt0 = 0.001;

[Rt,Rv] = ode45(@(t,y) quartercar(t, y, k1, k2, c1, c2, m1, m2, H, L, v),tspan,v_vec0,options);

%Olika dt
[e0Z,e0ZT] = qc_inv_trap(tspan, dt0, v_vec0, k1, k2, c1, c2, m1, m2, H, L, v);
[e1Z,e1ZT] = qc_inv_trap(tspan, dt0/2, v_vec0, k1, k2, c1, c2, m1, m2, H, L, v);
[e2Z,e2ZT] = qc_inv_trap(tspan, dt0/4, v_vec0, k1, k2, c1, c2, m1, m2, H, L, v);
[e3Z,e3ZT] = qc_inv_trap(tspan, dt0/8, v_vec0, k1, k2, c1, c2, m1, m2, H, L, v);

%Beräknar fel
e0 = max(abs(e0Z(2,end)-Rv(2,end)))
e1 = max(abs(e1Z(2,end)-Rv(2,end)))
e2 = max(abs(e2Z(2,end)-Rv(2,end)))
e3 = max(abs(e3Z(2,end)-Rv(2,end)))
%Konvergensordning
K=[];
E=[e3,e2,e1,e0];
for i=3:length(E)
    p = log(E(i-1)/E(i)) / log(E(i-2)/E(i-1));
    K=[K,p];
end

for i=1:length(K)
disp("Konvergens för " + i + ":a iterationen: " + K(i))
end

%Konvergens ~2, kvadratisk konvergens


% referenslösning med ode45, reltol abstol = 10^-9. ändra tspan för att få
% lösningen i sökta tidspunkter

% välj ett deltat0 för trapetsmetoden. kör koden med deltat = alpha*deltat0
%för alpha = 1, 0.5, 0.25, 0.125. beräkna felen och noggrannhetsordning
% inte rätt noggrannhetsordning? testa mindre deltat0
