clear all
N = 100; % Spacial resolution
t_vec = 0:0.05:10; % Points in time to evaluate
% u(x,t) = sin(exp(t) * x);

x = linspace(0,1,N);

% Calculate
u = [];
for t = t_vec
    u_vec = sin(exp(t) * x');

    u = cat(2,u,u_vec);
end

%Visualize
surf(u)


% Temperaturprofilen utvecklas mycket orealistiskt genom att med tiden,
% då den liknar mera en gas än en temperaturprofil. För att detta skulle
% uppnås i en riktig skav måste den värmas och kylas mycket snabbt.
% Man kan säga att istället för att utjämnas som allmänt förväntas gör den
% motsaten och blir mera kaotisk.
