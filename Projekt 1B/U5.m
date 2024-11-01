%% gör en implementation av den implicita trapetsmetoden (theta = 0.5)

% a) använd värdena från U4d, kör koden med tidssteget alpha*deltatmax med
% deltatmax som beräknats för euler framåt, för alpha = 1, 10, 100. vad
% händer? stämmer det med teorin?

%% b) konvergensstudie fram till t = 0.05 s, k2 = 100*k2,ref
% normen av felet = max absolutbelopp av felet i z2 i intervallet

% referenslösning med ode45, reltol abstol = 10^-9. ändra tspan för att få
% lösningen i sökta tidspunkter

% välj ett deltat0 för trapetsmetoden. kör koden med deltat = alpha*deltat0
%för alpha = 1, 0.5, 0.25, 0.125. beräkna felen och noggrannhetsordning
% inte rätt noggrannhetsordning? testa mindre deltat0
