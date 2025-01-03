alpha=400; 
xbar=pi/5;

m = 4:8;
err1 = zeros(1,length(m));
err2 = zeros(1,length(m));

for i=1:length(m)
    N=2^m(i);
    x=(0:N-1)/N;
    f=exp(-alpha*(x-xbar).^2); 
    % Analytiska lösningar
    fprim_a=-2*alpha*(x-xbar).*exp(-alpha*(x-xbar).^2);
    fbis_a=((-2*alpha*(x-xbar)).^2-2*alpha).*exp(-alpha*(x-xbar).^2);

    fk=fftshift(fft(f)); 
    kvec=-N/2:N/2-1; 
    % Derivatan av Fouriertransform
    fkprim = (1i * kvec*(2*pi)) .* fk;
    fkbis = ((1i * kvec*(2*pi)).^2) .* fk;
    
    fprim = real(ifft(ifftshift(fkprim)));
    fbis = real(ifft(ifftshift(fkbis)));
    % Diskreta 2-normen
    err1(i) = sqrt(sum(abs(fprim_a-fprim))/N);
    err2(i) = sqrt(sum(abs(fbis_a-fbis))/N);
end

loglog(2.^m,err1); hold on; loglog(2.^m,err2)
legend("f'","f''")
title("Felet e som funktion av antalet diskretiseringar N")
ylabel("e(N)"); xlabel("N")