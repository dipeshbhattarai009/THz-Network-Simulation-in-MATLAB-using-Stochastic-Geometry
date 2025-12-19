% Plotting and comparing CCDF of SNR (Coverage Probability) using simulation and analytical
% derivation

clear;
lambda=0.001; %node density
alpha=2; %path loss exponent
xmin=-100;
xmax=100;
ymin=-100;
ymax=100;

%modelcular absorption coeffiecient (cm^-1)
kf1=0.1;
kf2=0.15;
kf3=0.2;



No=0.0001; %arbitrary value for Noise power
iteration=5000;
n=linspace(0.0001,10,50); %SNR threshold values for comparision

%calculation of coverage probability distribution using simulation
for i=1:iteration
    %number of points
    N=poissrnd(lambda.*(xmax-xmin).*(ymax-ymin));

    %locating transmitting nodes
    X1=unifrnd(xmin,xmax,1,N);
    Y1=unifrnd(xmin,xmax,1,N);

    %calculating distance of each point from reciever at origin and sorting
    %them in ascending order
    distances=sqrt(X1.^2+Y1.^2);
    distances=sort(distances);

    %Rayleigh fading coefficient
    u=exprnd(1);

    %Formula for calculating CCDF of SNR for RF
    F1(i,:)=n<((u*distances(1)^(-alpha))/No);

    %Formula for calculating CCDF of SNR for THz with different absorption
    %coefficient
    F2(i,:)=n<((distances(1)^(-alpha)*exp(-kf1*distances(1)))/No);
    F3(i,:)=n<((distances(1)^(-alpha)*exp(-kf2*distances(1)))/No);
    F4(i,:)=n<((distances(1)^(-alpha)*exp(-kf3*distances(1)))/No);
end

simulated_CCDF_SNR1=mean(F1);
simulated_CCDF_SNR2=mean(F2);
simulated_CCDF_SNR3=mean(F3);
simulated_CCDF_SNR4=mean(F4);

%computation of coverage probability distribution using mathematical
%equations
for m=1:length(n)
    f1 = @(r) exp(-n(m)*No*r.^(alpha)).*(2*pi*lambda*r.*exp(-lambda*pi*r.^2));
    computed_CCDF_SNR1(m) = integral(f1,0,inf);

    computed_CCDF_SNR2(m)=1-exp(-lambda.*pi.*(alpha/kf1.*lambertw(kf1/alpha.*((n(m).*No).^(-1/alpha)))).^2);
    computed_CCDF_SNR3(m)=1-exp(-lambda.*pi.*(alpha/kf2.*lambertw(kf2/alpha.*((n(m).*No).^(-1/alpha)))).^2);
    computed_CCDF_SNR4(m)=1-exp(-lambda.*pi.*(alpha/kf3.*lambertw(kf3/alpha.*((n(m).*No).^(-1/alpha)))).^2);
    
end


%plotting the simulated and computed coverage probability distribution for 
% RF and Thz (i.e. with and without considering molecular absorption 
% affect and Rayleigh fading) for different values of absorption
% coefficient

figure;
hold on;
plot(n,simulated_CCDF_SNR1,'ro','DisplayName','simulated without absorption (RF)');
plot(n,computed_CCDF_SNR1,'r-','DisplayName', "computed without absorption (RF)");
plot(n,simulated_CCDF_SNR2,'bo','DisplayName',"simulated (kf=0.1)");
plot(n,computed_CCDF_SNR2,'b-','DisplayName',"computed (kf=0.1");
plot(n,simulated_CCDF_SNR3,'go','DisplayName',"simulated (kf=0.15");
plot(n,computed_CCDF_SNR3,'g-','DisplayName',"computed (kf=0.15");
plot(n,simulated_CCDF_SNR4,'mo','DisplayName',"simulated (kf=0.2");
plot(n,computed_CCDF_SNR4,'m-','DisplayName',"computed (kf=0.2");
legend('show');
xlabel('SNR value');
ylabel("Coverage Probability");
title("Coverage Probability (SNR) for RF and Thz");
hold off;

