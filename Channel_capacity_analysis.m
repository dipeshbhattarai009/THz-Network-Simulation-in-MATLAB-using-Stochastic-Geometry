% Plotting and comparing Channel Capacity vs node density using simulation and analytical
% derivation

clear;
lambda=0.001; %node density
alpha=2; %path loss exponent

%area for simulation
xmin=-100; 
xmax=100;
ymin=-100;
ymax=100;

kf1=0.1; % absorption coefficient
kf2=0.15;
kf3=0.2;

No=linspace(0.0001,0.001,100); %Noise Power
iteration=1000;

B=10^6; %Channel Bandwidth

%Simulation  
for j=1:length(No)

    for i=1:iteration
        N=poissrnd(lambda.*(xmax-xmin).*(ymax-ymin)); %Number of nodes
        X1=unifrnd(xmin,xmax,1,N); % X coordinates
        Y1=unifrnd(xmin,xmax,1,N); % Y coordinates

        %distance from the reciever (at origin)
        distances=sqrt(X1.^2+Y1.^2); 
        distances=sort(distances);

        %remove distances smaller than unity
        ind=find(distances<1); 
        distances(ind)=[];

        %Recieved Power without molecular absorption
        Recieved_snr_rf=distances(1)^-alpha/No(j); 

        %Recieved power with molecular absorption
        Recieved_snr_thz1=distances(1)^-alpha*exp(-kf1*distances(1))/No(j);
        Recieved_snr_thz2=distances(1)^-alpha*exp(-kf2*distances(1))/No(j);
        Recieved_snr_thz3=distances(1)^-alpha*exp(-kf3*distances(1))/No(j);
        %Channel capacity for RF and Thz
        Capacity_rf(i)= B*log2(1+Recieved_snr_rf);
        Capacity_thz1(i)=B*log2(1+Recieved_snr_thz1);
        Capacity_thz2(i)=B*log2(1+Recieved_snr_thz2);
        Capacity_thz3(i)=B*log2(1+Recieved_snr_thz3);

    end
    %Mean channel capacities
    mean_cap_rf(j)=mean(Capacity_rf);
    mean_cap_thz1(j)=mean(Capacity_thz1);
    mean_cap_thz2(j)=mean(Capacity_thz2);
    mean_cap_thz3(j)=mean(Capacity_thz3);
end

%Computation of SNR  with analytically derived formulas
for j=1:length(No)
    f1 = @(r) (r.^-alpha).*(2*pi*lambda*r.*(exp(-lambda*pi*r.^2)))/No(j);
    f2 = @(r) (r.^-alpha).*(exp(-kf1*r)).*(2*pi*lambda*r.*(exp(-lambda*pi*r.^2)))/No(j);
    f3 = @(r) (r.^-alpha).*(exp(-kf2*r)).*(2*pi*lambda*r.*(exp(-lambda*pi*r.^2)))/No(j);
    f4 = @(r) (r.^-alpha).*(exp(-kf3*r)).*(2*pi*lambda*r.*(exp(-lambda*pi*r.^2)))/No(j);

    comp_mean_SNR1(j) =  integral(f1,1,inf);
    comp_mean_SNR2(j) =  integral(f2,1,inf);
    comp_mean_SNR3(j) =  integral(f3,1,inf);
    comp_mean_SNR4(j) =  integral(f4,1,inf);
end

%computation of capacities using mean SNR values
cap_mean_rf_comp=B*log2(1+comp_mean_SNR1);
cap_mean_thz_comp1=B*log2(1+comp_mean_SNR2);
cap_mean_thz_comp2=B*log2(1+comp_mean_SNR3);
cap_mean_thz_comp3=B*log2(1+comp_mean_SNR4);


%Plots of mean capacity variation with noise power for RF and THz system
%with different molecular absorption coefficients for noise limited case
figure;
plot(No,mean_cap_rf,'ro','DisplayName','RF simulated capacity');
hold on
plot(No,cap_mean_rf_comp,"r-",'DisplayName','RF computed mean capacity');

plot(No,mean_cap_thz1,'bo','DisplayName','Thz simulated mean capacity (Kfa=0.1');
plot(No,cap_mean_thz_comp1,"b-",'DisplayName','Thz computed mean capacity (Kfa=0.1)');

plot(No,mean_cap_thz2,'go','DisplayName','Thz simulated mean capacity (Kfa=0.15');
plot(No,cap_mean_thz_comp2,"g-",'DisplayName','Thz computed mean capacity (Kfa=0.15)');

plot(No,mean_cap_thz3,'mo','DisplayName','Thz simulated mean capacity (Kfa=0.2');
plot(No,cap_mean_thz_comp3,"m-",'DisplayName','Thz computed mean capacity (Kfa=0.2)');

title('Mean channel capacity variation for RF and THz wrt noise power');
xlabel('Noise power (No)');
ylabel('Mean capacity (bps)');
legend('show');
hold off;
