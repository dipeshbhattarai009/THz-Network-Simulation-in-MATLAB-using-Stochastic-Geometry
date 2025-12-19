
% Plotting and Comparing mean of SNR wrt using simulation and
% analytical derivation

clear;
lambda=0.001; % node density
alpha=2 ; % path loss exponent


% parameters to define area of simulation
xmin=-100; 
xmax=100;
ymin=-100;
ymax=100;

kf=0.1;  % molecular absorption coefficient (cm^-1)
No=linspace(0.0001,0.001,100); % Noise Power
iteration=10000;

%finding means of SNR using simulation
for j=1:length(No)

    for i=1:iteration
        N=poissrnd(lambda.*(xmax-xmin).*(ymax-ymin)); % number of nodes
        X1=unifrnd(xmin,xmax,1,N); % x coordinates
        Y1=unifrnd(xmin,xmax,1,N); % y coordinates

        %distance from reciever (at origin)
        distances=sqrt(X1.^2+Y1.^2); 
        distances=sort(distances);

        %remove nodes with distance smaller than unity
        ind=find(distances<1); 
        distances(ind)=[];

        %recieved power without molecular absorption (RF)
        Recieved_pr1(i)=distances(1)^-alpha/No(j); 

        %recieved power with molecular absorption (Thz)
        Recieved_pr2(i)=distances(1)^-alpha*exp(-kf*distances(1))/No(j); 
    end

    %mean SNR
    mean_SNR1(j)=mean(Recieved_pr1); 
    mean_SNR2(j)=mean(Recieved_pr2);
end


%Computation of mean of SNR with analytical derived formulas
for j=1:length(No)
    f1= @(r) (r.^-alpha).*(2*pi*lambda*r.*(exp(-lambda*pi*r.^2)))/No(j);
    f2 = @(r) (r.^-alpha).*(exp(-kf*r)).*(2*pi*lambda*r.*(exp(-lambda*pi*r.^2)))/No(j);

    comp_mean_SNR1(j) =  integral(f1,1,inf);
    comp_mean_SNR2(j) =  integral(f2,1,inf);
end

%Plots of variation of mean SNR with respect to noise power for RF and Thz 

figure;
plot(No,mean_SNR1,'ro','DisplayName','RF simulated mean SNR');
hold on
plot(No,mean_SNR2,'bo','DisplayName','Thz simulated mean SNR');

plot(No,comp_mean_SNR1,'-r','DisplayName','RF computed mean SNR');
plot(No,comp_mean_SNR2,'-b','DisplayName','Thz computed mean SNR');
xlabel('Noise power (Watt)');
ylabel('mean SNR');
legend('show');
title('Mean SNR vs Noise power for RF and Thz');
hold off

