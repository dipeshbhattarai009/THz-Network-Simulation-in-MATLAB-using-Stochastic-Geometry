% Plotting and comparing mean inteference with simulation and computation
% method wrt node densities for isotropic antennas

clear;

lambda = linspace(0.0005,0.001,10); %node densities
side = 200;
Area = side.^2; %area of simulation
Pt = 1; %Transmitted Power
Gt = 1; %Transmitter gain
Gr = 1; %Reciever gain
Kfa =0.1; % molecular absorption coefficient (cm^-1)
alpha = 2; %path loss exponent
no_of_iterations =10000;
Total_interference = zeros(no_of_iterations,1);
Average_Interference=zeros(length(alpha),1);

%simulating and finding mean interference power
for j=1:length(lambda)

    for i = 1:no_of_iterations
        poisson_input = lambda(j) * Area;
        
        %Number of nodes
        no_of_points = poissrnd(poisson_input); 

        %placing nodes using poission point process
        [x_cd,y_cd,distances,new_size]=place_randomly(no_of_points,side);

        % distances = sqrt(x_cd.^2 + y_cd.^2);
        ind=find(distances<1); % remove nodes closer than unit distance
        distances(ind)=[];

        %Recieved power and total interference
        received_power=Pt*Gt*Gr*distances.^(-alpha).*exp(-Kfa*distances);
        Total_interference(i) = sum(received_power(1:end));
    end
    Average_Interference(j) = mean(Total_interference); %Mean interference for each node density
end

compute_integral=zeros(length(lambda),1);

%computation of mean interference using analytically derived formulas
for m=1:length(lambda)
    f = @(x) 2*pi*lambda(m)*x.^(1 - alpha).*exp(-Kfa*x);
    compute_integral(m) =  integral(f,1,inf);
end

%Plots of mean interferences calculated through simulation and analytical
%derivaiton
figure;
plot(lambda,Average_Interference,'r+','DisplayName','Simulated Interference');
hold on
plot(lambda,compute_integral,'b-','DisplayName','Computed Interference');
title("Mean interference for different node densities");
xlabel('Node density (lambda)');
ylabel('Mean Interference (Watt)');
legend('show');

