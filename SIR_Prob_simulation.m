
%Calculate and plot coverage probability distribution in a interference
%limited wireless communication system in THz range for different value of
%beamwidth of antennas

clear;

% defining all constants
lambda=0.001;  % intensity of transmitters /cm^2
side=200;      % side of square area in cm
Area=side^2;      % area where nodes lie in
Pt=1;           %Transmitting power
Kfa=0.1;        %absorbtion coefficient
alpha=2;          %path loss exponent
no_SIR_threshlods=50;  %no of SIR values to compare with
SIR_threshold=linspace(0,30,no_SIR_threshlods);
no_of_iterations=10000;
Gr_max=20;          %major lobe gain of reciever
Gr_min=0.3;         %minor lobe gain of reciever
Gt_max=20;          %major lobe gain of transmitter
Gt_min=0.3;         %minor lobe gain of transmitter
hpbw=[pi,pi/4,pi/6,pi/8];

% getting the array of arrays of comparison values to find the coverage
% probability
SIR_for_statistics=zeros(no_of_iterations,no_SIR_threshlods);

% array of interferences with each element as Interference Power at
% corresponding iteration
Interference_array=zeros(no_of_iterations,1);
Coverage_Probability=zeros(length(hpbw),no_SIR_threshlods);
CDF_interference=zeros(length(hpbw),10000);
% t = clock;  % Get current time as a vector [year, month, day, hour, minute, second]
% seed_value = sum(t);  % Sum the time components to create a seed
% rng(seed_value);  % Set the seed

PDF_interferences = []; % To store PDF values
bin_centers = []; % To store x-axis values

for n=1:length(hpbw)
    for iterations=1:no_of_iterations
        % getting the no. of nodes using poission distribution
        poission_input=Area*lambda;
        no_of_points=poissrnd(poission_input);

        if no_of_points==0
            i=i-1;
            continue;
        else


            % scattering the nodes in the area randomly
            [x_cd,y_cd,distances,new_size]=place_randomly(no_of_points,side/2);


            %applying probabilistic antenna model for reciever
            array1=binornd(1,hpbw(n)/pi,new_size,1);
            Gr_array=array1.*Gr_max;
            Gr_array(find(Gr_array==0))=Gr_min;

            %applying probabilistic antenna model for transmitter
            array2=binornd(1,hpbw(n)/pi,new_size,1);
            Gt_array=array2.*Gt_max;
            Gt_array(find(Gt_array==0))=Gt_min;

            %nearest node faces the transmitter so it has maximum gain
            %always
            Gt_array(1)=Gt_max;   

            % recieved powers
            recieved_powers=Pt*(Gt_array.*Gr_array.*(distances.^(-alpha)).*exp(-Kfa*distances));

            Pr=recieved_powers(1); % Power of desired node

            Interference_power=sum(recieved_powers(2:end)); %Interference Power

            SIR=Pr/Interference_power;  %SIR
            % SIR=pow2db(SIR);

            compare_result=SIR>SIR_threshold;
            SIR_for_statistics(iterations,:)=compare_result;
            Interference_array(iterations)=Interference_power;
            
        end
    end

    % array of no. of cases where SIR is greater than corresponding value in
    % array of SIR thresholds
    no_of_ones = sum(SIR_for_statistics == 1, 1);

    Coverage_Probability(n,:)=no_of_ones/no_of_iterations;

    %Interference in decibels
    Interference_array=pow2db(Interference_array); 

    %finding probability distributions of Interferences
    Interference_values=linspace(min(Interference_array),max(Interference_array),10000);
    CCDF_Interference=zeros(10000,1);

    for i=1:10000
        CCDF_Interference(i)=sum(Interference_array>=Interference_values(i))/no_of_iterations;
    end

    CDF_interference(n,:)=1-CCDF_Interference;

    f(n,:)= diff(CDF_interference(n,:))./diff(Interference_values);

end

%Coverage distribution plot

figure;
grid on;
plot(SIR_threshold,Coverage_Probability(1,:),'ro');
hold on
plot(SIR_threshold,Coverage_Probability(2,:),'bo');
plot(SIR_threshold,Coverage_Probability(3,:),'go');
plot(SIR_threshold,Coverage_Probability(4,:),'mo');
legend("isotropic","hpbw pi/4","hpbw pi/6", "hpbw pi/8");
xlabel("SIR value");
ylabel("Coverage Probability");
title("Coverage Probability distribution for different beam width of antenna");
hold off

%CDF of interference
figure;
plot(Interference_values,CDF_interference(1,:));
hold on
plot(Interference_values,CDF_interference(2,:));
plot(Interference_values,CDF_interference(3,:));
plot(Interference_values,CDF_interference(4,:));
legend("isotropic","hpbw pi/4","hpbw pi/6", "hpbw pi/8");
xlabel("Interference in db");
ylabel("Cumulative Probability");
title("CDF of Interference");
hold off

%PDF of interference
figure;
plot(Interference_values(1:end-1),f(1,:),"b-");
hold on
plot(Interference_values(1:end-1),f(2,:),'r-');
plot(Interference_values(1:end-1),f(3,:),'y-');
plot(Interference_values(1:end-1),f(4,:),'g-');
legend("isotropic","hpbw pi/4","hpbw pi/6", "hpbw pi/8");
xlabel("Interference in db");
ylabel("Probability");
title("PDF of Interferences");
hold off
