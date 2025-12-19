%Coverage Probabililty distribution for different values of molecular
%absorption coefficient;
clear;

% defining all constants
lambda=0.001;    % intensity of transmitters /cm^2
side=200;              % side of square area in cm
Area=side^2;           % area where nodes lie in
Pt=1;                  %Transmitting power
Kfa=[0.05,0.1,0.15,0.2];               %absorbtion coefficient
alpha=2;               %path loss exponent
no_SIR_threshlods=50;  %no of SIR values to compare with
SIR_threshold=linspace(0,30,no_SIR_threshlods);
no_of_iterations=10000;
Gr_max=1;           %major lobe gain of reciever
Gr_min=0.3;         %minor lobe gain of reciever
Gt_max=1;           %major lobe gain of transmitter
Gt_min=0.3;         %minor lobe gain of transmitter
hpbw=pi/4;

% getting the array of arrays of comparison values to find the coverage
% probability
SIR_for_statistics=zeros(no_of_iterations,no_SIR_threshlods);

% array of interferences with each element as Interference Power at
% corresponding iteration
Interference_array=zeros(no_of_iterations,1);


for l=1:length(Kfa)
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
        % ind=find(distances<1);
        % distances(ind)=[];

        %Array to store gain due to each
        array1=binornd(1,hpbw/(2*pi),new_size,1);
        Gr_array=array1.*Gr_max;
        Gr_array(find(Gr_array==0))=Gr_min;



        array2=binornd(1,hpbw/(2*pi),new_size,1);
        Gt_array=array2.*Gt_max;
        Gt_array(find(Gt_array==0))=Gr_min;

        Gt_array(1)=Gt_max;                 % tx-rx are directed to each other

        % recieved powers
        recieved_powers=Pt*(Gt_array.*Gr_array.*(distances.^(-alpha)).*exp(-Kfa(l)*distances));

        Pr=recieved_powers(1); % Power of desired node

        Interference_power=sum(recieved_powers(2:end)); %Interference Power

        SIR=Pr/Interference_power;  %SIR
        % SIR=pow2db(SIR);

        compare_result=SIR>SIR_threshold;
        SIR_for_statistics(iterations,:)=compare_result;
        Interference_array(iterations)=Interference_power;
    end
end
no_of_ones=sum(SIR_for_statistics==1,1);
Coverage_Probability(l,:)=no_of_ones/no_of_iterations;
end

Interference_array=pow2db(Interference_array);

%finding probability distributions of Interferences
Interference_values=linspace(min(Interference_array),max(Interference_array),10000);
CCDF_Interference=zeros(1,10000);

for i=1:10000
    CCDF_Interference(i)=sum(Interference_array>=Interference_values(i))/no_of_iterations;
end

CDF_interference=1-CCDF_Interference;

CDF_SIR=1-Coverage_Probability;

f= diff(CDF_interference)./diff(Interference_values);
f=[f,0];

% for l=1:length(Kfa)
% for m=1:length(SIR_threshold)
%     term1=abs(sqrt(Kfa*SIR_threshold(m)*Interference_values/(alpha*Gt_max*Gr_max)));
%     term2=(alpha.*lambertw(1./term1)./Kfa(l)).^2;
%     Cov=(1-exp(-pi*lambda.*term2)).*f;
%     Cov_comp(l,m)=trapz(Interference_values,Cov);
% end
% end

figure;
% plot(SIR_threshold,Cov_comp(1,:),'r-');
hold on
plot(SIR_threshold,Coverage_Probability(1,:),'ro',"DisplayName","Kfa=0.05");
% plot(SIR_threshold,Cov_comp(2,:),'b-');
plot(SIR_threshold,Coverage_Probability(2,:),'bo',"DisplayName","Kfa=0.1");
% plot(SIR_threshold,Cov_comp(3,:),'y-');
plot(SIR_threshold,Coverage_Probability(3,:),'yo',"DisplayName","Kfa=0.15");
% plot(SIR_threshold,Cov_comp(4,:),'g-');
plot(SIR_threshold,Coverage_Probability(4,:),'go',"DisplayName","Kfa=0.2");
xlabel("SIR value");
legend("show");
ylabel("Coverage probability");
title("Coverage probability distribution for different values of absorption coefficient");
hold off

figure
plot(Interference_values,f,'b');
xlabel("Interference power (dB)");
ylabel("Probability density");
title("Probability density function of interference in THz wireless system");


