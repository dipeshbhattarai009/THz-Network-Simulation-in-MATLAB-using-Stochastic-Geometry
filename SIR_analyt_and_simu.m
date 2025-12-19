clear;

% defining all constants
lambda=[0.001,0.002,0.003,0.004];      % intensity of transmitters /cm^2
side=200;              % side of square area in cm
Area=side^2;           % area where nodes lie in
Pt=1;                  %Transmitting power
Kfa=0.1;               %absorbtion coefficient
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


for l=1:length(lambda)
for iterations=1:no_of_iterations
    % getting the no. of nodes using poission distribution
    poission_input=Area*lambda(l);
    no_of_points=poissrnd(poission_input);

    if no_of_points==0
        i=i-1;
        continue;
    else


        % scattering the nodes in the area randomly
        [x_cd,y_cd,distances,new_size]=place_randomly(no_of_points,side/2);
        % ind=find(distances<1);
        % distances(ind)=[];
        % new_size=length(distances);
        % x_cd(ind)=[];
        % y_cd(ind)=[];

        % angles=atan2(y_cd,x_cd); %angle made by each transmitter
        %
        % % finding the points that are within -hpbw/2 to hpbw/2 range of x axis
        % indices_los=find((angles>=-hpbw)&(angles<=hpbw));
        % points_in_los=[x_cd(indices_los),y_cd(indices_los)];
        %
        % Gr_array=ones(new_size,1)*Gr_min;
        % Gr_array(indices_los)=Gr_max;   %maximum gain for transmitters which are in major lobe of reciever

        array1=binornd(1,hpbw/(2*pi),new_size,1);
        Gr_array=array1.*Gr_max;
        Gr_array(find(Gr_array==0))=Gr_min;

        % x_fr_nodes=-x_cd;   %coordinate of origin seen by tx nodes
        % y_fr_nodes=-y_cd;
        %
        % rand_coeff_for_phi=2*rand(new_size,1)-1; %random angles for rotating the beams randomly
        % phi=rand_coeff_for_phi*pi;
        %
        % x_fr_nodes2=x_fr_nodes.*cos(phi)-y_fr_nodes.*sin(phi); %rotating the beams randomly
        % y_fr_nodes2=x_fr_nodes.*sin(phi)+y_fr_nodes.*cos(phi);
        % orientations=atan2(y_fr_nodes2,x_fr_nodes2);
        %
        % %finding the nodes which find origin in thier major lobe
        % tx_max_indices=find((orientations>=-hpbw)&(orientations<=hpbw));
        %
        % Gt_array=ones(new_size,1)*Gt_min;  %declaring array of transmitting gains
        % Gt_array(tx_max_indices)=Gt_max;    %gain of tx is maximum for those which find origin at major lobe


        array2=binornd(1,hpbw/(2*pi),new_size,1);
        Gt_array=array2.*Gt_max;
        Gt_array(find(Gt_array==0))=Gr_min;

        Gt_array(1)=Gt_max;                 % tx-rx are directed to each other

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

for l=1:length(lambda)
for m=1:length(SIR_threshold)
    term1=abs(sqrt(Kfa*SIR_threshold(m)*Interference_values/(alpha*Gt_max*Gr_max)));
    term2=(alpha.*lambertw(1./term1)./Kfa).^2;
    Cov=(1-exp(-pi*lambda(l).*term2)).*f;
    Cov_comp(l,m)=trapz(Interference_values,Cov);
end
end

figure;
% plot(SIR_threshold,Cov_comp(1,:),'r-');
hold on
plot(SIR_threshold,Coverage_Probability(1,:),'ro',"DisplayName","lambda=0.001");
% plot(SIR_threshold,Cov_comp(2,:),'b-');
plot(SIR_threshold,Coverage_Probability(2,:),'bo',"DisplayName","lambda=0.002");
% plot(SIR_threshold,Cov_comp(3,:),'y-');
plot(SIR_threshold,Coverage_Probability(3,:),'yo',"DisplayName","lambda=0.003");
% plot(SIR_threshold,Cov_comp(4,:),'g-');
plot(SIR_threshold,Coverage_Probability(4,:),'go',"DisplayName","lambda=0.004");
xlabel("SIR value");
legend("show");
ylabel("Coverage probability")
title("Coverage probability distribution for different node densities (lambda)");
hold off



