clear;
close all;


% defining all constants
lambda=0.001; %no of base stations 
side=500;      %side of square area taken
Area=side^2;   %area of the square
Pt=1;          %transmission power
Kfa=0.1;       %absorbtion coefficients
alpha=3;       %path loss exponent
SIR_values=linspace(0,30,50);  
no_of_iterations=10000;
Gr_max=30;
Gr_min=1.3;
Gt_min=1.3;
Gt_max=30;

% getting the array of arrays of comparison values to find the probability
% of SIR

SIR_for_statistics=zeros(no_of_iterations,50);
SIR_for_statistics_isotropic=zeros(no_of_iterations,50);

Interference_array=zeros(no_of_iterations,1);
Interference_array_isotropic=zeros(no_of_iterations,1);

for t=1:no_of_iterations
    % getting the no. of nodes using poission distribution
    poission_input=Area*lambda;
    no_of_points=poissrnd(poission_input);
    
    % tx_prob=binornd(1,0.5,1,no_of_points);
    % 
    % tx_max_indices=find(tx_prob==1);
     

    if no_of_points==0
        i=i-1;
        continue;
    else


        % scattering the nodes in the area randomly
        [x_cd,y_cd,distances,new_size]=place_randomly(no_of_points,side);
        no_of_points=new_size;

        Gr_array=ones(no_of_points,1)*Gr_min;
        Gt_array=ones(no_of_points,1)*Gt_min;

        x_fr_nodes=-x_cd;
        y_fr_nodes=-y_cd;

        rand_coeff_for_phi=rand(no_of_points,1);
        phi=rand_coeff_for_phi*2*pi;

        x_fr_nodes2=x_fr_nodes.*cos(phi)+y_fr_nodes.*sin(phi);
        y_fr_nodes2=-x_fr_nodes.*sin(phi)+y_fr_nodes.*cos(phi);

        orientations=atan2(y_fr_nodes2,x_fr_nodes2);
        tx_max_indices=find((orientations>=-(pi/6))&(orientations<=pi/6));

        Gt_array(tx_max_indices)=Gt_max; %finding the antennas on line of sight;


        Gt_array(1)=Gt_max;
        Gr_array(1)=Gr_max;

        %finding the points that are within -pi/4 to pi/4 range of x axis
        angles=atan2(y_cd,x_cd);
        indices_los=find((angles>=-(pi/6))&(angles<=(pi/6)));
        
        points_in_los=[x_cd(indices_los),y_cd(indices_los)];
        Gr_array(indices_los)=Gr_max;


        % recieved powers
        recieved_powers=Pt*Gt_array.*(Gr_array.*(distances.^(-alpha)).*exp(-Kfa*distances));
        recieved_powers_isotropic=distances.^(-alpha).*exp(-Kfa*distances);

        Pr=recieved_powers(1);
        Pr_isotropic=recieved_powers_isotropic(1);
        Interference_power=sum(recieved_powers(2:end));
        Interference_power_isotropic=sum(recieved_powers_isotropic(2:end));

        SIR=Pr/Interference_power;
        SIR=10*log10(SIR+0.0001);
        SIR_isotropic=Pr_isotropic/Interference_power_isotropic;
        SIR_isotropic=10*log10(SIR_isotropic+0.0001);
       
        compare_result=SIR>SIR_values;
        compare_result_isotropic=SIR_isotropic>SIR_values;

        SIR_for_statistics(t,:)=compare_result;
        SIR_for_statistics_isotropic(t,:)=compare_result_isotropic;

        Interference_array(t)=Interference_power;
        Interference_array_isotropic(t)=Interference_power_isotropic;
    end
end

% finding the no. of conditions in which certain SIR is good
no_of_ones = sum(SIR_for_statistics == 1, 1);
no_of_ones_iso=sum(SIR_for_statistics_isotropic==1,1);

Probability_of_SIR_values=no_of_ones/no_of_iterations;
Probability_of_SIR_values_iso=no_of_ones_iso/no_of_iterations;

%Pdf_SIR=gradient((1-Probability_of_SIR_values))./gradient(SIR_values);
%Pdf_SIR_isotropic=gradient((1-Probability_of_SIR_values_iso))./gradient(SIR_values);

Interference_array=10*log10(Interference_array+0.0001);
Interference_array_isotropic=10*log10(Interference_array+0.0001);

Interference_values=linspace(min(Interference_array),max(Interference_array),100);

Prob_Interference=zeros(100,1);
Prob_Interference_isotropic=zeros(100,1);

for i=1:100
    Prob_Interference(i)=sum(Interference_array>=Interference_values(i))/no_of_iterations;
end

for i=1:100
    Prob_Interference_isotropic(i)=sum(Interference_array_isotropic>=Interference_values(i))/no_of_iterations;
end

Pdf_direct_interference=gradient(1-Prob_Interference)./gradient(Interference_values);
Pdf_isotropic_interference=gradient(1-Prob_Interference_isotropic)./gradient(Interference_values);



%Plottings

figure;
scatter(x_cd,y_cd,'bo');
hold on;
plot(0,0,'r+');

title('location of nodes');
axis equal;
hold off;

figure;
plot(SIR_values,Probability_of_SIR_values);
hold on
plot(SIR_values,Probability_of_SIR_values_iso,'--');
xlabel('SIR');
ylabel('Probability of SIR');
title('SIR values vs Probability of actual SIR being greater than them');
legend('Directional','isotropic','Location','southwest');
hold off

figure;
plot(Interference_values,Prob_Interference);
hold on
plot(Interference_values,Prob_Interference_isotropic,'.');
xlabel('Interference');
ylabel('Probability of Interference');
title('Interference Distribution');
legend('Directional','isotropic');
hold off

figure;
plot(Interference_values,Pdf_direct_interference,'b-');
hold on
plot(Interference_values,Pdf_isotropic_interference,'o-');
xlabel('Interference');
ylabel('Probability');
title('PDF of Interference');
legend('Directional','Isotropic');
hold off