% Plotting and Comparing mean of Interference using simulation and
% analytical derivation

clear;

lambda = linspace(0.0005,0.005,10); %node densities
side = 200;
Area = side.^2;
Kfa =0.1; %molecular absorption coefficient
alpha =2; %path loss exponent
no_of_iterations =10000;

Pt=1;
Gr_max=1;          %major lobe gain of reciever
Gr_min=0.3;         %minor lobe gain of reciever
Gt_max=1;          %major lobe gain of transmitter
Gt_min=0.3;         %minor lobe gain of transmitter
hpbw= pi/4;         %Half power beam width of antenna


for n=1:length(lambda) 
    for iterations=1:no_of_iterations
        poisson_input=Area*lambda(n);
        no_of_points=poissrnd(poisson_input); %no. of points

        if no_of_points==0
            i=i-1;
            continue;
        else

            [x_cd,y_cd,distances,new_size]=place_randomly(no_of_points,side); %placement of nodes

            %remove nodes closer than distance of 1 as path loss model isn't accurate for them
            ind=find(distances<1);
            distances(ind)=[];
            new_size=length(distances);
            x_cd(ind)=[];
            y_cd(ind)=[];

              % angles at which transmitters are
            % rx_max_angles=atan2(y_cd,x_cd);
            % indices_los=find((rx_max_angles>=-hpbw)&(rx_max_angles<=hpbw));
            % points_in_los=[x_cd(indices_los),y_cd(indices_los)];

            array1=binornd(1,0.25,new_size,1);
            Gr_array=array1.*Gr_max;
            Gr_array(find(Gr_array==0))=Gr_min;

            % % location of reciever wrt transmitter (at origin)
            % x_fr_nodes=-x_cd;
            % y_fr_nodes=-y_cd;
            % 
            % tx_angles=unifrnd(-pi,pi,new_size,1); %angle of origin wrt transmitter
            % 
            % %rotating the location of reciever with random angles for each
            % %transmitter
            % x_fr_nodes2=x_fr_nodes.*cos(tx_angles)-y_fr_nodes.*sin(tx_angles);
            % y_fr_nodes2=x_fr_nodes.*sin(tx_angles)+y_fr_nodes.*cos(tx_angles);
            % orientations=atan2(y_fr_nodes2,x_fr_nodes2);
            % 
            % %nodes which are directed towards reciever
            % tx_max_indices=find((orientations>=-hpbw)&(orientations<=hpbw));
            % Gt_array=ones(new_size,1)*Gt_min;
            % Gt_array(tx_max_indices)=Gt_max;

            array2=binornd(1,0.25,new_size,1);
            Gt_array=array2.*Gt_max;
            Gt_array(find(Gt_array==0))=Gr_min;

            % Gt_array(1)=Gt_max;
            % Gr_array(1)=Gr_max;

            %Recieved powers at reciever
            recieved_powers=Pt*Gt_array.*Gr_array.*(distances.^(-alpha)).*exp(-Kfa*distances);
            recieved_powers=sort(recieved_powers,'descend');

            Interference_power=sum(recieved_powers(1:end));
            Interference_array(iterations,n)=Interference_power;
            mean1=Gt_min*0.75+Gt_max*0.25;
            mean2=Gr_min*0.75+Gr_max*0.25;
        end

    end
end

mean_interference=mean(Interference_array); %Array of mean interference for each density

%computation of mean Interference for each node density
for n=1:length(lambda)
    f1 = @(x) mean1*mean2*2*pi*lambda(n)*x.^(1 - alpha).*exp(-Kfa*x);
    I(n)=  integral(f1,1,inf); %mean interference
end

figure;
plot(lambda,I,'b-');
hold on
plot(lambda,mean_interference,'ro');
title("Mean Interference vs node density");
xlabel("node density");
ylabel("Mean interference");
legend("Analytical computation","simulation");
hold off
