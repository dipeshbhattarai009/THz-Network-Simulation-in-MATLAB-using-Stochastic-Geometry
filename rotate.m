clear;

% defining all constants
lambda=0.15;
Area=100;
Pt=1;
Kfa=0.1;
alpha=2;
SIR_values=linspace(0,5,15);
no_of_iterations=10000;
Gr_max=100;
Gr_min=2;
Gt=1;
% getting the array of arrays of comparison values to find the probability
% of SIR
SIR_for_statistics=zeros(no_of_iterations,15);

for th=1:no_of_iterations
    % getting the no. of nodes using poission distribution
    poission_input=Area*lambda;
    no_of_points=poissrnd(poission_input);
    Gr_array=ones(no_of_points,1)*Gr_min;

    if no_of_points==0
        i=i-1;
        continue;
    else


        % scattering the nodes in the area randomly
        random_coefficients_for_xcoordinates=2*rand(no_of_points,1)-1;
        random_coefficients_for_ycoordinates=2*rand(no_of_points,1)-1;
        x_coordinates=random_coefficients_for_xcoordinates.*sqrt(Area);
        y_coordinates=random_coefficients_for_ycoordinates.*sqrt(Area);

        % distances
        distances=sqrt(x_coordinates.^2+y_coordinates.^2);
        [distances2,index]=sort(distances);

        %finding nearest point
        nearest_point=[x_coordinates(index(1)),y_coordinates(index(1))];

        %finding angle to rotate so that the nearest point lies in x-axis
        angle_rotate=atan2(nearest_point(2),nearest_point(1));
        %rotating the points
        x_coordinates2 = x_coordinates*cos(angle_rotate)+y_coordinates*sin(angle_rotate);
        y_coordinates2 = -x_coordinates*(sin(angle_rotate))+y_coordinates*cos(angle_rotate);

        x_coordinates2=x_coordinates2(index);
        y_coordinates2=y_coordinates2(index);
        
        %finding the points that are within -pi/4 to pi/4 range of x axis
        angles=atan2(y_coordinates2,x_coordinates2);
        indices_los=find((angles>=-(pi/4))&(angles<=(pi/4)));
        
        points_in_los=[x_coordinates2(indices_los),y_coordinates2(indices_los)];
        Gr_array(indices_los)=Gr_max;


        % recieved powers
        recieved_powers=Pt*Gt*(Gr_array.*(distances2.^(-alpha)).*exp(-Kfa*distances2));

        Pr=recieved_powers(1);
        Interference_power=sum(recieved_powers(2:end));

        SIR=Pr/Interference_power;
       
        compare_result=SIR>SIR_values;
        SIR_for_statistics(th,:)=compare_result;
    end
end

% finding the no. of conditions in which certain SIR is good
no_of_ones = sum(SIR_for_statistics == 1, 1);

Probability_of_SIR_values=no_of_ones/no_of_iterations;

%plotting the nodes

% figure;
% scatter(x_coordinates,y_coordinates,'bo');
% hold on;
% plot(0,0,'rx');
% title('location of nodes');
% axis equal;

x = 0:0.1:10;
y1 = tan(-pi/4) * x;
y2 = tan(pi/4) * x;

figure;
scatter(x_coordinates2,y_coordinates2,'bo');
hold on;
plot(0,0,'r+');
plot(x, y1, 'g-', 'LineWidth', 1);
plot(x, y2, 'g-', 'LineWidth', 1);
title('location of nodes');
axis equal;

figure;
plot(SIR_values,Probability_of_SIR_values);
xlabel('SIR');
ylabel('Probability of SIR');
title('SIR values vs Probability of actual SIR being greater than them');

