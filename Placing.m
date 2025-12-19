
no_of_points=30;
max_size=100;

random_coefficients_for_x=2*rand(no_of_points,1)-1;
random_coefficients_for_y=2*rand(no_of_points,1)-1;

x_coordinates=random_coefficients_for_x.*max_size;
y_coordinates=random_coefficients_for_y.*max_size;

beta=0.002;
block_prob=binornd(1,1-exp(-beta*sqrt(x_coordinates.^2+y_coordinates.^2)),no_of_points,1);

block_indice=find(block_prob==1);

new_size=no_of_points-length(block_indice);

x_coordinates(block_indice)=[];
y_coordinates(block_indice)=[];

[distances,index]=sort(sqrt(x_coordinates.^2+y_coordinates.^2));

% x_cd=x_coordinates(index);
% y_cd=y_coordinates(index);
% 
% theta=atan2(y_cd(1),x_cd(1));
% 
% x_coordinates2=x_cd.*cos(theta(1))+y_cd.*sin(theta(1));
% y_coordinates2=-x_cd.*sin(theta(1))+y_cd.*cos(theta(1));

scatter(x_coordinates,y_coordinates,'bo');
hold on;
plot(0,0,'+');
