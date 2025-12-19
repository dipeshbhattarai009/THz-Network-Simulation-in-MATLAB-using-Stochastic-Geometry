function [x_cd,y_cd,distances,new_size]=place_randomly(no_of_points,max_size)

% x and y coordinates using uniform random distribution
x_coordinates=unifrnd(-max_size/2,max_size/2,no_of_points,1);
y_coordinates=unifrnd(-max_size/2,max_size/2,no_of_points,1);

%implementing blockage using boolean blockage model
beta=0.002;
block_prob=binornd(1,1-exp(-beta*sqrt(x_coordinates.^2+y_coordinates.^2)),no_of_points,1);
block_indice=find(block_prob==1);
new_size=no_of_points-length(block_indice);
x_coordinates(block_indice)=[];
y_coordinates(block_indice)=[];

%sorting points according to distance
[distances,index]=sort(sqrt(x_coordinates.^2+y_coordinates.^2));
ind=find(distances<1);
x_cd=x_coordinates(index);
y_cd=y_coordinates(index);
x_cd(ind)=[];
y_cd(ind)=[];
end   