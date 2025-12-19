clear;
lambda=0.001;
side=1000;
area=side^2;
alpha=4;
no_of_iterations=10000;
interferences=zeros(no_of_iterations,1);
for i=1:no_of_iterations
    poission_input=area*lambda;
    no_of_points=poissrnd(poission_input);
    [x_cd,y_cd,distances,new_size]=place_randomly(no_of_points,side/2);
    rel_distances_sqrd=(distances/min(distances)).^2;
    power_recieved=rel_distances_sqrd.^(-alpha);
    interferences(i)=sum(power_recieved);
end
average_inteference=mean(interferences);
disp(average_inteference);