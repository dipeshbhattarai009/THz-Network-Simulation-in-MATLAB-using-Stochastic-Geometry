%CDF of Poission point process with simulation and analytical

clear;
lambda=0.001;
xmin=-100;
xmax=100;
ymin=-100;
ymax=100;
iteration=3000;
r=linspace(0,250,100);
for i=1:iteration
N=poissrnd(lambda.*(xmax-xmin).*(ymax-ymin));
X1=unifrnd(xmin,xmax,1,N);
Y1=unifrnd(xmin,xmax,1,N);

distances=sqrt(X1.^2+Y1.^2);
distances=sort(distances);
F(i,:)=distances(1)<r;

end

CDF=mean(F);

plot(r,CDF,'bo');
hold on

plot(r,1-exp(-lambda.*pi.*r.^2),'-r')




% compare_matrix=zeros(no_of_pts,no_of_pts);
% for i=1:no_of_pts
%     compare_matrix(i,:)=distances>distances(i);
% end
% mean_array=mean(compare_matrix);
% 
% x=linspace(0,side,no_of_pts);
% CDF_calculated=1-exp(-lambda*pi*x.^2);
% 
% figure;
% plot(x_coord,y_coord,'b.');
% 
% figure;
% plot(x,CDF_calculated,'-');
% hold on
% plot(x,mean_array,'ro');