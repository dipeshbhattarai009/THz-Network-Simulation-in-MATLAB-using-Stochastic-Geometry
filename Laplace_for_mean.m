clear;
r=linspace(1,10,1000);
y=exp(-0.1*r).*r.^-2;
plot(r,y);
xlabel('distance (cm)');
ylabel('Power (Watt)');
title('Path Loss')