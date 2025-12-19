r=linspace(0,10);
lambda=[0.1,0.5,0.3,5];
figure;
n=1:10;
for n=1:length(n)
    for i=1:length(lambda)
        l=lambda(i);
        y=(exp(-l*pi*r.^2).*(l*pi*r.^2).^n)/factorial(n);
        plot(r,y);
        hold on;

        label=sprintf("\\lambda=%g",lambda);
        legend(label);
    end
end