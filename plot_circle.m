function plot_circle(radius)
x = 0:0.1:radius*cos(pi/6);
y1 = tan(11*pi/6)*x;
y2 = tan(pi/6) *x;
angles=linspace(0,2*pi,100);
circle_x=radius*cos(angles);
circle_y=radius*sin(angles);

plot(circle_x,circle_y,'r-');
plot(x, y1, 'g-', 'LineWidth', 1);
plot(x, y2, 'g-', 'LineWidth', 1);
end