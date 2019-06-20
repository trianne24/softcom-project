function drawCircle(x, y, r)

theta = 0:pi/360:2*pi;
v1 = x + r*sin(theta);
v2 = y + r*cos(theta);

plot(v1, v2, ':');




