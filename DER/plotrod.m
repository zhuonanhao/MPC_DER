function plotrod(x, d1, d2, m1, m2, ctime, saveImage, imageDirectory)
global dt
nv = (length(x)+1)/4;
x1 = x(1:4:end);
x2 = x(2:4:end);
x3 = x(3:4:end);
L = sum(sqrt( (x1(2:end) - x1(1:end-1)).^2 + ...
    (x2(2:end) - x2(1:end-1)).^2 + ...
    (x3(2:end) - x3(1:end-1)).^2));
d1 = 0.1*L * d1;
d2 = 0.1*L * d2;
m1 = 0.1*L * m1;
m2 = 0.1*L * m2;

h1 = figure(1);
% set(h1, 'visible', 'off');
clf()
plot3(x1,x2,x3, 'k-', 'LineWidth',2);
hold on
plot3(x1(1),x2(1),x3(1), 'k^', 'LineWidth',2);
plot3(x1(end),x2(end),x3(end), 'k^', 'LineWidth',2);

middleIndex = floor(nv/2);
plot3(x1(middleIndex),x2(middleIndex),x3(middleIndex), 'ro', 'LineWidth',2);
% for c=1:nv-1
%     xa = x(4*c-3:4*c-1);
%     xb = x(4*c+1:4*c+3);
%     xp = (xa+xb)/2;
%     p1 = plot3( [xp(1), xp(1) + d1(c,1)], [xp(2), xp(2) + d1(c,2)], ...
%         [xp(3), xp(3) + d1(c,3)], 'b--', 'LineWidth', 2);
%     p2 = plot3( [xp(1), xp(1) + d2(c,1)], [xp(2), xp(2) + d2(c,2)], ...
%         [xp(3), xp(3) + d2(c,3)], 'c--', 'LineWidth', 2);
%     p3 = plot3( [xp(1), xp(1) + m1(c,1)], [xp(2), xp(2) + m1(c,2)], ...
%         [xp(3), xp(3) + m1(c,3)], 'r-');
%     p4 = plot3( [xp(1), xp(1) + m2(c,1)], [xp(2), xp(2) + m2(c,2)], ...
%         [xp(3), xp(3) + m2(c,3)], 'g-');
% end
hold off
% legend([p1,p2,p3,p4], 'd_1','d_2','m_1','m_2');
title(num2str(ctime, 't=%.2f'));
% axis equal
view(0,0);
xlabel('x');
ylabel('y');
zlabel('z');
xlim([-0.1, 0.1])
ylim([-0.05, 0.05])
zlim([-0.1, 0.1])
box on
daspect([1 1 1])
axis off
drawnow
if (saveImage~=0)
    saveas(h1, num2str(ctime/dt, [imageDirectory, '/t=%09.0f.png']));
end
