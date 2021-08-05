function[] = err_plot(t,x,std,color)
t = t';
x = x';
std = std';
if color == "blue"
    cl = 'b';
    cp = [51 153 255]./255;
elseif color == "red"
    cl = 'r';
    cp = [205 133 63]./255;
elseif color == "green"
    cl = 'g';
    cp = [34 139 34]./255;
end
x_min = x - std;
x_max = x + std;
hold on
plot(t,x,cl,'Linewidth',2)

patch_x = fill([t,fliplr(t)],[x_min,fliplr(x_max)],cp);
% plot(t,x_min,cl,'Linewidth',2)
% plot(t,x_max,cl,'Linewidth',2)
set(patch_x, 'edgecolor', 'none');
set(patch_x, 'FaceAlpha', 0.5);
hold off
end