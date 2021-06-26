function[R]=rotation(alpha,achse)
% 3-D rotation, alpha is the ratationsangle in radiant
% aches are axis, aound which the rotation happens
% the output R is the matrix for the rotation. 
if achse == 'x'
    R = [1,0,0;0,cosd(alpha),sind(alpha);0,-sind(alpha),cosd(alpha)];
elseif achse == 'y'
    R = [cosd(alpha),0,-sind(alpha);0,1,0;sind(alpha),0,cosd(alpha)];
elseif achse == 'z'
    R = [cosd(alpha),sind(alpha),0;-sind(alpha),cosd(alpha),0;0,0,1];
else
    error('Fehler');
end
end