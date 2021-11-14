function retval = ShiftRot2d(x,y,Dx,Dy,Dphi)
% x,y - vectors to transform (can be matrices)
% Dx,Dy - shift vectors (can be matrices)
% Dphi - rotation angles (can be a matrix)

[imax jmax] = size(x);
% ~= is the same as !=
if [imax jmax] ~= size(y)  | [imax jmax] ~= size(Dphi) | ...
   [imax jmax] ~= size(Dx) | [imax jmax] ~= size(Dy)
  error('stop ShiftRot2d: sizes of input do not match');
end

xs =  x.*cos(-Dphi) + y.*sin(-Dphi);
ys = -x.*sin(-Dphi) + y.*cos(-Dphi);

xp = xs + Dx;
yp = ys + Dy;

retval = [ xp, yp ];

