%% out = A93__94_test_data (5, 7, 17, 25);
%% for ii = 1 : size (out, 1)
%%   figure
%%   plot3 (out{ii, 1, 1}(:), out{ii, 1, 2}(:), out{ii, 1, 3}(:), 'x')
%%   hold on
%%   mesh (out{ii, 2, 1}, out{ii, 2, 2}, out{ii, 2, 3})
%% end

function out = A93_94_95_test_data (nin, min, nplot, mplot)

  [uin, vin]     = ndgrid (linspace (.1, 1, nin), linspace (.1, 1, min));
  [uplot, vplot] = ndgrid (linspace (.1, 1, nplot), linspace (0, 1, mplot));
  [out{1, 1, 1}, out{1, 1, 2}, out{1, 1, 3}] = cone (uin, vin);
  [out{1, 2, 1}, out{1, 2, 2}, out{1, 2, 3}] = cone (uplot, vplot);

  [uin, vin]     = ndgrid (linspace (.1, .9, nin), linspace (-.7, .7, min));
  [uplot, vplot] = ndgrid (linspace (.1, .9, nplot), linspace (-.7, .7, mplot));
  [out{2, 1, 1}, out{2, 1, 2}, out{2, 1, 3}] = ellipsoid (uin, vin);
  [out{2, 2, 1}, out{2, 2, 2}, out{2, 2, 3}] = ellipsoid (uplot, vplot);

  [uin, vin]     = ndgrid (linspace (.1, .9, nin), linspace (-.7, .7, min));
  [uplot, vplot] = ndgrid (linspace (.1, .9, nplot), linspace (-.7, .7, mplot));
  [out{3, 1, 1}, out{3, 1, 2}, out{3, 1, 3}] = zplane (uin, vin);
  [out{3, 2, 1}, out{3, 2, 2}, out{3, 2, 3}] = zplane (uplot, vplot);

  [uin, vin]     = ndgrid (linspace (.1, .9, nin), linspace (-.7, .7, min));
  [uplot, vplot] = ndgrid (linspace (.1, .9, nplot), linspace (-.7, .7, mplot));
  [out{4, 1, 1}, out{4, 1, 2}, out{4, 1, 3}] = yplane (uin, vin);
  [out{4, 2, 1}, out{4, 2, 2}, out{4, 2, 3}] = yplane (uplot, vplot);

  [uin, vin]     = ndgrid (linspace (.1, 1, nin), linspace (.1, 1, min));
  [uplot, vplot] = ndgrid (linspace (.1, 1, nplot), linspace (0, 1, mplot));
  [out{5, 1, 1}, out{5, 1, 2}, out{5, 1, 3}] = par (uin, vin);
  [out{5, 2, 1}, out{5, 2, 2}, out{5, 2, 3}] = par (uplot, vplot);

end


function [x, y, z] = cone (u, v)
  R = 10;
  H = 10;
  x = (R*v) .* cos (2*pi*u);
  y = (R*v) .* sin (2*pi*u);
  z = v*H;
end

function [x, y, z] = ellipsoid (u, v)
  A = 10;
  B = 4;
  C = 6;
  x = A * cos (2*pi*u) .* C .* cos (pi*v/2);
  y = B * sin (2*pi*u) .* C .* cos (pi*v/2);
  z = C * sin (pi*v/2);
end

function [x, y, z] = zplane (u, v)
  A = 10;
  B = 7;
  C = 6;
  x = A * u;
  y = B * v;
  z = C * ones (size (u));
end

function [x, y, z] = yplane (u, v)
  A = 10;
  B = 7;
  C = 6;
  x = A * u;
  z = B * v;
  y = C * ones (size (u));
end

function [x, y, z] = par (u, v)
  R = 10;
  H = 10;
  z = (H*u);
  y = (R*v);
  x = H*u.^2;
end
