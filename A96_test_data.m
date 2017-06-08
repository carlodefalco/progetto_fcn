


function out = A96_test_data (nin, nd, nplot)

  tin       = linspace (-5, 5, nin);
  tplot     = linspace (-5, 5, nplot);

  [x, y, z] = crown (tin);
  out{1, 1} = [x; y; z];
  for ii = 1 : nd
    [x, y, z] = crownder (tin, ii);
    out{1, ii+1} = [x; y; z];
  end

  [x, y, z] = crown (tplot);
  out{1, nd+2} = [x; y; z];

end


function [x, y, z] = crown (theta)

  R     = 10;
  P     = 2;
  omega = 9;

  x = R * cos (2 * pi * theta);
  y = R * sin (2 * pi * theta);
  z = sin (2*pi*omega*theta);

end

function [x, y, z] = crownder (theta, n)

  R     = 10;
  P     = 2;
  omega = 9;

  x = R * (2 * pi)^n * cos (2 * pi * theta - n * pi/2);
  y = R * (2 * pi)^n * sin (2 * pi * theta + n * pi/2);
  z = (2 * pi * omega)^n * sin (2 * pi * omega * theta + n * pi/2); 

end
