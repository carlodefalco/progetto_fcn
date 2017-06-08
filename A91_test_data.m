
%% per vedere i dati dell'esempio usare il comando
%% out = A91_test_data (nin, nplot)
%% for ii = 1 : size (out, 1)
%%   figure
%%   plot3 (out{ii, 1}(1, :), out{ii, 1}(2, :), out{ii, 1}(3, :), 'o',
%%          out{ii, 2}(1, :), out{ii, 2}(2, :), out{ii, 2}(3, :));
%% end


function out = A91_test_data (nin, nplot)

  xin = linspace (-5, 5, nin);
  xplot = linspace (-5, 5, nplot);
  out{1, 1} = [xin; runge(xin); zeros(size (xin))];
  out{1, 2} = [xplot; runge(xplot); zeros(size (xplot))];

  xin = linspace (0, 1, nin);
  xplot = linspace (0, 1, nplot);
  [x, y, z] = circle (xin);
  out{2, 1} = [x; y; z];
  [x, y, z] = circle (xplot);
  out{2, 2} = [x; y; z];

  xin = linspace (0, 5, nin);
  xplot = linspace (0, 5, nplot);
  [x, y, z] = helix (xin);
  out{3, 1} = [x; y; z];
  [x, y, z] = helix (xplot);
  out{3, 2} = [x; y; z];

  xin = linspace (0, 1, nin);
  xplot = linspace (0, 1, nplot);
  [x, y, z] = crown (xin);
  out{4, 1} = [x; y; z];
  [x, y, z] = crown (xplot);
  out{4, 2} = [x; y; z];

end


function y = runge (x)
  y = 1 ./ (1 + x.^2);
end

function [x, y, z] = circle (theta)
  R = 10;
  x = R * cos (2 * pi * theta);
  y = R * sin (2 * pi * theta);
  z = zeros (size (theta));
end

function [x, y, z] = helix (theta)
  R = 10;
  P = 2;
  x = R * cos (2 * pi * theta);
  y = R * sin (2 * pi * theta);
  z = P * theta;
end

function [x, y, z] = crown (theta)
  R = 10;
  P = 2;
  omega = 9;
  x = R * cos (2 * pi * theta);
  y = R * sin (2 * pi * theta);
  z = sin (2*pi*omega*theta);
end
