
%% per vedere i dati dell'esempio usare il comando
%% out = A92_test_data (10, 100)
%% for ii = 1 : size (out, 1)
%%   figure
%%   plot3 (out{ii, 1}(1, :), out{ii, 1}(2, :), out{ii, 1}(3, :), 'o',
%%          out{ii, 2}(1, :), out{ii, 2}(2, :), out{ii, 2}(3, :));
%%   hold on
%%   quiver3 (out{ii, 1}(1, [1 end]), out{ii, 1}(2, [1 end]), out{ii,1}(3, [1 end]),
%%            out{ii, 3}(1, :), out{ii, 3}(2, :), out{ii, 3}(3, :));
%%   hold off
%% end


function out = A92_test_data (nin, nplot)

  xin = linspace (-5, 5, nin);
  xplot = linspace (-5, 5, nplot);
  out{1, 1} = [xin; runge(xin); zeros(size (xin))];
  out{1, 2} = [xplot; runge(xplot); zeros(size (xplot))];
  out{1, 3} = [ones(1, 2); rungeder(xin([1, end])); zeros(1, 2)];

  xin = linspace (0, .9, nin);
  xplot = linspace (0, .9, nplot);
  [x, y, z] = circle (xin);              out{2, 1} = [x; y; z];
  [x, y, z] = circle (xplot);            out{2, 2} = [x; y; z];
  [x, y, z] = circleder (xin([1, end])); out{2, 3} = [x; y; z];
  
  xin = linspace (0, 5, nin);
  xplot = linspace (0, 5, nplot);
  [x, y, z] = helix (xin);               out{3, 1} = [x; y; z];
  [x, y, z] = helix (xplot);             out{3, 2} = [x; y; z];
  [x, y, z] = helixder (xin([1, end]));  out{3, 3} = [x; y; z];

  xin = linspace (0, .9, nin);
  xplot = linspace (0, .9, nplot);
  [x, y, z] = crown (xin);               out{4, 1} = [x; y; z];
  [x, y, z] = crown (xplot);             out{4, 2} = [x; y; z];
  [x, y, z] = crownder (xin([1, end]));  out{4, 3} = [x; y; z];

end


function y = runge (x)
  y = 1 ./ (1 + x.^2);
end

function y = rungeder (x)
  y = - 2*x ./ (1 + x.^2).^2;
end

function [x, y, z] = circle (theta)
  R = 10;
  x = R * cos (2 * pi * theta);
  y = R * sin (2 * pi * theta);
  z = zeros (size (theta));
end

function [x, y, z] = circleder (theta)
  R = 10;
  x = - 2 * pi * R * sin (2 * pi * theta);
  y =   2 * pi * R * cos (2 * pi * theta);
  z = zeros (size (theta));
end

function [x, y, z] = helix (theta)
  R = 10;
  P = 2;
  x = R * cos (2 * pi * theta);
  y = R * sin (2 * pi * theta);
  z = P * theta;
end

function [x, y, z] = helixder (theta)
  R = 10;
  P = 2;
  x = - 2 * pi * R * sin (2 * pi * theta);
  y =   2 * pi * R * cos (2 * pi * theta);
  z = P * ones (size (theta));
end

function [x, y, z] = crown (theta)
  R = 10;
  P = 2;
  omega = 9;
  x = R * cos (2 * pi * theta);
  y = R * sin (2 * pi * theta);
  z = sin (2*pi*omega*theta);
end

function [x, y, z] = crownder (theta)
  R = 10;
  P = 2;
  omega = 9;
  x = - 2 * pi * R * sin (2 * pi * theta);
  y =   2 * pi * R * cos (2 * pi * theta);
  z = 2 * pi * omega * cos (2*pi*omega*theta);
end
