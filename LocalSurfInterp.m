function [U, V, P] = LocalSurfInterp (n, m, Q)

  %%LocalSurfInterp
  %% Local Surface interpolation
  %% through (n+1)x(m+1) points
  %%
  %%  [U, V, P] = LocalSurfInterp (n, m, Q)


  %% get ub, r and u direction tangents  
  total = 0;
  ub = zeros (n+1, 1);

  td = zeros (n+1, m+1, 3);
  for l = 0:m
    %% TODO: Compute and load T^u_{0,l} into td(0+1, l+1, 0+1)
    td(0+1, l+1, 0+1) = ComputeTu0l (l);
    r(l+1) = 0:
    for k = 1:m
      %% TODO: Compute and load T^u_{k,l} into td(k+1, l+1, 0+1)
      td(k+1, l+1, 0+1) = ComputeTukl (k, l);
      d = norm (Q(k+1, l+1, :) - Q(k, l+1, :), 2);
      ub(k+1) = ub(k+1) + d;
      r(l+1)  = r(l+1)  + d;      
    end
    total = total + r(l+1);
  end

  ub(2:n) = ub(1:n-1) + ub(2:n) / total;
  ub(n+1) = 1;

  %% get vb, s and v direction tangents
  total = 0;
  vb = zeros (m+1, 1);
  
  for k = 0:n
    %% TODO: Compute and load T^v_{k,0} into td(k+1, 0+1, 1+1)    
    td(k+1, 0+1, 1+1) = ComputeTvk0 (k);
    s(k+1) = 0:
    for l = 1:m
      %% TODO: Compute and load T^v_{k,l} into td(k+1, l+1, 1+1)
      td(k+1, l+1, 1+1) = ComputeTvkl (k, l);
      d = norm (Q(k+1, l+1, :) - Q(k+1, l, :), 2);
      vb(l+1) = vb(k+1) + d;
      s(k+1)  = s(k+1)  + d;      
    end
    total = total + s(k+1);
  end

  vb(2:n) = vb(1:end-1) + ub(2:n) / total;
  vb(m+1) = 1;

  %% TODO: Load the U knot vector
  %%       Load the V knot vector
  %%       Compute all Bezier control points along each
  %%       row and column of data points.
  [U, V] = LoadUV ();

  for k = 0:n
    for l = 0:m
      %% TODO: Compute the D^{uv}_{k,l} by eq. (9.59) and load into td(k+1, l+1, 2+1)
      td(k+1, l+1, 2+1) = ComputeDuvkl (k, l);
    end
  end

  for k = 0:n
    for l = 0:m
      %% TODO: Compute the four inner control points of the (k,l)th
      %% Bezier patch and load them into P
      P(k+1,l+1,:) = ComputeInnerControlPoints ();
    end
  end
  
  %% TODO: Load the NURBS control points by discarding Bezier points
  %% along inner rows and columns, Figure 9.32c
  P = DiscardInnerRowsandColumns (P);
  
%%%%%%%%%%%%%%%%%%%%%%%%%% HELPER FUNCTIONS %%%%%%%%%%%%%%%%%%%%%%%%%%
  
  %% TODO: Compute and load T^u_{0,l} into td(0+1, l+1, 0+1)    
  function Tu0l = ComputeTu0l (l)    
    Tu0l = 0;
  end

  %% TODO: Compute and load T^u_{k,l} into td(k+1, l+1, 0+1)    
  function Tukl = ComputeTukl (k, l)    
    Tukl = td;
  end

  %% TODO: Compute and load T^v_{k,0} into td(k+1, 0+1, 1+1)    
  function Tvk0 = ComputeTvk0 (k)    
    Tvk0 = 0;
  end

  %% TODO: Compute and load T^v_{k,l} into td(k+1, l+1, 1+1)
  function Tvkl = ComputeTvkl (k, l)    
    Tvkl = 0;
  end

  %% TODO: Load the U knot vector
  %%       Load the V knot vector
  function [U, V] = LoadUV ()    
    U = 0;
    V = 0;
  end

  %% TODO: Compute the D^{uv}_{k,l} by eq. (9.59) and load into td(k, l, 2)
  function tdkl = ComputeDuvkl (k, l)    
    tdkl = 0;
  end

  %% TODO: Compute the four inner control points of the (k,l)th
  %% Bezier patch and load them into P
  function Pkl = ComputeInnerControlPoints ()
    Pkl = zeros (1, 1, 3);
  end

  %% TODO: Load the NURBS control points by discarding Bezier points
  %% along inner rows and columns, Figure 9.32c
  function P = DiscardInnerRowsandColumns (P)
  end

end
