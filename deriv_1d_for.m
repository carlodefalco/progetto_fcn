function Dk = deriv_1d_for (uk, Qk)

  n = numel (uk) - 1;
  
  for k = 0 : n
    if (k == 0)

      Duk = diff (uk(1:3));
      ak  = Duk(1) / (Duk(1) + Duk(2)); 
      qk  = diff (Qk(:, 1:3), 1, 2);
      dk  = qk ./ Duk;
      dk(:,2)  = (1 - ak) * dk(:,1) + ak * dk(:,2);
      c1  =  2;
      c2  = -1;
      
    elseif (k == n)

      Duk = diff (uk(n-1:n+1));
      ak  = Duk(1) / (Duk(1) + Duk(2)); 
      qk  = diff (Qk(:, n-1:n+1), 1, 2);
      dk  = qk ./ Duk;
      dk(:,1) = (1 - ak) * dk(:,1) + ak * dk(:,2);
      c1  = -1;
      c2  =  2;

    else
      
      Duk = diff (uk(k:k+2));
      ak  = Duk(1) / (Duk(1) + Duk(2)); 
      qk  = diff (Qk(:, k:k+2), 1, 2);
      dk  = qk ./ Duk;
      c1  = (1 - ak);
      c2  = ak;

    end
    
    Dk(:, k+1) = c1 * dk(:,1) + c2 * dk(:,2);

  end
end
