function error = mahalanobis_error_em(realparamsfile,estparamsfile)
  load(realparamsfile);
  load(estparamsfile);
  len = length(realn1.mean_wm);
  real{1} = realn1;
  real{2} = realn5;
  real{3} = realn9;

  est{1} = emn1;
  est{2} = emn5;
  est{3} = emn9;
  error.em = mahalanobis_single_error(real,est);

  est{1} = initn1;
  est{2} = initn5;
  est{3} = initn9;
  error.eminit = mahalanobis_single_error(real,est);

  
function err = mahalanobis_single_error(real,est)
   
   for j = 1:3
     err(j) = sqrt((est{j}.mean_wm - real{j}.mean_wm)*inv(real{j}.var_wm)...
                   *(est{j}.mean_wm - real{j}.mean_wm)') + ...
              sqrt((est{j}.mean_gm - real{j}.mean_gm)*inv(real{j}.var_gm)...
                   *(est{j}.mean_gm - real{j}.mean_gm)') + ...
              sqrt((est{j}.mean_csf - real{j}.mean_csf)*inv(real{j}.var_csf)...
                   *(est{j}.mean_csf - real{j}.mean_csf)');
   end
