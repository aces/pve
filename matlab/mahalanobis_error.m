function error = mahalanobis_error(realparamsfile,estparamsfile)
  load(realparamsfile);
  load(estparamsfile);
  len = length(realn1.mean_wm);
  real{1} = realn1;
  real{2} = realn5;
  real{3} = realn9;

  est{1} = mln1;
  est{2} = mln5;
  est{3} = mln9;
  error.ml = mahalanobis_single_error(real,est);

  est{1} = mven1;
  est{2} = mven5;
  est{3} = mven9;
  error.mve = mahalanobis_single_error(real,est);

  est{1} = mcdn1;
  est{2} = mcdn5;
  est{3} = mcdn9;
  error.mcd = mahalanobis_single_error(real,est);

  est{1} = tmln1;
  est{2} = tmln5;
  est{3} = tmln9;
  error.tml = mahalanobis_single_error(real,est);

  est{1} = tmven1;
  est{2} = tmven5;
  est{3} = tmven9;
  error.tmve = mahalanobis_single_error(real,est);

  est{1} = tmcdn1;
  est{2} = tmcdn5;
  est{3} = tmcdn9;
  error.tmcd = mahalanobis_single_error(real,est);
function err = mahalanobis_single_error(real,est)
   
   for j = 1:3
     err(j) = sqrt((est{j}.mean_wm - real{j}.mean_wm)*inv(real{j}.var_wm)...
                   *(est{j}.mean_wm - real{j}.mean_wm)') + ...
              sqrt((est{j}.mean_gm - real{j}.mean_gm)*inv(real{j}.var_gm)...
                   *(est{j}.mean_gm - real{j}.mean_gm)') + ...
              sqrt((est{j}.mean_csf - real{j}.mean_csf)*inv(real{j}.var_csf)...
                   *(est{j}.mean_csf - real{j}.mean_csf)');
   end
