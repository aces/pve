% --------------------------------------------------------------------
% Parameter estimation for the partial volume classification
% --------------------------------------------------------------------
% March 5 2002 Jussi Tohka:

% A function for nuisance parameter estimation. 
% Given an image, a segmented image , and a brainmask it estimates
% nuisance parameters based on different estimators:
% 1. A simple maximum likelihood estimation where no outlier rejection is performed.
% 2. A maximum likelihood estimation where outliers are discovered (and rejected) based
%    on spatial information. i.e. voxels laying on boundaries of two tissue types are not
%    included in estimation.
% 3. Estimation based on Minimum Variance Ellipsoid estimator. Uses function minvarellipse.
% 4. Estimation based on Minimum Covariance Determinant estimator. Uses function fastmcdm.m
% 5. Minimum variance ellipsoid estimation based on trimmed segmentation like in 2.
% 6. Minimum covariance determinant estimation based trimmed segmentation like in 2.


%  Note that 3 and 4 are performed by approximatetive algorithms. This is for speed.
 
% Writes results also to a file(s). Assumes that the background mean intensity is 
% zero.
% Processes the whole volume at the time, so the memory usage can be a problem...

% Matlab function cov() produces minimum variance unbiased estimates and
% this function produces maximum likelihood estimates, which, however, do not 
% differ much if the number of data points is large.

% Input : segfile  : filename of segmented image, check that labels are correct!
%       : maskfile : filename for brainmask
%       : outfile  : filename for output
%       : varargin(1) : filename for image data
%       : varargin(2 - 3) [optional] : filenames for two other files of image data for 
%                                      multispectral case

% Output:
%         structs containing required parameter estimates


function [ml,trimmed_ml,mve, mcd, trimmed_mve, trimmed_mcd] = ...
param_estim_multi(segfile,maskfile,outfile,varargin);


  WM_NUMBER = 3;
  GM_NUMBER = 2;
  CSF_NUMBER = 1;
  TR = 0.5;
  DATA_TYPE  = 4096;
  SAMPLE_SIZE = 10000;

  tissue_type(1) = WM_NUMBER;
  tissue_type(2) = GM_NUMBER;
  tissue_type(3) = CSF_NUMBER;  

  structuring_element = zeros(3,3,3);
  structuring_element(2,2,1:3) = 1;
  structuring_element(1:3,2,2) = 1;
  structuring_element(2,1:3,2) = 1;

  nof_images = length(varargin);
  if nof_images == 1 | nof_images == 3 
    hin1 = openimage(varargin{1});
    hmask = openimage(maskfile);
    hseg  = openimage(segfile);
    if nof_images == 3
      hin2 = openimage(varargin{2}); 
      hin3 = openimage(varargin{3});
    end

    nslices = getimageinfo(hin1,'NumSlices');
    sz = getimageinfo(hin1,'ImageSize');
 
    for i = 1:6
      pcell{i}.means{1} = zeros(nof_images,1);
      pcell{i}.means{2} = zeros(nof_images,1);
      pcell{i}.means{3} = zeros(nof_images,1);
      pcell{i}.vars{1} = zeros(nof_images);
      pcell{i}.vars{2} = zeros(nof_images);
      pcell{i}.vars{3} = zeros(nof_images);
    end

    imgmask = reshape(getimages(hmask,1:nslices),prod(sz)*nslices,1);
    imgseg  = round(reshape(getimages(hseg,1:nslices),prod(sz)*nslices,1));
    if nof_images == 1
      imgin = reshape(getimages(hin1,1:nslices),prod(sz)*nslices,1);
    else
      imgin = zeros(prod(sz)*nslices,3);
      imgin(:,1) = reshape(getimages(hin1,1:nslices),prod(sz)*nslices,1);
      imgin(:,2) = reshape(getimages(hin2,1:nslices),prod(sz)*nslices,1);
      imgin(:,3) = reshape(getimages(hin3,1:nslices),prod(sz)*nslices,1);
    end 

    for tt = 1:3
      tt
      img_tmp = (imgseg == tissue_type(tt)) & (imgmask > TR);
      nvoxels = sum(img_tmp);
      pcell{1}.means{tt} = sum(repmat(img_tmp,1,nof_images).*imgin)/nvoxels;

      for i = 1:nof_images
        for j = 1:i     
          pcell{1}.vars{tt}(i,j) = sum((imgin(:,i) - pcell{1}.means{tt}(i)).*(imgin(:,j) - pcell{1}.means{tt}(j))...
                      .*img_tmp)/nvoxels;
          pcell{1}.vars{tt}(j,i) = pcell{1}.vars{tt}(i,j);
        end
      end
     
      
      data = zeros(nvoxels,nof_images);  
      for i = 1:nof_images
        data(:,i) = nonzeros(img_tmp.*imgin(:,i));
      end

      if nvoxels > 50000
        data = data(unidrnd(nvoxels,50000,1),:);
      end
      
      for i = 1:nof_images
        mi = min(data(:,i));
        ma = max(data(:,i));
        data(:,i) = data(:,i) + ((ma - mi)/(2*DATA_TYPE))*(rand(length(data(:,i)),1) - 0.5); 
      end
      tt
      [pcell{3}.means{tt} pcell{3}.vars{tt}] = minvarellipse(data);
      
      opt.alpha = 0.5;
      opt.lts = 1;
      tmpres = fastmcdm(data,opt);
      pcell{4}.means{tt} = tmpres.center;
      pcell{4}.vars{tt} = tmpres.cov;
      
      img_tmp = reshape(img_tmp,sz(2),sz(1),nslices);
      img_tmp = imerode(img_tmp,structuring_element); 
      img_tmp = double(reshape(img_tmp,sz(1)*sz(2)*nslices,1));

      nvoxels = sum(img_tmp);
      pcell{2}.means{tt} = sum(repmat(img_tmp,1,nof_images).*imgin)/nvoxels;
      
      for i = 1:nof_images
        for j = 1:i     
          pcell{2}.vars{tt}(i,j) = sum((imgin(:,i) - pcell{2}.means{tt}(i)).*(imgin(:,j) - pcell{2}.means{tt}(j))...
                      .*img_tmp)/nvoxels;
          pcell{2}.vars{tt}(j,i) = pcell{2}.vars{tt}(i,j);
        end
      end 
    

      data = zeros(nvoxels,nof_images);  
      for i = 1:nof_images
        data(:,i) = nonzeros(img_tmp.*imgin(:,i));
      end

      if nvoxels > 50000
        data = data(unidrnd(nvoxels,50000,1),:);
      end
      tt
      for i = 1:nof_images
        mi = min(data(:,i));
        ma = max(data(:,i));
        data(:,i) = data(:,i) + ((ma - mi)/(2*DATA_TYPE))*(rand(length(data(:,i)),1) - 0.5); 
      end
      tt
      [pcell{5}.means{tt} pcell{5}.vars{tt}] = minvarellipse(data);
      tt
      opt.alpha = 0.5;
      opt.lts = 1;
      tmpres = fastmcdm(data,opt);
      pcell{6}.means{tt} = tmpres.center;
      pcell{6}.vars{tt} = tmpres.cov;

    end
  
    var_bg = 10*eye(nof_images);
    var_mea = zeros(nof_images);
    mean_bg = zeros(nof_images,1);
    beta = 0.1;

    ml.mean_wm = pcell{1}.means{1};
    ml.mean_gm = pcell{1}.means{2};
    ml.mean_csf = pcell{1}.means{3};   
    ml.var_wm = pcell{1}.vars{1};
    ml.var_gm = pcell{1}.vars{2};
    ml.var_csf = pcell{1}.vars{3};    
    
    trimmed_ml.mean_wm = pcell{2}.means{1};
    trimmed_ml.mean_gm = pcell{2}.means{2};
    trimmed_ml.mean_csf = pcell{2}.means{3};   
    trimmed_ml.var_wm = pcell{2}.vars{1};
    trimmed_ml.var_gm = pcell{2}.vars{2};
    trimmed_ml.var_csf = pcell{2}.vars{3}; 

    mve.mean_wm = pcell{3}.means{1};
    mve.mean_gm = pcell{3}.means{2};
    mve.mean_csf = pcell{3}.means{3};   
    mve.var_wm = pcell{3}.vars{1};
    mve.var_gm = pcell{3}.vars{2};
    mve.var_csf = pcell{3}.vars{3};       

    mcd.mean_wm = pcell{4}.means{1};
    mcd.mean_gm = pcell{4}.means{2};
    mcd.mean_csf = pcell{4}.means{3};   
    mcd.var_wm = pcell{4}.vars{1};
    mcd.var_gm = pcell{4}.vars{2};
    mcd.var_csf = pcell{4}.vars{3};    

    trimmed_mve.mean_wm = pcell{5}.means{1};
    trimmed_mve.mean_gm = pcell{5}.means{2};
    trimmed_mve.mean_csf = pcell{5}.means{3};   
    trimmed_mve.var_wm = pcell{5}.vars{1};
    trimmed_mve.var_gm = pcell{5}.vars{2};
    trimmed_mve.var_csf = pcell{5}.vars{3};       

    trimmed_mcd.mean_wm = pcell{6}.means{1};
    trimmed_mcd.mean_gm = pcell{6}.means{2};
    trimmed_mcd.mean_csf = pcell{6}.means{3};   
    trimmed_mcd.var_wm = pcell{6}.vars{1};
    trimmed_mcd.var_gm = pcell{6}.vars{2};
    trimmed_mcd.var_csf = pcell{6}.vars{3};    



    printparamsfile(strcat(outfile,'.ml'),ml.mean_wm,ml.mean_gm,ml.mean_csf,mean_bg,...
                  ml.var_wm,ml.var_gm,ml.var_csf,var_bg,var_mea);
    printparamsfile(strcat(outfile,'.trimmed_ml'),trimmed_ml.mean_wm,trimmed_ml.mean_gm,...
                                                  trimmed_ml.mean_csf,mean_bg,...
                  trimmed_ml.var_wm,trimmed_ml.var_gm,trimmed_ml.var_csf,var_bg,var_mea);
    printparamsfile(strcat(outfile,'.mcd'),mcd.mean_wm,mcd.mean_gm,mcd.mean_csf,mean_bg,...
                  mcd.var_wm,mcd.var_gm,mcd.var_csf,var_bg,var_mea);
    printparamsfile(strcat(outfile,'.mve'),mve.mean_wm,mve.mean_gm,mve.mean_csf,mean_bg,...
                  mve.var_wm,mve.var_gm,mve.var_csf,var_bg,var_mea);

     printparamsfile(strcat(outfile,'.trimmed_mve'),trimmed_mve.mean_wm,trimmed_mve.mean_gm,...
                                                  trimmed_mve.mean_csf,mean_bg,...
                  trimmed_mve.var_wm,trimmed_mve.var_gm,trimmed_mve.var_csf,var_bg,var_mea);

     printparamsfile(strcat(outfile,'.trimmed_mcd'),trimmed_mcd.mean_wm,trimmed_mcd.mean_gm,...
                                                 trimmed_mcd.mean_csf,mean_bg,...
                  trimmed_mcd.var_wm,trimmed_mcd.var_gm,trimmed_mcd.var_csf,var_bg,var_mea);

  else  % There was a wrong number of input parameters  
    fprintf('There was a wrong number of input parameters');
  end
 
