% pve_error.m : Jussi Tohka 12 March 2002
% Computes the error measure for partial volume estimation based on L1-norm
% error = 1/N*(Sum_{i = 1}^N ||e_i||_1),
% where ||e_i||_1 is the L1-norm of the residual between ith "real" partial
% volume vector and the estimated one. N is the number of voxels in the
%  brain area. Calculates also the dice coefficient and misclasiification 
% rate (for convinience).


function [error,dice,mcr] = pve_error(filename,mask_filename,varargin);

   MASK_TR = 0.5;
   CSFLABEL = 1;
   GMLABEL = 2;
   WMLABEL = 3; 
     
   csf_filename = strcat(filename,'_csf.mnc');
   gm_filename  = strcat(filename,'_gm.mnc');
   wm_filename  = strcat(filename,'_wm.mnc');

   h_csf = openimage(csf_filename);
   h_gm  = openimage(gm_filename);
   h_wm  = openimage(wm_filename);
 
   if length(varargin) == 0
     h_ref_csf = openimage('/data/sbd/phantoms/phantom_1.0mm_normal_csf.mnc.gz');
     h_ref_gm  = openimage('/data/sbd/phantoms/phantom_1.0mm_normal_gry.mnc.gz');
     h_ref_gl  = openimage('/data/sbd/phantoms/phantom_1.0mm_normal_gli.mnc.gz');
     h_ref_wm  = openimage('/data/sbd/phantoms/phantom_1.0mm_normal_wht.mnc.gz');
     h_ref_crisp = openimage('/data/sendai/sendai1/temp/jupeto/brainwebimages/phantom_1.0mm_brainonly.mnc');
   else
     h_ref_csf = openimage(varargin{1});
     h_ref_gm  = openimage(varargin{2});
     h_ref_wm  = openimage(varargin{3});
     h_ref_gl  = openimage(varargin{4});
     h_ref_crisp = openimage(varargin{5});
   end

   h_mask = openimage(mask_filename);

   nslices = getimageinfo(h_mask,'NumSlices');
   error = 0;

   mask_img = getimages(h_mask,1:nslices);
   img      = getimages(h_csf,1:nslices);
   ref_img  = getimages(h_ref_csf,1:nslices);
   
   brain_voxels = sum(sum(mask_img > MASK_TR));
   
   seg_img = (img > 0.5) & (mask_img > 0.5);   

   img = abs(img - ref_img); 
   img = img.*(mask_img > MASK_TR);
   error = error + ... 
        (1/brain_voxels)*sum(img(:));

   ref_img = getimages(h_ref_crisp,1:nslices);
   ref_img = (round(ref_img) == CSFLABEL) & (mask_img > MASK_TR); 

   dice(CSFLABEL) = 2*length(nonzeros(ref_img(:) & seg_img(:)))/...
   (length(nonzeros(ref_img(:))) + length(nonzeros(seg_img(:))));

   mcr = length(nonzeros(ref_img(:) & ~seg_img(:)))/brain_voxels;  

   %Gray matter

   img      = getimages(h_gm,1:nslices);
   ref_img  = getimages(h_ref_gm,1:nslices) + getimages(h_ref_gl,1:nslices);;
   
   seg_img = (img > 0.5) & (mask_img > 0.5);  

   img = abs(img - ref_img); 
   img = img.*(mask_img > MASK_TR);
   error = error + ... 
        (1/brain_voxels)*sum(img(:));

   ref_img = getimages(h_ref_crisp,1:nslices);
   ref_img = (round(ref_img) == GMLABEL) & (mask_img > MASK_TR); 

   dice(GMLABEL) = 2*length(nonzeros(ref_img(:) & seg_img(:)))/...
   (length(nonzeros(ref_img(:))) + length(nonzeros(seg_img(:))));

   mcr = mcr + length(nonzeros(ref_img(:) & ~seg_img(:)))/brain_voxels;    

   % White matter 
   img      = getimages(h_wm,1:nslices);
   ref_img  = getimages(h_ref_wm,1:nslices);
   
   seg_img = (img > 0.5) & (mask_img > 0.5);  
   img = abs(img - ref_img); 
   img = img.*(mask_img > MASK_TR);
   error = error + ... 
        (1/brain_voxels)*sum(img(:));  

   ref_img = getimages(h_ref_crisp,1:nslices);
   ref_img = (round(ref_img) == WMLABEL) & (mask_img > MASK_TR); 

   dice(WMLABEL) = 2*length(nonzeros(ref_img(:) & seg_img(:)))/...
   (length(nonzeros(ref_img(:))) + length(nonzeros(seg_img(:))));
   
   mcr = mcr + length(nonzeros(ref_img(:) & ~seg_img(:)))/brain_voxels;

   closeimage(h_csf);
   closeimage(h_gm);
   closeimage(h_wm);
   closeimage(h_mask);
   closeimage(h_ref_csf);
   closeimage(h_ref_wm);
   closeimage(h_ref_gm);


  










