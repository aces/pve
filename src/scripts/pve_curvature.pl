#!xPERLx -w

use strict;

use MNI::Startup;
use MNI::Spawn;
use Getopt::Tabular;
use File::Basename;

# global variables
my ($Help, $Usage);
my ($spect,$class,$mask,$output);


self_announce if $Verbose;

# Set defaults for the global variables.
$Verbose      = 1;
$Execute      = 1;
$Clobber      = 0;
$KeepTmp      = 0;

&CreateInfoText;

my (@argTbl) = (@DefaultArgs);

my(@leftOverArgs);
GetOptions(\@argTbl, \@ARGV, \@leftOverArgs) || die "\n";

if ($#leftOverArgs != 3) {
  warn $Usage;
  die "Incorrect number of arguments\n"; 
}

$spect = shift @leftOverArgs;
$class = shift @leftOverArgs;
$mask = shift @leftOverArgs;
$output = shift @leftOverArgs;
$output .= "_cg.mnc";

if (-e "${output}_cg.mnc" && !$Clobber) {
    die "Output file $output exists and -clobber not specified\n";
}


RegisterPrograms( ['minccalc','dilate_volume','mincmorph','mincreshape',
		  'mincchamfer','geo_smooth','mincmath','mincblur',
		  'make_curvature_volume'] );

AddDefaultArgs('mincmath',['-clobber']) if ($Clobber);
AddDefaultArgs('minccalc',['-expression']);


system("mkdir $TmpDir") unless ( -d $TmpDir );

my $maskdil = "$TmpDir/dil.mnc";
my $classmasked = "$TmpDir/classmasked.mnc";
my $classdil = "$TmpDir/classdil.mnc";
my $white = "$TmpDir/white.mnc";
my $greycsf = "$TmpDir/greycsf.mnc";
my $whitegrp = "$TmpDir/whitegrp.mnc";
my $greycsfgrp = "$TmpDir/greycsfgrp.mnc";
my $greycsf1 = "$TmpDir/greycsf1.mnc";
my $greycsfe = "$TmpDir/greycsfe.mnc";
my $greycsf50 = "$TmpDir/greycsf50.mnc";
my $chamf = "$TmpDir/chamf.mnc";
my $dx = "$TmpDir/dx.mnc";
my $dy = "$TmpDir/dy.mnc";
my $dz = "$TmpDir/dz.mnc";
my $dxyz ="$TmpDir/dxyz.mnc";
my $smooth = "$TmpDir/smooth.mnc";
my $smoothmask = "$TmpDir/smoothmask.mnc";
my $smooths1 = "$TmpDir/smooths1.mnc";
my $smoothm1 = "$TmpDir/smoothm1.mnc";
my $smoothlt0 = "$TmpDir/smoothlt0.mnc";
my $smoothlt0m3 = "$TmpDir/grad.mnc";
my $inv = "$TmpDir/inv.mnc";
my $blurpref = "$TmpDir/temp";
my $blur = "$TmpDir/temp_blur.mnc";
my $curvepref = "$TmpDir/temp";
my $curve = "$TmpDir/temp_k1.mnc";
my $curvelt0 = "$TmpDir/curve.mnc";
my $cplusg = "$TmpDir/cplusg.mnc";

my $kerndx = "$TmpDir/dx.kern";
my $kerndy = "$TmpDir/dy.kern";
my $kerndz = "$TmpDir/dz.kern";
&CreateKernelFiles;

######create the gradient volume#####

#This dilation of the CSF is used to ensure that the
#gradient image created is well outside the GM/CSF interface
Spawn(["dilate_volume",$mask,$maskdil,"1","26","4"]);
Spawn(["minccalc",'if (A[1] > 0) A[0] else 0',
       $class,$mask,$classmasked]);
#creation of volume with dilated CSF
Spawn(["minccalc",'if (A[0] > 0) A[0] else if (A[1] > 0) 1  else 0',
       $classmasked,$maskdil,$classdil]);

#These steps are used to create a GM/CSF without interference
#from isolated white matter voxels in the GM, which are removed
#by grouping clusters in both the WM and GMCSF and replacing
#WM and GMCSF outside the largest cluster with the opposite grouping.
Spawn(["minccalc",'A[0] == 3',$classdil,$white]);
Spawn(["minccalc",'A[0] > 0.5 && A[0] < 2.5',$classdil,$greycsf]);
Spawn(["mincmorph","-group",$white,$whitegrp]);
Spawn(["mincmorph","-group",$greycsf,$greycsfgrp]);
Spawn(["minccalc",'A[0] > 1.5 || A[1] > 0.9 && A[1] < 1.1',
       $whitegrp,$greycsfgrp,$greycsf1]);
Spawn(["mincmorph","-erosion","$greycsf1","$greycsfe"]);

#the GMCSF is inverted and then a distance transform is calculated
Spawn(["minccalc",'!(A[0])',$greycsf1,$inv]);
Spawn(["mincreshape","-image_range","0","50",$inv,$greycsf50]);
Spawn(["mincchamfer",$greycsf50,$chamf]);

#derivatives of the distance transform are taken and a
#derivative magnitude volume is created.  (Note that this
#procedure is NOT equivalent to use make_gradient_volume from
#conglomerate
Spawn(["mincmorph","-convolve","-kern",$kerndx,$chamf,$dx]);
Spawn(["mincmorph","-convolve","-kern",$kerndy,$chamf,$dy]);
Spawn(["mincmorph","-convolve","-kern",$kerndz,$chamf,$dz]);
Spawn(["minccalc","sqrt(A[0]^2+A[1]^2+A[2]^2)",$dx,$dy,$dz,$dxyz]);

#Geometric smoothing is performed to remove noise and increase the 
#width of the areas where we expect find CSF
Spawn(["geo_smooth","0.1","4",$dxyz,$smooth]);

#remove some noisy voxels from the grad image, and convert to a value
#useful to the pve.
Spawn(["minccalc",'if (A[1]) { if ((A[0]-1) > 0) 0 else 3*(A[0]-1); } else 0',
       $smooth,$greycsf1,$smoothlt0m3]);

# create the curvature volume after blurring the input image
Spawn(["mincblur","-fwhm","4",$spect,$blurpref]);
Spawn(["make_curvature_volume",$blur,$curvepref]);
Spawn(["minccalc",'if (A[0] < 0) A[0] else 0',$curve,$curvelt0]);

# create the combination grad+curve volume used for pve and mask out non GM/CSF
Spawn(["minccalc",'A[0]+0.5*A[1]',$curvelt0,$smoothlt0m3,$cplusg]);
Spawn(["minccalc",'if (A[1]) A[0] else 0',$cplusg,$greycsfe,$output]);

sub CreateInfoText
{
  $Usage = <<USAGE;
usage:
       $ProgramName [options] t1_image.mnc classified.mnc mask.mnc output_pref
       $ProgramName -help to list options
USAGE
  
  $Help = <<HELP;
  
$ProgramName takes an intensity image, a classified image and a mask
as input and generates an image representing the degree to which
the pve/pve3 classifiers should be biased to finding CSF through
the modulation of a Markovian Random Field.
HELP

   Getopt::Tabular::SetHelp ($Help, $Usage);
}

#Create the kernel files used in taking the image derivative
sub CreateKernelFiles
{
 
  open DXFILE, "> $kerndx";
  print DXFILE "MNI Morphology Kernel File\n";
  print DXFILE "Kernel_Type = Normal_Kernel;\n";
  print DXFILE "Kernel =\n";
  print DXFILE "  -1.0  0.0  0.0  0.0  0.0     -0.5\n";
  print DXFILE "   1.0  0.0  0.0  0.0  0.0      0.5;\n";
  close DXFILE;


  open DYFILE, "> $kerndy";
  print DYFILE "MNI Morphology Kernel File\n";
  print DYFILE "Kernel_Type = Normal_Kernel;\n";
  print DYFILE "Kernel =\n";
  print DYFILE "   0.0 -1.0  0.0  0.0  0.0     -0.5\n";
  print DYFILE "   0.0  1.0  0.0  0.0  0.0      0.5;\n";
  close DYFILE;

  open DZFILE, "> $kerndz";
  print DZFILE "MNI Morphology Kernel File\n";
  print DZFILE "Kernel_Type = Normal_Kernel;\n";
  print DZFILE "Kernel =\n";
  print DZFILE "   0.0  0.0 -1.0  0.0  0.0     -0.5\n";
  print DZFILE "   0.0  0.0  1.0  0.0  0.0      0.5;\n";
  close DZFILE;

}
