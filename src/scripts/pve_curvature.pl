#!xPERLx

use strict;
use warnings "all";
use Getopt::Tabular;
use File::Basename;
use File::Temp qw/ tempdir /;

my($Help, $Usage, $me);
my(@opt_table, %opt );

$me = &basename($0);
%opt = (
   'verbose'   => 1,
   'clobber'   => 0,
   );

@opt_table = (
   ["-verbose", "boolean", 0, \$opt{verbose},
      "be verbose" ],
   ["-clobber", "boolean", 0, \$opt{clobber},
      "clobber existing files" ]
);

GetOptions( \@opt_table, \@ARGV )
  or exit 1;

&CreateInfoText;

# define input variables:

my $spect = $ARGV[0];
my $class = $ARGV[1];
my $mask = $ARGV[2];
my $output = $ARGV[3];
$output .= "_cg.mnc";

my @clobber_opt;
if ($opt{clobber}) {
    push @clobber_opt,"-clobber";
}
if (-e $output && !$opt{clobber} ) {
    print "Output file $output exists and -clobber not specified\n";
    die "Output file $output exists and -clobber not specified\n";
}

my $temp_dir = &tempdir( "$me-XXXXXXXX", TMPDIR => 1, CLEANUP => 1 );

my $maskdil = "$temp_dir/dil.mnc";
my $classmasked = "$temp_dir/classmasked.mnc";
my $classdil = "$temp_dir/classdil.mnc";
my $white = "$temp_dir/white.mnc";
my $greycsf = "$temp_dir/greycsf.mnc";
my $whitegrp = "$temp_dir/whitegrp.mnc";
my $greycsfgrp = "$temp_dir/greycsfgrp.mnc";
my $greycsf1 = "$temp_dir/greycsf1.mnc";
my $greycsfe = "$temp_dir/greycsfe.mnc";
my $greycsf50 = "$temp_dir/greycsf50.mnc";
my $chamf = "$temp_dir/chamf.mnc";
my $dx = "$temp_dir/dx.mnc";
my $dy = "$temp_dir/dy.mnc";
my $dz = "$temp_dir/dz.mnc";
my $dxyz ="$temp_dir/dxyz.mnc";
my $smooth = "$temp_dir/smooth.mnc";
my $smoothmask = "$temp_dir/smoothmask.mnc";
my $smooths1 = "$temp_dir/smooths1.mnc";
my $smoothm1 = "$temp_dir/smoothm1.mnc";
my $smoothlt0 = "$temp_dir/smoothlt0.mnc";
my $smoothlt0m3 = "$temp_dir/grad.mnc";
my $inv = "$temp_dir/inv.mnc";
my $blurpref = "$temp_dir/temp";
my $blur = "$temp_dir/temp_blur.mnc";
my $curvepref = "$temp_dir/temp";
my $curve = "$temp_dir/temp_k1.mnc";
my $curvelt0 = "$temp_dir/curve.mnc";
my $cplusg = "$temp_dir/cplusg.mnc";

my $kerndx = "$temp_dir/dx.kern";
my $kerndy = "$temp_dir/dy.kern";
my $kerndz = "$temp_dir/dz.kern";
&CreateKernelFiles;

######create the gradient volume#####

#This dilation of the CSF is used to ensure that the
#gradient image created is well outside the GM/CSF interface

&run( "dilate_volume",$mask,$maskdil,"1","26","4" );
&run( "minccalc","-expression", 'if (A[1] > 0) A[0] else 0', $class,$mask,$classmasked );

#creation of volume with dilated CSF
&run( "minccalc","-expression", 'if (A[0] > 0) A[0] else if (A[1] > 0) 1  else 0',
      $classmasked,$maskdil,$classdil );

#These steps are used to create a GM/CSF without interference
#from isolated white matter voxels in the GM, which are removed
#by grouping clusters in both the WM and GMCSF and replacing
#WM and GMCSF outside the largest cluster with the opposite grouping.

&run( "minccalc","-expression", 'A[0] == 3',$classdil,$white );
&run( "minccalc","-expression", 'A[0] > 0.5 && A[0] < 2.5',$classdil,$greycsf );
&run( "mincmorph","-group",$white,$whitegrp );
&run( "mincmorph","-group",$greycsf,$greycsfgrp );
&run( "minccalc","-expression", 'A[0] > 1.5 || A[1] > 0.9 && A[1] < 1.1',
       $whitegrp,$greycsfgrp,$greycsf1 );
&run( "mincmorph","-erosion","$greycsf1","$greycsfe" );

#the GMCSF is inverted and then a distance transform is calculated
&run( "minccalc","-expression", '!(A[0])',$greycsf1,$inv );
&run( "mincreshape","-image_range","0","50",$inv,$greycsf50 );
&run( "mincchamfer",$greycsf50,$chamf );

#derivatives of the distance transform are taken and a
#derivative magnitude volume is created.  (Note that this
#procedure is NOT equivalent to use make_gradient_volume from
#conglomerate
&run( "mincmorph","-convolve","-kern",$kerndx,$chamf,$dx );
&run( "mincmorph","-convolve","-kern",$kerndy,$chamf,$dy );
&run( "mincmorph","-convolve","-kern",$kerndz,$chamf,$dz );
&run( "minccalc","-expression", "sqrt(A[0]^2+A[1]^2+A[2]^2)",$dx,$dy,$dz,$dxyz );

#Geometric smoothing is performed to remove noise and increase the 
#width of the areas where we expect find CSF
&run( "geo_smooth","0.1","4",$dxyz,$smooth );

#remove some noisy voxels from the grad image, and convert to a value
#useful to the pve.
&run( "minccalc","-expression", 'if (A[1]) { if ((A[0]-1) > 0) 0 else 3*(A[0]-1); } else 0',
       $smooth,$greycsf1,$smoothlt0m3 );

# create the curvature volume after blurring the input image
&run( "mincblur","-fwhm","4",$spect,$blurpref );
&run( "make_curvature_volume",$blur,$curvepref );
&run( "minccalc","-expression", 'if (A[0] < 0) A[0] else 0',$curve,$curvelt0 );

# create the combination grad+curve volume used for pve and mask out non GM/CSF
&run( "minccalc","-expression", 'A[0]+0.5*A[1]',$curvelt0,$smoothlt0m3,$cplusg );
&run( "minccalc", @clobber_opt, "-expression", 'if (A[1]) A[0] else 0',
      $cplusg, $greycsfe, $output );

sub CreateInfoText
{
  $Usage = <<USAGE;
usage:
       $me [options] t1_image.mnc classified.mnc mask.mnc output_pref
       $me -help to list options
USAGE
  
  $Help = <<HELP;
  
$me takes an intensity image, a classified image and a mask
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

#Execute a system call.

sub run {
  print "@_\n";
  system(@_)==0 or die "Command @_ failed with status: $?";
}


