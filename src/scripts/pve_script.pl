#!xPERLx -w

use strict;

use MNI::Startup;
use MNI::Spawn;
use Getopt::Tabular;
use File::Basename; 

my $help = <<HELP;

This is 'pve_script.pl'  It is essentially a wrapper around
the pve/pve3 programs that automatically determines which of
these 2 programs you would like to run based on the \# of
input files and also includes a subcortical class in the
classification if the '-subcortical' option is specified, among
other things.  If you would like to use a subcortical volume other 
than the default you can call one of the binaries directly with a
'-subcortical <file.mnc>' flag.

HELP

my $usage = <<USAGE;
usage: 
       $ProgramName [options] input.mnc output_prefix -mask <mask.mnc> -image <classified.mnc>
       $ProgramName -help to list options

USAGE

Getopt::Tabular::SetHelp ($help, $usage);

my $vol_classified = undef; 
my $vol_curvature = undef;
my $vol_mask = undef;
my $use_subcort = 1;
my $zip_files = 0;
my $clobber = 0;

my $sc_volume = "xPKGDATADIRx/gipt.mnc";
my $pve;
my @pve_ops;

my @options = (["-mask", "string", 1, \$vol_mask,
		"specify a mask volume (required)", "<mask.mnc>"],
	       ["-curve", "string",1, \$vol_curvature,
		"specify a curvature volume (optional)","<curvature.mnc>"],
	       ["-image", "string",1, \$vol_classified,
		"specify a classified volume (required)","<classified.mnc>"],
	       ['-subcortical|-nosubcortical','boolean', 1, \$use_subcort,
		'use subcortical tissue class [default]'],
	       ['-zip|-nozip','boolean',1, \$zip_files,
		'gzip the outputted files [default: -nozip]'],
	       ['-clobber|-noclobber','boolean',1,\$clobber,
		'clobber output files if the exist [default: -noclobber]']
	      );

GetOptions( \@options, \@ARGV ) || exit 1;

RegisterPrograms(['pve','pve3','gzip']);

if ( @ARGV == 2 ){
  $pve = "pve";
}
elsif ( @ARGV == 4) {
  $pve = "pve3";
}
else {
  die "$usage";
}

if (!defined($vol_mask) || !defined($vol_classified)) {
  die "$usage";
}

my $outpref = pop @ARGV;
my @invols = @ARGV;
my @outfiles = ("${outpref}_disc.mnc","${outpref}_wm.mnc","${outpref}_gm.mnc","${outpref}_csf.mnc");

if (defined($vol_curvature)) {
  push @pve_ops,"-curve";
  push @pve_ops,$vol_curvature;
}

if ($use_subcort) {
  push @pve_ops,"-subcortical";
  push @pve_ops,$sc_volume;
  push @outfiles,"${outpref}_sc.mnc";
}

if ($clobber) {
  push @pve_ops,"-clobber";
}


Spawn([$pve,@invols,$outpref,"-mask",$vol_mask,"-image",$vol_classified,@pve_ops]);

if ($zip_files) {

  Spawn(["gzip",@outfiles]);

}
