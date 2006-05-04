#!xPERLx -w

use strict;

use MNI::Startup;
use MNI::Spawn;
use Getopt::Tabular;
use File::Basename; 

my $help = <<HELP;

This is 'pve_4classes.pl'. It create a 4 class tissue classification,
combining sub-cortical gray into white. It is essentially a wrapper 
around the minccalc program in order to call it from CIVET. (PMP cannot
handle quotes in the minccalc expression while doing an echo of the
command string. Silly, but true!)

HELP

my $usage = <<USAGE;
usage: 
       $ProgramName [options] cls_correct.mnc cls_4classes.mnc
       $ProgramName -help to list options

USAGE

Getopt::Tabular::SetHelp ($help, $usage);

RegisterPrograms(['minccalc']);

if ( @ARGV != 2 ){
  die "$usage";
}

my $cls_correct = shift @ARGV;
my $cls_4classes = shift @ARGV;

my $expression = 'if (A[0] > 3.5 && A[0] < 4.5) { out = 2; } else { out = A[0]; }';
Spawn( [ "minccalc", "-clobber", "-expression", $expression, $cls_correct, $cls_4classes ] );

