#!/usr/bin/perl -w
use strict;
use Time::HiRes qw(usleep ualarm gettimeofday tv_interval);

my $old_fh = select(STDOUT);
$| = 1;
select($old_fh);

unless(scalar(@ARGV) >= 2){
  warn "Usage:  $0 <L> <N> [--CMD <dir/sfs_code>][-t 0.001] [-r 0] [-W 0] [-B 0.5] [-v]\n";
  warn "where:\tL = total desired locus length to simulate\n";
  warn "\tN = total ancestral population size to simulate\n";
  warn "\nOther parameters are optional and unordered, with defaults indicated above.  Update following SFS_CODE documentation.\n";
  warn "\t--CMD <dir/sfs_code> = Use this flag and specify location of sfs_code if not cwd\n";
  warn "\t-t theta = population-scaled mutation rate per base pair\n";
  warn "\t-r rho = population-scaled recombination rate per base pair\n";
  warn "\t-W <args> = selection parameters (see SFS_CODE documentation).\n";
  warn "\t-B burnin = the burn-in time for the optimization, only a small fraction of 2N gens needed.\n";
  warn "\t-v = verbose mode\n\n";

  warn "SFS_CODE is an efficient simulator, but when a long simulated sequence is desired (and/or large population 
size), breaking it up into smaller linked loci can sometimes be more efficient.  This program provides a 
heuristic method for optimizing the locus length.  Because a majority of a forward simulation's computation 
time is usually spent during the burn-in, optimizeLL.pl seeks to reduce this computational burden by 
finding a combination of linked loci that reduces computational time.  optimizeLL.pl focuses on the 
initial 0.5*2N generations, because by this point it usually evident which combination will be most 
efficient.  However, if you are doing a lot of simulations, it is recommended that you run at least 
10 iterations of the full-scale simulation across a range of different locus lengths and compare their 
average run times.  optimizeLL.pl should only act as a rough guide.\n\n";
  exit;
}

#PARTS is an array of the number of loci to consider.  Only evaluate until the run time gets worse.
my @PARTS = (1, 2, 5, 10, 20, 50, 100, 200, 500, 1000, 2000, 5000, 10000, 20000, 50000);

#define default parameter values
my $L = $ARGV[0];
my $N = $ARGV[1];
my $CMD = "./sfs_code";
my $THETA = 0.001;
my $RHO = 0;
my $W = "0";
my $B = 0.5;
my $VERBOSE = 0;
my $PRINTGEN = 0;

#parse command-line to get parameters and SFS_CODE options
for(my $i=2; $i<scalar(@ARGV); $i++){
  if($ARGV[$i] eq "--CMD" && scalar(@ARGV) >= $i+1){
    $CMD = $ARGV[$i+1];
    $i++;
  }
  elsif($ARGV[$i] eq "-t" && scalar(@ARGV) >= $i+1){
    $THETA = $ARGV[$i+1];
    $i++;
  }
  elsif($ARGV[$i] eq "-r" && scalar(@ARGV) >= $i+1){
    $RHO = $ARGV[$i+1];
    $i++;
  }
  elsif($ARGV[$i] eq "-B" && scalar(@ARGV) >= $i+1){
    $B = $ARGV[$i+1];
    $i++;
  }
  elsif($ARGV[$i] eq "-W"){
    if($ARGV[$i+1] == 0){
      $i++;
    }
    elsif($ARGV[$i+1] == 1 && scalar(@ARGV) >= $i+4){
      $W = "1 $ARGV[$i+2] $ARGV[$i+3] $ARGV[$i+4]";
      $i += 4;
    }
    elsif($ARGV[$i+1] == 2 && scalar(@ARGV) >= $i+6){
      $W = "2";
      for(my $j=2; $j<=6; $j++){
	$W .= " $ARGV[$i+$j]";
      }
      $i += 6;
    }
    elsif($ARGV[$i+1] == 3 && scalar(@ARGV) >= $i+3){
      $W = "3 $ARGV[$i+2] $ARGV[$i+3]";
      $i += 3;
    }
    else{
      die "error specifying selection types.  Please choose type 0-3, and check the arguments\n";
    }
  }
  elsif($ARGV[$i] eq "-v"){
    $VERBOSE = 1;
  }
  elsif($ARGV[$i] eq "--printGen"){
    $PRINTGEN = 1;
  }
  else{
    die "unrecognized switch \"$ARGV[$i]\", or missing parameters.\n";
  }
}

#this is the basic command to run
my $RUN = "$CMD 1 1 -N $N -t $THETA -r $RHO -W $W -B $B ";
if($VERBOSE){
  print "$RUN\n";
}
else{
  if($PRINTGEN){
    die "Can only use --printGen option in verbose mode (-v).\n";
  }
}
$RUN = $RUN." -o /dev/null"; #this is to suppress output

if($PRINTGEN){
  $RUN = $RUN." --printGen";
}
else{
  $RUN = $RUN." -e /dev/null";
}

#march through partitions until time gets worse, return optimal
my @TIMES = ();
my $Opart = 1;
my $Olength = $L;
my $Otime = 1e10;
for(my $i=0; $i<scalar(@PARTS); $i++){
  my $part = $PARTS[$i];
  my $length = int($L/$part);
  if($VERBOSE){
    print "-L $part $length: ";
  }
  my ($startS, $startM) = gettimeofday();
  system("$RUN -L $part $length");
  my ($stopS, $stopM) = gettimeofday();
  $TIMES[$i] = $stopS+$stopM/1e6-($startS+$startM/1e6);
  
  if($VERBOSE){
    print " $TIMES[$i]";
  }
  
  if($TIMES[$i] < $Otime){
    $Opart = $part;
    $Olength = $length;
    $Otime = $TIMES[$i];
    if($VERBOSE){
      print " **";
    }
  }
  if($VERBOSE){
    print "\n";
  }
  #check to see if the locus length is getting too small, or if there is consistent increase in runtime.
  if($length <= 100 || ($i>=5 && $TIMES[$i] > $TIMES[$i-1] && $TIMES[$i-1] > $TIMES[$i-2] && $TIMES[$i-2] > $TIMES[$i-3])){
    if($VERBOSE){
      print "~optimal locus configuration is:\n";
    }
    print "-L $Opart $Olength\n";
    last;
  }
}

exit(1);
