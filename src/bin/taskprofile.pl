#!/usr/bin/perl

#
#  This file is part of MADNESS.
#
#  Copyright (C) 2012 Jinmei Zhang
#
#  This program is free software; you can redistribute it and/or modify
#  it under the terms of the GNU General Public License as published by
#  the Free Software Foundation; either version 2 of the License, or
#  (at your option) any later version.
#
#  This program is distributed in the hope that it will be useful,
#  but WITHOUT ANY WARRANTY; without even the implied warranty of
#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
#  GNU General Public License for more details.
#
#  You should have received a copy of the GNU General Public License
#  along with this program; if not, write to the Free Software
#  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
#
#  For more information please contact:
#
#  Robert J. Harrison
#  Oak Ridge National Laboratory
#  One Bethel Valley Road
#  P.O. Box 2008, MS-6367
#
#  email: harrisonrj@ornl.gov
#  tel:   865-241-3937
#  fax:   865-572-0680
#
#  $Id $
#

&usage() if $#ARGV == -1;

my %func_counts = ();
my %func_tottime = ();
my %func_lattime = ();

my @funcs;
my $run_time = 0;
my $tot_tasktime = 0;
my $nthreads = -1;

foreach $file (@ARGV)
{
  # Open file
  open(INFILE, "$file") || die "cannot open file $file";
  
  # Collect data for each task in the file
  while (<INFILE>) {
    s/\n//;
    my @line = split("\t");
    
    # get task data
    my $thread = $line[0];        # Thread that the task ran on
    my $func_address = $line[1];  # Address of the task function
    my $func_name = $line[2];     # Name of task function
    my $func_nthreads = $line[3]; # Number of threads used by the task
    my $func_submit = $line[4];   # Task submit time
    my $func_start = $line[5];    # Task start time
    my $func_finish = $line[6];   # Task finish time
    #if ($max_nthread < $nthread) {
    #  $max_nthread = $nthread;
    #}
    #printf "%s\n", $func_name,
    
    # Calculate task run time and delay time
    my $func_time = ($func_finish - $func_start) * $func_nthreads;
    my $func_delaytime = $func_start - $func_submit;
    
    # Accumulate function statistics
    if ($func_counts{$func_name} == 0) {
      push (@funcs, $func_name);
      
      $func_counts{$func_name} = 1;
      $func_tottime{$func_name} = $func_time;
      $func_lattime{$func_name} = $func_delaytime;
      
    } else {
      ++$func_counts{$func_name};
      $func_tottime{$func_name} += $func_time;
      $func_lattime{$func_name} += $func_delaytime;
    }
    
    # Get the run time
    if ($run_time < $func_finish) {
      $run_time = $func_finish;
    }
    
    # Accumulate total function run time
    $tot_functime += $func_time;
  }
  close INFILE;
}

printf "\nWall time:       %10.6f (s)\n", $run_time;
printf "Total task time: %10.6f (s)\n\n", $tot_functime;

printf "%40s\n", "Time (s)";
printf "%-5s %7s %7s %10s  %5s   %10s\n", "name", "counts", "total", "average", "%", "latency";

foreach $func (sort by_time (keys(%func_tottime)))
{
  my $counts = $func_counts{$func};
  my $tottime = $func_tottime{$func};
  
  my $avgtime = $tottime / $counts;
  my $lattime = $func_lattime{$func} / $counts;
  my $perctime = $tottime / $run_time * 100;
  printf "%s\n %10d  %10.7f %10.7f %5d   %10.7f\n", $func, $counts, $tottime, $avgtime, $perctime, $lattime;
}

print "\n";
exit 0;

sub usage {
  printf STDERR "taskprofile.pl <input file(s)>\n";
  exit 0;
}

sub by_time {
  $func_tottime{$b} <=> $func_tottime{$a};
}



 
