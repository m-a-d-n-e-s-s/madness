#!/usr/bin/perl

#
#  This file is part of MADNESS.
#
#  Copyright (C) 2012 Jinmei Zhang, Edward Valeev
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
#  $Id$
#

&usage() if $#ARGV == -1;

my %func_counts  = ();
my %func_tottime = ();
my %func_lattime = ();

my @funcs;
my $run_time     = 0;
my $tot_tasktime = 0;
my %have_thread     = {};

# to produce time-histograms set this to 1 and adjust the number of timeslices, if needed
my $make_histo        = 1;
my $nslices           = 100;
my $skip_shorter_than =
  0.00;    # skip if task shorter than this fraction of the timeslice

# use this search pattern to only include particular tasks in the running histogram
my %runninghisto_includes = (
  "all" => ".*"   # account all tasks 
# add additional patterns below, e.g. for ta_dgemm these are useful
  , "contract" => "contract"
  , "bcasth" => "bcast_.*_handler"
  , "bcastt" => "BcastTask"
);

foreach $file (@ARGV) {

	# Open file
	open( INFILE, "<$file" ) || die "cannot open file $file";

	# Collect data for each task in the file
	while (<INFILE>) {
		s/\n//g;
		my @line = split("\t");

		# get task data
		my $thread        = $line[0];    # Thread that the task ran on
		my $func_address  = $line[1];    # Address of the task function
		my $func_name     = $line[2];    # Name of task function
		my $func_nthreads = $line[3];    # Number of threads used by the task
		my $func_submit   = $line[4];    # Task submit time
		my $func_start    = $line[5];    # Task start time
		my $func_finish   = $line[6];    # Task finish time

        # mark presence of this thread
        $have_thread{$thread} = 1;

		# Calculate task run time and delay time
		my $func_time      = ( $func_finish - $func_start ) * $func_nthreads;
		my $func_delaytime = $func_start - $func_submit;

		# Accumulate function statistics
		if ( $func_counts{$func_name} == 0 ) {
			push( @funcs, $func_name );

			$func_counts{$func_name}  = 1;
			$func_tottime{$func_name} = $func_time;
			$func_lattime{$func_name} = $func_delaytime;

		}
		else {
			++$func_counts{$func_name};
			$func_tottime{$func_name} += $func_time;
			$func_lattime{$func_name} += $func_delaytime;
		}

		# Get the run time
		if ( $run_time < $func_finish ) {
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
printf "%-5s %7s %7s %10s  %5s   %10s\n", "name", "counts", "total", "average",
  "%", "latency";

foreach $func ( sort by_time ( keys(%func_tottime) ) ) {
	my $counts  = $func_counts{$func};
	my $tottime = $func_tottime{$func};

	my $avgtime  = $tottime / $counts;
	my $lattime  = $func_lattime{$func} / $counts;
	my $perctime = $tottime / $run_time * 100;
	printf "%s\n %10d  %10.7f %10.7f %5d   %10.7f\n", $func, $counts, $tottime,
	  $avgtime, $perctime, $lattime;
}

print "\n";

#
# time histograms
#
if ($make_histo) {
	use POSIX;
	my $timestep = $run_time / $nslices;

	# time-weighted
	my @pendinghisto = ();
	foreach $ts ( 0 .. $nslices - 1 ) { push @pendinghisto, 0.0; }
	my @inprogresshisto = ();
	foreach $ts ( 0 .. $nslices - 1 ) { push @inprogresshisto, 0.0; }
	my %runninghistos = {};
    foreach my $key (keys %runninghisto_includes) {
    	my $values = [];
        foreach $ts ( 0 .. $nslices - 1 ) { push @$values, 0.0; }
        $runninghistos{$key} = $values;
    }

	# not time-weighted
	my @spawnedhisto = ();
	foreach $ts ( 0 .. $nslices - 1 ) { push @spawnedhisto, 0.0; }
	my @startedhisto = ();
	foreach $ts ( 0 .. $nslices - 1 ) { push @startedhisto, 0.0; }
	my @retiredhisto = ();
	foreach $ts ( 0 .. $nslices - 1 ) { push @retiredhisto, 0.0; }

    # this is the list of tasks active for each time slice
    # each entry is (ref to) an array of tasks active (not necessarily actively running, i.e. could be waiting for
    # subtasks to finish)
    my @inprogresstasks = ();
    foreach my $slice (0 .. $nslices-1) {
    	$inprogresstasks[$slice] = [];
    }

	foreach $file (@ARGV) {

		# Open file
		open( INFILE, "<$file" ) || die "cannot open file $file";

		# this will keep running tasks for each thread for deeper processing
		my %runningtasks = {};
		foreach my $threadid (keys %have_thread) {
			$runningtasks{$threadid} = [];
		}

		# Collect data for each task in the file
		while (<INFILE>) {
			s/\n//g;
			my @line = split("\t");

            my $thread        = $line[0];    # Thread that the task ran on
            my $func_name     = $line[2];    # Name of task function
			my $func_nthreads = $line[3];   # Number of threads used by the task
#			die("don't know how to handle a multithreaded task")
#			  if ( $func_nthreads != 1 );
			my $time_submit = $line[4];     # Task submit time
			my $time_start  = $line[5];     # Task start time
			my $time_finish = $line[6];     # Task finish time

            my $tasklist = $runningtasks{$thread};
			my $task = { "start" => $time_start, "finish" => $time_finish,
				         "name" => $func_name, "thread" => $thread };
			push @$tasklist, $task;

			if ( $time_start - $time_submit > $timestep * $skip_shorter_than ) {
				my $submit_ts_0 = floor( $time_submit / $timestep );
				my $submit_ts_1 = floor( $time_start / $timestep );
				if ( $submit_ts_0 != $submit_ts_1 ) {
					$pendinghisto[$submit_ts_0] +=
					  ( $submit_ts_0 + 1 ) - ( $time_submit / $timestep );
					foreach my $ts ( $submit_ts_0 + 1 .. $submit_ts_1 - 1 ) {
						$pendinghisto[$ts] += 1;
					}
					$pendinghisto[$submit_ts_1] +=
					  ( $time_start / $timestep ) - $submit_ts_1;
				}
				else {
					$pendinghisto[$submit_ts_0] +=
					  ( $time_start - $time_submit ) / $timestep;
				}
			}

			if ( $time_finish - $time_start > $timestep * $skip_shorter_than ) {
				my $start_ts_0 = floor( $time_start / $timestep );
				my $start_ts_1 = floor( $time_finish / $timestep );
				if ( $start_ts_0 != $start_ts_1 ) {
					$inprogresshisto[$start_ts_0] +=
					  ( $start_ts_0 + 1 ) - ( $time_start / $timestep );
					push @{ $$inprogresstasks[$start_ts_0] }, $task;
					foreach my $ts ( $start_ts_0 + 1 .. $start_ts_1 - 1 ) {
						$inprogresshisto[$ts] += 1;
                        push @{ $$inprogresstasks[$ts] }, $task;
					}
					$inprogresshisto[$start_ts_1] +=
					  ( $time_finish / $timestep ) - $start_ts_1;
                        push @{ $$inprogresstasks[$start_ts_1] }, $task;
				}
				else {
					$inprogresshisto[$start_ts_0] +=
					  ( $time_finish - $time_start ) / $timestep;
                    push @{ $$inprogresstasks[$start_ts_0] }, $task;
				}
			}

			my $spawned_ts = floor( $time_submit / $timestep );
			$spawnedhisto[$spawned_ts] += 1;
			my $started_ts = floor( $time_start / $timestep );
			$startedhisto[$started_ts] += 1;
			my $stop_ts = floor( $time_finish / $timestep );
			$retiredhisto[$stop_ts] += 1;
		}
		close INFILE;

		# compute the time-averaged # of running
		foreach my $threadid (keys %runningtasks) {
#			printf STDOUT "threadid = %d\n", $threadid;
			my $tasklistref = $runningtasks{$threadid};
			my @tasklist = sort { $$a{"start"} <=> $$b{"start"} } @$tasklistref;
				                  
		    my $time_finish_last = 0.0;
			foreach my $taskptr (@tasklist) {
				my %task = %$taskptr;
                my $time_start = $task{"start"};
                my $time_finish = $task{"finish"};
                my $func_name = $task{"name"};
#				printf STDOUT "%lf %lf\n", $task{"start"}, $task{"finish"};

                if ($time_finish > $time_finish_last) {
                	
                	if ($time_start < $time_finish_last) {
                		$time_start = $time_finish_last;
                	}
                	
                    my $start_ts_0 = floor( $time_start / $timestep );
                    my $start_ts_1 = floor( $time_finish / $timestep );
                    
                    if ( $start_ts_0 != $start_ts_1 ) {
                    	foreach my $key (keys %runninghisto_includes) {
                    		my $func_name_pattern = $runninghisto_includes{$key};
                            if ($func_name =~ /$func_name_pattern/) {
                              ${ $runninghistos{$key} }[$start_ts_0] +=
                                  ( $start_ts_0 + 1 ) - ( $time_start / $timestep );
                              foreach my $ts ( $start_ts_0 + 1 .. $start_ts_1 - 1 ) {
                                  ${ $runninghistos{$key} }[$ts] += 1;
                            }
                            ${ $runninghistos{$key} }[$start_ts_1] +=
                              ( $time_finish / $timestep ) - $start_ts_1;
                            }                    		
                    	}
                    }
                    else {
                        foreach my $key (keys %runninghisto_includes) {
                            my $func_name_pattern = $runninghisto_includes{$key};
                            if ($func_name =~ /$func_name_pattern/) {
                                ${ $runninghistos{$key} }[$start_ts_0] +=
                                  ( $time_finish - $time_start ) / $timestep;
                         	}
                        }
                    }
                	
                	$time_finish_last = $time_finish;
                }
			}
		}

		# print out the histograms for this trace
		printf STDOUT "file = %s\n",      $file;
		printf STDOUT "timestep = %lf\n", $timestep;

		printf STDOUT "# ts nspawned nstarted nretired\n";
		foreach my $ts ( 0 .. $nslices - 1 ) {
			printf STDOUT "%d\t%lf\t%lf\t%lf\n", $ts, $spawnedhisto[$ts], $startedhisto[$ts], $retiredhisto[$ts];
		}
#		printf STDOUT "time-weighted # of waiting\n";
#		foreach my $ts ( 0 .. $nslices - 1 ) {
#			printf STDOUT "%d %lf\n", $ts, $pendinghisto[$ts];
#		}
#		printf STDOUT "time-weighted # of in-progress\n";
#		foreach my $ts ( 0 .. $nslices - 1 ) {
#			printf STDOUT "%d %lf\n", $ts, $inprogresshisto[$ts];
#		}

		# list keys in a reasonable order (keys returns keys in a weird order)
		# "all" first, then the rest
		my @running_keys = ("all");
		foreach my $key (keys %runninghisto_includes) {
			push @running_keys, $key if $key ne "all";
		}
        printf STDOUT "# ts, pend";		
        foreach my $key (@running_keys) {
        	printf STDOUT ", run_%s", $key;
        }
        printf STDOUT "\n";
		foreach my $ts ( 0 .. $nslices - 1 ) {
			printf STDOUT "%d", $ts;
			printf STDOUT ", %lf", $pendinghisto[$ts];
			foreach my $key (@running_keys) {
  			    printf STDOUT ", %lf", ${ $runninghistos{$key} }[$ts];
			}
			printf STDOUT "\n";
		}
		
		my $active_tasks_fname = "$file.active_tasks.log";
		open ATFILE, ">$active_tasks_fname" || die "could not open $active_tasks_fname";
		foreach my $ts (0 .. $nslices-1) {
			printf ATFILE "============================ ts %d: [%lf,%lf] =============================\n",
			  $ts, $ts*$timestep, ($ts+1)*$timestep;
			foreach my $taskref (sort { $$a{"thread"} <=> $$b{"thread"} } @{ $$inprogresstasks[$ts] }) {
				printf ATFILE "%d %s %lf %lf\n",
				  $$taskref{"thread"}, $$taskref{"name"},
				  $$taskref{"start"}, $$taskref{"finish"};
			}
		}
		close ATFILE;

	}

}

exit 0;

sub usage {
	printf STDERR
	  "taskprofile.pl converts trace files info a summary profile. Usage:\n";
	printf STDERR
	  "  taskprofile.pl <trace file1> [<trace file2> <trace file3>... ]\n";
	exit 0;
}

sub by_time {
	$func_tottime{$b} <=> $func_tottime{$a};
}
