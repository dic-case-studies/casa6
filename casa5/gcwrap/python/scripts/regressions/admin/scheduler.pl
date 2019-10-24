#!/usr/bin/perl -sW
#
# Run a test when load average is low
#
# Note:
#    WARNING: The given workdir is recursively deleted!
#    $SHELL should be /bin/bash
#    If $DISPLAY is defined, that X connection is used
#    If $DISPLAY is undefined, a virtual Xvfb server is
#    started and stopped

# Distribution: END_USER_BINARY

use FindBin;
use Sys::Hostname;

###$#ARGV == 3 or die "Usage: $0 [options] reg_dir data_dir load_limit timeout";
$#ARGV == 1 or die "Usage: $0 [options] load_limit timeout";

$la_limit = $ARGV[0];
$timeout  = $ARGV[1];

###
### Avoid perl error for command line options
###
$pack_all = "" if (0);   # create a result-all.tgz from all the tar files?

if ( ! $prefix ) {
    $prefix = "/casa/regression/01";
}
if ( ! $work_dir ) {
    $work_dir = "$prefix/work/";
}
if (! $reg_dir) {
    $reg_dir = "$prefix/result/";
}
if (! $res_dir) {
    $res_dir = "$prefix/sresult/";
}

if (! $install_dir) {
    $install_dir = "$prefix/install/";
}
if (! $data_dir) {
    $data_dir = "$prefix/data/";
}

(-d $res_dir) or mkdir($res_dir) or die;

$cmd="date -u +%Y%m%d%H%M";
$date=`$cmd`;
if ($?) {
    print STDERR "$cmd: $! ". "bailing out";
    exit 1;
}
chomp($date);

$admin_dir = $FindBin::Bin;
$hostname = hostname();
$file_next  = $reg_dir . "/tests_next.txt";
$file_tests = $reg_dir . "/tests_list.txt";
(-d $reg_dir) or mkdir($reg_dir) or die "$reg_dir: $!";
(-e $file_tests) or $file_tests = $admin_dir . "/tests_list.txt";

if ($all) {
    unlink($file_next);
}

repeat:

open FILE, "<$file_tests" or die $! . $file_tests;

@tests=();
while(<FILE>) {
    s/^\s*//;
    s/\s*$//;
    s/\#.*$//;
    if (/^\S+$/) {
	push @tests, $_;
    }
    elsif (/^\S+\s+\S+\s+\S+$/) {
	@t = split;
	push @tests, $t[0];
    }
}

@tests = sort @tests;

#print $#tests+1, " tests:\n";
#print join("\n", @tests), "\n";

close FILE;
# force a run $la_limit=999;
`uptime 2>/dev/null` =~ /load averages?: (\S+),?/ or die $!;
$load_average_15 = $1;
chomp($load_average_15);
$load_average_15 =~ s/,/ /g;
print gettime();
printf "la %.2f ", $load_average_15;
printf "lim %.2f\n", $la_limit;

$reached_the_end = 0;
$executeproc = `ps -Af | egrep -v '^jmlarsen' | grep execute.py | grep -vw grep | wc -l`;
if ($load_average_15 < $la_limit) {
    if ($executeproc > 0) {
        print gettime(), "execute.py already running\n";
    }
    elsif (!(-d $data_dir)) {
        print gettime(), "Data directory $data_dir does not exist!\n";
    }
    else {
        # increment next
        (-e $file_next) or system("echo 0 0 > $file_next") == 0 or die $!;
        open FILE, "<$file_next" or die $!;
        ($next, $base) = split('\s', <FILE>);
        close FILE;
        # the file is expected to contain either
        #     <next>
        # or
        #     <next>  <base>
        # where <next> is the number of the next test to be run,
        # and <base> is an offset (defaulting to 0) to be added to <next>
        if (! $base) { $base = 0; }
	
	if (! $noloop && $next > $#tests) {
	    print STDERR gettime();
	    print STDERR
		"Warning: Next scheduled test is no. $next, " .
		"there are only " . ($#tests+1) . " tests in $file_tests. " .
		"Reset counter to 0\n";
	    $next = 0;
	}

	if ($next <= $#tests) {
 
	    open FILE, ">$file_next" or die;
            if (! $noloop) {
                print FILE (($next+1) % ($#tests+1)) . " " . $base . "\n";
            } else {
                print FILE ($next+1) . " " . $base . "\n";
            }
            close FILE;


            $test_number = ($next + $base) % ($#tests+1);	    
	    $testname = $tests[$test_number];
	    $runcasa_log     = "$work_dir/Log/run-$testname-$hostname.log";
	    
	    system("/bin/rm -rf $work_dir") == 0 or die $!;
	    mkdir("$work_dir") or die;
	    mkdir("$work_dir/admin") or die;
	    mkdir("$work_dir/Log") or die;
	    chdir("$work_dir") or die;    
	    system("cp $admin_dir/*.py $admin_dir/*.pl $work_dir/admin/") == 0 or die $!;
	    
	    $xdisplay = 7;
	    $p = "0";
	    if ($profile) {
		$p = "1";
	    }
	    $cmd = 
		"$admin_dir/process_manager.pl $timeout " .
		"$admin_dir/runcasa_from_shell.sh $xdisplay $work_dir/admin/execute.py $data_dir $work_dir $testname $p";
	    
	    print gettime(), "Run test $test_number $testname: $cmd > $runcasa_log 2>&1\n";
	    system("$cmd > $runcasa_log 2>&1"); # no check on return code
            if ($profile) {
		# be sure to stop the oprofile deamon if the casapy session did not
                system("sudo opcontrol --stop") == 0 or die $!;
            }
	    if (-d "$work_dir/result") {
		rename("$work_dir/result/", "$work_dir/Result/") or die $!;
		$result_files = "Result";
	    }
	    else {
		$result_files = "";
	    }
	    chdir("$work_dir") or die;
	    # can no longer write to logfile from here
	    system("/bin/rm -rf $work_dir/work") == 0 or die $!;
	    if (! $noclean) {
		$cmd = "tar zc $result_files Log/ > $res_dir/result-$hostname$date.tar.gz";
		if (system($cmd) != 0) {
		    print STDERR "$cmd: $!";
		}
		else {
		    print "Created $res_dir/result-$hostname$date.tar.gz\n";
		}
	    }
	    else {
		$cmd = "cp -R $result_files Log $res_dir";
		if (system($cmd) != 0) {
		    print STDERR "$cmd: $!";
		}
		else {
		    print "Updated $res_dir\n";
		}
	    }

	    if (! $noclean) {		
		# remove old results  (see get_results.pl)

		chdir("$res_dir") or die;
		if (-e "./retrieved") {
		    system("/bin/ls -1tr . | grep -B9999 retrieved | grep tar.gz | xargs /bin/rm -f");
		}
            }
            if ($pack_all) {
		$cmd = "tar zc result-*.tar.gz > result-all.tgz_temp";
		if (system($cmd) != 0) {
		    print STDERR "$cmd: $!";
		}
		else {
		    rename("result-all.tgz_temp", "result-all.tgz");
		    print "Updated $res_dir/result-all.tgz\n";
		}
	    }

	} # endif   next < $#tests
	else {
	    $reached_the_end = 1;
	}
    }
}
if ($all and !$reached_the_end)
{
    goto repeat; 
}


sub gettime
{
    ($sec,$min,$hour,$mday,$mon,$year,$wday,  $yday,$isdst)=localtime(time);
    $wday=$isdst=$yday=0;

    $ss = sprintf "%4d-%02d-%02d %02d:%02d:%02d -- ",  $year+1900,$mon+1,$mday,$hour,$min,$sec;
    return $ss;
}

exit 0;
