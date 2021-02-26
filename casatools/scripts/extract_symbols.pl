use File::Basename;
use File::Find;
use File::Path 2.07 qw( make_path remove_tree );
use File::Spec;
use Cwd qw(abs_path getcwd);
use List::Util qw(max);
use POSIX qw(uname);
use File::Copy;
use Term::ANSIColor;
#use URI::Split qw(uri_split);
use feature 'state';
use File::Copy qw(copy);
use File::Copy qw(move);
use Cwd qw(abs_path);

foreach ( @ARGV ) {
    m|^ws=(\S+)|i && -d "$1" && ( $workspace = $1, next );
    m|^binaries=(\S+)|i && ( $binaries = $1, next );
    m|^version=(\S+)|i && ( $version = $1, next );
}

process_debug_symbols($binaries, $workspace);

sub process_debug_symbols {
    
    my $appdir = shift( @_ );
    my $workspace = shift( @_ );
    print "Processing debug symbols\n";
    print "appdir: $appdir\n";
    print "workspace: $workspace\n";
    print "version: $version\n";

    $osname = lc((uname( ))[0]);

    my $script_dir = dirname(abs_path(__FILE__));
    print("Script dir: $script_dir\n");
    
    print("$appdir\n"); 

    $dump_syms_link="https://casa.nrao.edu/download/devel/dump_syms/$osname/dump_syms";
    print("Getting dump_syms from: $dump_syms_link");
    `curl -o  $workspace/dump_syms -O -L https://casa.nrao.edu/download/devel/dump_syms/$osname/dump_syms`;
    $dump_syms="$workspace/dump_syms";
    print "dump_sym cmd: $dump_syms\n";
    $mode = 0744;
    chmod $mode, $dump_syms;  

    $symbols_dir_path = $workspace;
    $symbols_dir = "casa-$version-symbols";
    $syms_zipfile = "casa-$version-debugsymbols.zip";
    print("Scanning symbols\n");
    
    scan_syms($appdir, "$symbols_dir_path/$symbols_dir", $dump_syms);
    print("Packing symbols\n");
    pack_syms($symbols_dir_path, $symbols_dir, $syms_zipfile);
    
}

sub pack_syms {
    my $symbols_dir_path = shift( @_ );
    my $dir = shift( @_ );
    my $syms_zipfile = shift( @_ );
    print "Packaging directory $dir to package $syms_zipfile\n";
    print "symbols_dir_path: $symbols_dir_path";
    `cd $symbols_dir_path; zip -r $syms_zipfile $dir; cd -`
}

sub submit_syms {
    my $syms_zipfile = shift( @_ );
    my $syms_destination = shift( @_ );
    print "Sending $symbols_dir_path/$syms_zipfile to $syms_destination\n";
}

sub scan_syms {
    my $dir = shift( @_ );
    my $symbols_dir = shift( @_ );
    my $dump_syms = shift( @_ );
    my @executables = ( );
    print "scan_syms dir: $dir\n";

    my $find_func = sub {  my $dir = $File::Find::dir;
                           ## skip subversion files, cmake-created applications, and cmake build-tree files...
                           #   if ( $dir =~ m|/\.svn| || $dir =~ m|/MacOS| || $dir =~ m|/CMakeFiles| ) { return }
                           ## skip subversion files and cmake build-tree files... (scanning our application tree)
                           if ( $dir =~ m|/\.svn| || $dir =~ m|/CMakeFiles| ) { return }
                           if ( -f $_ ) {
                               my $file_type = resolve_file_type($_);
                               $file_type eq 'exe' && push( @executables, "$File::Find::dir/$_" );
                               $file_type eq 'lib' && push( @executables, "$File::Find::dir/$_" );
                               $file_type eq 'ldso' && push( @executables, "$File::Find::dir/$_" );
                           }
                        };

    find( { wanted => $find_func }, $dir );
    #unless ( $quiet ) { print "\tfound ",sprintf($spf,scalar(@executables))," binaries\n" }
    mkdir ($symbols_dir);
    foreach ( @executables ) {
        print("$_\n");
        if ( $osname eq "linux" ) {
          `readelf -S $_ | grep "no sections in this file"`;
          if ($? == 0) {
            print("$_ does not have sections. This will crash dump_syms Skipping...\n");
            next;
          }
        }
        $fileheader = `$dump_syms $_ | head -n1`;
        chomp($fileheader);
        print "fileheader: $fileheader\n";
        $filename = (split / /, $fileheader)[4] ;
        print "filename: $filename\n";
        $filehash = (split / /, $fileheader)[3];
        print "filehash: $filehash\n";
        print "Creating directory $symbols_dir/$filename\n";
        mkdir ("$symbols_dir/$filename");
        print "Creating directory $symbols_dir/$filename/$filehash\n";
        mkdir ("$symbols_dir/$filename/$filehash");
        print "Writing headers to:";
        $symfile = "$symbols_dir/$filename/$filehash/$filename.sym";
        print "$symfile\n";
        $dumpcmd = "$dump_syms $_ > $symfile";
        print "Dump cmd: $dumpcmd\n";
        `$dumpcmd`;
        print "Stripping $_\n";
        `strip -S $_`;
    }
    
    # Find cmake object files
    my $find_func2 = sub {  my $dir = $File::Find::dir;
                           ## skip subversion files.
                           if ( $dir =~ m|/\.svn| ) { return }
                           if ( -f $_ ) {
                               my $file_type = resolve_file_type($_);
                               $file_type eq 'cmakeobjectfile' && push( @executables, "$File::Find::dir/$_" );
                               $file_type eq 'exe' && push( @executables, "$File::Find::dir/$_" );
                               $file_type eq 'lib' && push( @executables, "$File::Find::dir/$_" );
                               $file_type eq 'ldso' && push( @executables, "$File::Find::dir/$_" );
                           }
                        };
    # Strip cmake object files
    @executables = ( );
    $ws="$dir/../";
    print "Stripping files in $ws. Including the files in cmake directories\n";
    find( { wanted => $find_func2 }, $ws );
    foreach ( @executables ) {
        print "Stripping $_\n";
        `strip -S $_`;
    }

}


sub resolve_file_type {
    state %cache;
    my $type = '';
    my $f = shift(@_);
    if ( $osname eq "linux" && -l $f ) {
        ### internal symlinks (like all of the mpi binaries) should be disregarded
        ### OSX looks for python library linked to Python (i.e. a symlink)... so
        ### we only do this break for linux (there's probably a better way)
        return "?";
    }
    if ( -f $f ) {
        my $file = abs_path($f);
        return $cache{$file} if exists $cache{$file};
        $type = '?';
        # assign type based upon filename
        $f =~ m|\.h(?:pp)?$| && ($type = 'h');
        $f =~ m@\.c(?:pp|c)?$@ && ($type = 'c');
        $f =~ m|\.f$| && ($type = 'f');
        $f =~ m|\.cc\.o$| && ($type = 'cmakeobjectfile');
        $f eq 'CMakeLists.txt' && ($type = 'cmake');
        $f =~ m@^(?:m|GNUm|M)akefile$@ && ($type = 'make');
        if ( $type eq '?' ) {
            # resort to using file utility...
            my $current_file = `file '$file'`;
            if ($current_file =~ m|Mach-O.*?executable|) { $type = 'exe'; }
            if ($current_file =~ m|Bourne-Again shell script|)  {$type = 'sh'; }
            if ($current_file =~ m|Mach-O 64-bit bundle|) { $type = 'ldso'; }
            if ($current_file =~ m|Mach-O 64-bit x86_64 bundle|) { $type = 'ldso'; }
            if ($current_file =~ m|archive random library|)  {$type = 'a'; }
            if ($current_file =~ m|ELF.*executable|) { $type = 'exe'; }
            if ( $current_file =~ m@dynamically linked shared library|ELF.*shared object@ ) {
                my $nme = basename($f);
                ### for osx framework libraries can look like 'Python' or 'QtCore' so
                ### this exception keys off of "*.so" for osx... for linux, the
                ### 'lib...' check is probably enough...
                $type =  $nme=~ m|^lib| ? 'lib' : $nme =~ m|\.so$| ? 'ldso' : 'lib';
            }
            ##
            ## we absolutely have to check, because CASA's build system seems to think
            ## that every python script should be marked as executable...
            ##
            if ( -x $f && $current_file =~ m|python.*?script(?:, ASCII)? text executable| ) {
                open( my $py, "< $f" );
                my $first_line = <$py>;
                close($py);
                if ( $first_line =~ m|^#!| ) { $type = 'py'; }
                ## should create remove_execute_permission(...)
                else { system( "chmod -x $f" ) == 0 or die "could not chmod, $f" } ## fix permissions
            }

        }
        $cache{$file} = $type;
    }
    return $type;
}
