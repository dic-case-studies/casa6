#!/usr/bin/perl
use strict;
use warnings;

(my $major, my $minor, my $micro) = split(/\./, `sw_vers -productVersion`);
my $os_ver=$major . "." .$minor;

# Read libs

my $solibdir = "./build/lib.macosx-" . $os_ver . "-x86_64-3.6/casatools/__casac__/";
print ("libdir: $solibdir\n");
my $dylibdir = $solibdir . "/lib/";

# Fix .so
my $t = localtime();
print("$t Fixing solibs\n");

opendir(DIR, $solibdir) or die "Could not open $solibdir\n";
while (my $filename = readdir(DIR)) {
    if ( $filename =~ /.so$/) {
      print ("$filename\n");
      fixlibpaths($solibdir. "/" .$filename,"/opt/casa/0.?", "\@rpath");
  }
}
closedir(DIR);

$t = localtime();
print("$t Fixing dylibs\n");
# Fix .dylib
opendir(DIR, $dylibdir) or die "Could not open $dylibdir\n";
while (my $filename = readdir(DIR)) {
  if ( $filename =~ /.dylib$/) {
      print ("$filename\n");
      fixlibpaths($dylibdir. "/" .$filename,"/opt/casa/0.?", "\@rpath");
  }
}
closedir(DIR);

sub fixlibpaths() {

    my ($filename,$old_path_ish, $new_path_ish) = @_;
    my @deps=`otool -L $filename`;
    #open(LIB,$deps) or die("unable to open build.conf");
        foreach (@deps) {
            my $oldpath= (split(" ", $_))[0];
            chomp($oldpath);
            $oldpath =~ s/^\s+//;
            $oldpath =~ s/://;
            #print($oldpath);
            $_ =~ s/$old_path_ish/$new_path_ish/;
            my $r1="\tlibcfitsio.dylib";
            my $r2="\@rpath/lib/libcfitsio.dylib";
            if ( $filename =~ /.dylib$/) {
                $r2="\@rpath/libcfitsio.dylib";
            }
            $_ =~ s/$r1/$r2/;
            my $newpath = (split(" ",$_))[0];
            chomp($newpath);
            $newpath =~ s/^\s+//;
            $newpath =~ s/://;
            #flatten libgcc
            $newpath =~ s/libgcc\///;
            #print ( $oldpath . " " . $newpath);
            print("install_name_tool -change $oldpath $newpath $filename\n");
            my $output = `install_name_tool -change $oldpath $newpath $filename`;
            print($output);
        }
    #close(LIB)
}
