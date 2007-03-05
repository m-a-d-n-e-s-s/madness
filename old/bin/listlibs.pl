
@libraries = ();
@includes = ();
%defines = ();
%read_files = ();

$debug = 0;

$includes[++$#includes] = ".";
$filename = "";

foreach $arg (@ARGV) {
    if ($arg =~ /^-d$/) {
        $debug = 1;
    }
    elsif ($arg =~ /^-D(.*)$/) {
        local($def) = $1;
        local($symbol) = $1;
        $def =~ s/^.*=//;
        $symbol =~ s/=.*$//;
        $defines{$symbol} = $def;
    }
    elsif ($arg =~ /^-I(.*)$/) {
        $includes[++$#includes] = $1;
    }
    else {
        $filename = $arg;
    }
}

if ($filename eq "") {
    print STDERR "listlibs.pl: require a filename\n";
    exit 1;
}

if (-d "$filename") {
    %current_includes = ();
    @libraries = ();
    %known_libs = ();
    %known_includes = ();
    &process_directory($filename);
}
else {
    &process_file($filename);
    %current_includes = ();
    @libraries = ();
    %known_libs = ();
    %known_includes = ();
    &find_libraries($filename);
}

@libraries = reverse(@libraries);

print "got $#libraries of them\n" if ($debug);

&substitute_defines();

foreach $i (0..$#libraries) {
    printf "%s", $libraries[$i];
    if ($i < $#libraries) { printf " "; }
}
printf "\n";

###########################################################################

sub process_file {
    local($filename) = shift;
    if ($debug) {
        printf "process_file: filename: %s\n", $filename;
    }

    # find the file
    local($ifile) = "";
    if ($filename =~ /^\//) {
        $ifile = $filename;
    }
    else {
        foreach $include (@includes) {
            $ifile = "$include/$filename";
            if ($debug) {
                #printf "process_file: looking for: %s\n", $ifile;
            }
            if (-f $ifile) { last; }
        }
    }
    if ($ifile eq "" || ! -f $ifile) {
        print STDERR "listlibs.pl: couldn't find file $ifile\n";
        exit 1;
    }

    # read the file
    local($filecontents) = "";
    open(IFILE,"<$ifile");
    while (<IFILE>) {
        if (/^\s*$/) { next; }
        $filecontents = "$filecontents$_";
    }
    close(IFILE);
    $read_files{$filename} = $filecontents;
    # an empty file will look like a new file below so put in a newline
    if ($read_files{$filename} eq "") {
        $read_files{$filename} = "\n"
    }

    # read in other files referenced by this file
    foreach $line (&get_lines($filecontents)) {
        if ($line =~ /^\#\s*include\s*<(.+)>/) {
            local($newfile) = $1;
            if ($read_files{$newfile} eq "") {
                &process_file($newfile);
            }
        }
    }
}

sub get_lines {
    local($filecontents) = shift;
    local(@lines) = ();
    local($ifdepth) = 0;
    while ($filecontents ne "") {
        # get next line
        $filecontents =~ s/^(.*)\n//;
        local($line) = $1;
        # remove comments
        $line =~ s/\/\/.*$//;
        # remove leading trailing whitespace
        $line =~ s/^\s*//;
        $line =~ s/\s*$//;
        # this only handles ifdef's that are one level deep
        if ($line =~ /\#\s*ifdef\s+([a-zA-Z_]\w*)/) {
            local($symbol) = $1;
            if (! defined $defines{$symbol}) {
                while ($filecontents ne "") {
                    $filecontents =~ s/^(.*)\n//;
                    local($tline) = $1;
                    last if ($tline =~ /\#\s*endif/);
                }
            }
        }
        elsif ($line =~ /\#\s*endif/) {
            # eat this endif
        }
        else {
            $lines[++$#lines] = $line if ($line ne "");
        }
    }
    return @lines;
}

sub find_libraries {
    local($filename) = shift;
    if ($current_includes{$filename} == 1) {
        print STDERR "listlibs.pl: recursive include detected for $filename\n";
        exit 1;
    }
    $current_includes{$filename} = 1;
    foreach $line (reverse(&get_lines($read_files{$filename}))) {
        if ($line =~ /^\#\s*include\s*<(.+)>/) {
            local($newfile) = $1;
            if ($known_includes{$newfile} != 1) {
                $known_includes{$newfile} = 1;
                &find_libraries($newfile);
            }
        }
        elsif ($known_libs{$line} != 1) {
            $known_libs{$line} = 1;
            $libraries[++$#libraries] = $line;
        }
    }
    delete $current_includes{$filename};
}

sub substitute_defines {
    local($i);
    local($symbol);
    foreach $i (0..$#libraries) {
        foreach $symbol (keys(%defines)) {
            $libraries[$i] =~ s/$symbol/$defines{$symbol}/g;
        }
    }
}

sub process_directory {
    local ($dir) = shift;
    opendir(DIR,"$dir");
    local (@files) = readdir(DIR);
    closedir(DIR);
    local ($i);
    foreach $i (@files) {
        if ("$i" eq "." || "$i" eq "..") {
            # skip
        }
        elsif (-d "$dir/$i") {
            process_directory("$dir/$i");
        }
        elsif ("$i" eq "LIBS.h") {
            process_file("$dir/$i");
            &find_libraries("$dir/$i");
        }
    }
}
