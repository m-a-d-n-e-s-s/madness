# Emacs should use -*- Perl -*- mode.

# This script takes up to 2 parameters:
# 1) F77 name->symbol conversion descriptor, one of the following strings: symbol_, symbol, SYMBOL, SYMBOL_
# 2) (optional) F77 name to convert
# If second argument is given, script converts the F77 name, prints the corresponding F77 symbol
# to STDOUT and exits.
# If only one argument is given then the script reads names of F77 names from STDIN. For each F77 name XXX
# the output is '#define xxx_ <XXX symbol>' where <XXX symbol> is the F77 symbol which corresponds
# to name XXX. If <XXX symbol>==xxx_, however, the output is null. The output is written to STDOUT

die if ($#ARGV < 0 || $#ARGV > 1);

my $to_macro = ($#ARGV == 0);
my $f77symconv = $ARGV[0];

if (!$to_macro) {
  my $f77name = $ARGV[1];
  printf STDOUT name2symbol($f77symconv,$f77name);
}
else {
  while (<STDIN>) {
    chomp;
    s/\s*//g;
    next if ($_ eq "");
    my $f77symbol = name2symbol($f77symconv,$_);
    my $redef_f77symbol = name2symbol("symbol_",$_);
    if ($f77symbol ne $redef_f77symbol) {
      printf STDOUT "#define %s %s\n",$redef_f77symbol,$f77symbol;
    }
  }
}

exit 0;


sub name2symbol
{
  my ($f77symconv, $name) = @_;
  # by default, F77 symbol is same as F77 name
  my $symbol = $name;
  if ($f77symconv eq "symbol") {
    $symbol = lc($name);
  }
  elsif ($f77symconv eq "symbol_") {
    $symbol = lc($name)."_";
  }
  elsif ($f77symconv eq "SYMBOL") {
    $symbol = uc($name);
  }
  elsif ($f77symconv eq "SYMBOL_") {
    $symbol = uc($name)."_";
  }
  return $symbol;
}
