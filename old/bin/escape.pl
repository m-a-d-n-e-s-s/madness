
# read standard input, escape objectionable characters, and print to standard out
while(<STDIN>) {
  printf STDOUT escape($_);
}

exit 0;

sub escape {
  local($line) = shift;
  $line =~ s/ /\\ /g;

  return $line;
}

