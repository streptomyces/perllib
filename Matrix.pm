package Sco::Matrix;
use 5.14.0;

# {{{ Blurb

=head1 Name

Micro::Matrix

=head1 Blurb

This is to manage and use scoring matrices.

=head1 Example

See the script matrix.pl in /home/sco/customers/nick_tucker/

=head1 See also

Put something here later.

=head1 History of changes

=over 2

=item Wed May 16 2007

First writing

=back

=head1 Author

Govind Chandra E<lt>govind.chandra@bbsrc.ac.ukE<gt>

=head1 Subroutines


=cut

# }}}

#{{{ sub new()

=head2 new

=head3 Arguments

=over 2

=item -file

A filename from which the matrix is read and created.

=back

=head3 Returns

A Micro::Matrix object

=cut

sub new {
my($class, $self, %argv);
($class, %argv)=@_;
$self={
'file' => $argv{-file}
};
$self->{matrix}=&_make_matrix($self->{file});
&_set_length($self);
bless($self, $class);
return($self);
}

# }}}

# {{{ sub _set_length

=head2 _set_length

=head3 Arguments

=over 2

=item 1

A Micro::Matrix object.

=back

=head3 Returns

Nothing

=head3 Note

This is an internal function. It is used to store
the length of the matrix in the Micro::Matrix object
during its creation in C<new()>.

=cut


sub _set_length {
my $self=shift(@_);
my $matrix = $self->{matrix};
my @keys=keys(%{$matrix});
my $key=shift(@keys);
my $matlen=scalar @{$matrix->{$key}};
$self->{length}=$matlen;
}

# }}}

# {{{ sub _make_matrix

=head2 _make_matrix

=head3 Arguments

=over 2

=item 1

A filename

=back

=head3 Returns

A reference to the matrix created.

=head3 Note

This is an internal function called from C<new()>.
This is the function which actually reads the matrix
file and creates the matrix.

=cut


sub _make_matrix {
my $file=shift(@_);
open(INFILE, "<$file");
my $matref=undef;
while(<INFILE>) {
my $line=$_;
chomp($line);
if($line=~m/^\s*\#/ or $_=~m/^\s*$/) {next;}
my @llist=split(/\s+/, $line);
my $symbol=uc(shift(@llist));
my $cnt=0;
foreach my $prob (@llist) {
$matref->{$symbol}->[$cnt] = $prob;
$cnt+=1;
}
}
close(INFILE);
return $matref;
}

# }}}

#{{{ sub dump_matrix

=head2 dump_matrix

=head3 Arguments

None apart from $self

=head3 Returns

Nothing

=head3 Comment

Prints out the matrix in a formatted form. Use this to confirm
that the matrix created is indeed the same as that read from
the matrix file.

=cut

sub dump_matrix {
my $self=shift(@_);
my $matrix=$self->{matrix};
my $length=$self->{length};
my @symbols = sort(keys %{$matrix});
foreach my $symbol (@symbols) {
my @plist=map {sprintf("%4.3f", $_) } @{$matrix->{$symbol}};
print("$symbol   ",join(" ", @plist), "\n");
}

print(STDERR "    ");
foreach my $ind (0..$length-1) {
my $coltot=0;
  foreach my $sym (@symbols) {
    unless($sym eq '>') {
    $coltot += $matrix->{$sym}->[$ind];
  }
  }
my $pcoltot=sprintf("%4.3f", $coltot);
print(STDERR "$pcoltot ");
}

my $maxscore=&max_score($self);
my $minscore=&min_score($self);

print(STDERR "\nLength: $self->{length}  Max. score:  $maxscore  Min. score: $minscore\n");
}

# }}}

# {{{ sub score

=head2 score

=head3 Arguments

=over 2

=item 1

A sequence string

=back

=head3 Returns

=over 2

=item 1

The score for the given string.

=item 2

A (possibly empty) string containing an error message.

=back

=head3 Comment

If the length of the string provided is not the same as the
length of the matrix then 0 is returned as the score along
with an error message.

=cut

sub score {
my $self=shift(@_);
my $string=uc(shift(@_));
my $matrix=$self->{matrix};
my $score=0;
my $error='';
if(length($string) != $self->{length}) {
  my $error="matrix length and string length do not match";
  return(0,$error);
}
my $pos=0;
while(my $char=uc(substr($string, $pos, 1))) {
$score += $matrix->{$char}->[$pos] * $matrix->{'>'}->[$pos];
$pos+=1;
}
return($score,'');
}

# }}}

# {{{ sub max_score

=head2 max_score

=head3 Arguments

None

=head3 Returns

The maximum possible score a sequence string can achieve
where scored using this scoring matrix.

=cut


sub max_score {
my $self=shift(@_);
my $matrix=$self->{matrix};
my $matlen=$self->{length};
my @symbols = sort(keys %{$matrix});
my $maxscore=0;

foreach my $ind (0..$matlen-1) {
  my $maxwt=0;
  foreach my $symbol (@symbols) {
    if($symbol eq '>') {next;}
    my $wt=$matrix->{$symbol}->[$ind] * $matrix->{'>'}->[$ind];
    $maxwt = $wt > $maxwt ? $wt : $maxwt;
  }
  $maxscore+=$maxwt;
}
return($maxscore);
}

# }}}

# {{{ sub min_score

=head2 min_score

=head3 Arguments

None

=head3 Returns

The minimum possible score a sequence string can achieve
where scored using this scoring matrix.

=cut


sub min_score {
my $self=shift(@_);
my $matrix=$self->{matrix};
my $matlen=$self->{length};
my @symbols = sort(keys %{$matrix});
my $minscore=0;

foreach my $ind (0..$matlen-1) {
  my $minwt=1e50;
  foreach my $symbol (@symbols) {
    if($symbol eq '>') {next;}
    my $wt=$matrix->{$symbol}->[$ind] * $matrix->{'>'}->[$ind];
    $minwt = $wt < $minwt ? $wt : $minwt;
  }
  $minscore+=$minwt;
}
return($minscore);
}

# }}}

# {{{ sub score_string

=head2 score_string()

=head3 Arguments

=over 2

=item 1

A sequence string

=back

=head3 Returns

A reference to an array containing listrefs
[$lenstr, $distance, $score]
where $lenstr is the string for which this is the score,
$distance is the distance from here to the start of the
gene and $score is the score.

=cut


sub score_string {
  my $self=shift(@_);
  my $string=uc(shift(@_));
  my $stringlen=length($string);
  my $matrix=$self->{matrix};
  my @scores=();
#print(STDERR "$string\n");
  my $pos=0;
  while (my $lenstr=substr($string, $pos, $self->{length})) {
    if (length($lenstr) < $self->{length}) {
      last;
    } else {
      my $score=0;
      my $pos1=0;
      while (my $char=uc(substr($lenstr, $pos1, 1))) {
	$score += $matrix->{$char}->[$pos1] * $matrix->{'>'}->[$pos1];
	$pos1+=1;
      }
#print(STDERR "$lenstr  $score\n");
      my $distance = $stringlen-($pos + $self->{length});
      push(@scores, [$lenstr, $distance, $score]);
    }
    $pos+=1;
  }
  return(\@scores);
}

# }}}

return 1;


__END__

