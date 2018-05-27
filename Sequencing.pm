package Sco::Sequencing;
use 5.14.0;
use Bio::SeqIO;
use Bio::Seq;
our $AUTOLOAD;

# The offsets used to convert to Phred quality scores
our $illuminaQualityOffset = 64;
our $solexaQualityOffset = 33;

sub new {
  my($class, $self);
  $class=shift(@_);
  my %argv = @_;
  $self={};
  foreach my $key (keys %argv) {
    $self->{$key} = $argv{$key}
  }
  bless($self, $class);
  return($self);
}

# {{{ sub countMatchingReads hash(file, column, pattern, regex, numeric, test). 
# Returns hash(reads, matches);
# Set regex to true if pattern is a regular expression.
sub countMatchingReads {
  my $self = shift(@_);
  my %args = @_;
  my $testCnt = $args{test};
  my $col = $args{column}; # Zero indexed.
    my $pattern = $args{pattern};
  if($args{regex}) {
    my $temp = $pattern;
    $pattern = qr/$temp/;
  }
  open(my $ifh, "<", $args{file});
  my $lineCnt = 0;
  my $matchCnt = 0;
  while(<$ifh>) {
    my $line = $_;
    if($line =~ m/^@/) { next; } # Ignore SAM header lines.
      $lineCnt += 1;
      if($testCnt and ($lineCnt > $testCnt)) { last; }
    chomp($line);
    my @ll = split(/\t/, $line);
    if($args{regex}) {
      if($ll[$col] =~ m/$pattern/) {
        $matchCnt += 1;
      }
    }
    else {                        # Not a regex
      if($args{numeric}) {
        if($ll[$col] == $pattern) { $matchCnt += 1; }
      }
      else {
        if($ll[$col] eq $pattern) {
          $matchCnt += 1;
        }
      }
    }
  }
  return(reads => $lineCnt, matches => $matchCnt);
}
# }}}


# {{{ bowtieFlags (flag) returns (hashref, listref).
sub bowtieFlags {
  my $self = shift(@_);
  my $flag = shift(@_);
  my %rethash;
  my $binStr = sprintf("%0.8b", $flag);
  my @keys =  qw(paired alnPE noAln noAlnPE revStrand mateOnRev mate1 mate2);
  my @bf = split(//, $binStr);
  @bf = reverse(@bf);
  if($bf[0]) {$rethash{paired} = 1;} else {$rethash{paired} = 0;}
  if($bf[1]) {$rethash{alnPE} = 1;} else {$rethash{alnPE} = 0;}
  if($bf[2]) {$rethash{noAln} = 1;} else {$rethash{noAln} = 0;}
  if($bf[3]) {$rethash{noAlnPE} = 1;} else {$rethash{noAlnPE} = 0;}
  if($bf[4]) {$rethash{revStrand} = 1;} else {$rethash{revStrand} = 0;}
  if($bf[5]) {$rethash{mateOnRev} = 1;} else {$rethash{mateOnRev} = 0;}
  if($bf[6]) {$rethash{mate1} = 1;} else {$rethash{mate1} = 0;}
  if($bf[7]) {$rethash{mate2} = 1;} else {$rethash{mate2} = 0;}
  return(\%rethash, \@keys);
}
# }}}

# {{{ sub bowtieSeqHead
sub bowtieSeqHead {
# @SQ     SN:Sven LN:8226158
  my $self = shift(@_);
  my $hline = shift(@_);
  # print(STDERR "$hline\n");
  chomp($hline);
  $hline =~ s/^\@//;
  my %rethash;
  my ($type, $info) = split(/\t/, $hline, 2); # Limiting to 2 is important.
  if($type eq "SQ") {
    my @ilist = split(/\s+/, $info);
    foreach my $in (@ilist) {
      my ($key, $val) = split(/:/, $in);
      # print(STDERR "$key  $val\n");
      $rethash{$key} = $val;
    }
  }
  return(%rethash);
}
# }}}

# {{{ sub samSeqLens (samfilename) returns a hash.
sub samSeqLens {
  my $self = shift(@_);
  my $samfn = shift(@_);
  my %rethash;
  open(my $fh, "<", $samfn);
  while(<$fh>) {
    my $line = $_;
    chomp($line);
    unless($line =~m/^\@/) { last; }
    my @llist = split(/\t/, $line);
    if($llist[0] =~ m/SQ/) {
      shift(@llist);
      my $name; my $len;
      foreach my $sqf (@llist) {
        my @flist = split(/:/, $sqf);
        if($flist[0] eq "SN") { $name = $flist[1]; }
        if($flist[0] eq "LN") { $len = $flist[1]; }
      }
      if($name and $len) { $rethash{$name} = $len; }
    }
  }
  close($fh);
return(%rethash);
}
# }}}

# {{{ sub samLineToHash (ALineOfSamFile) returns a hash.
sub samLineToHash {
  my $self = shift(@_);
  my $samline = shift(@_); chomp($samline);
  if($samline =~ m/^\@/) { return undef; }
  my %rh;
  my @llist= split(/\t/, $samline, 12);
  $rh{readname} = $llist[0];
  $rh{flags} = $llist[1];
  $rh{refname} = $llist[2];
  $rh{alnstart} = $llist[3];
  $rh{mapqual} = $llist[4];
  $rh{cigar} = $llist[5];
  $rh{materef} = $llist[6];
  $rh{matealnstart} = $llist[7];
  $rh{fragsize} = $llist[8];
  $rh{readseq} = $llist[9];
  $rh{readlen} = length($llist[9]);
  $rh{alnend} = $llist[3] + (length($llist[9]) - 1);
  $rh{readqual} = $llist[10];
  my @optlist = split(/\t/, $llist[11]);
  $rh{optional}=[@optlist];
return(%rh);
}
# }}}

### more subs here ###

# {{{ fastqLengths (fastqfilename) returns(hash(1,2,3...)); Not implemented.
sub fastqLengths {
}
# }}}

# {{{ seqAndQuality (hash(seqobj, qtype)) returns($seqobj->seq(), [@qualities]);
sub seqAndQuality {
  my $self = shift(@_);
  my %args = @_;
my $seqobj = $args{seqobj};
my $temp = $seqobj->qual();
my @pq;
my @qtemp = @{$temp};
if($args{qtype}=~m/^i/i) {
@pq = map {$_ - 31} @qtemp;
}
else {
  @pq = @qtemp;
}
return($seqobj->seq(), [@pq]);
}
# }}}

# {{{ qualityTrim (hash(seqobj, qtype, qthresh, minlen))
sub qualityTrim {
  my $self = shift(@_);
  my %args = @_;
my $seqobj = $args{seqobj};
my $temp = $seqobj->qual();
my @pq;
my @qtemp = @{$temp};
if($args{qtype}=~m/^i/i) {
@pq = map {$_ - 31} @qtemp;
}
else {
  @pq = @qtemp;
}
my $rt = &_right_trim(\@pq, $args{qthresh});
if($rt eq "none") {
return($seqobj);
}
elsif($rt >= $args{minlen}) {
my $truncobj = $seqobj->trunc(1, $rt);
return($truncobj);
}
return(0);
}
# }}}

# {{{ internal sub _right_trim {
sub _right_trim {
  my $qref = shift(@_);
  my $qthresh = shift(@_);
  my $bpindex = 0; 
  foreach my $qu (@{$qref}) {
    $bpindex += 1;
    if($qu < $qthresh) {
      return($bpindex-1);
    }
  }
  return("none");
}
# }}}


# {{{ sub AUTOLOAD

=head2 sub AUTOLOAD()

Does its trick whenever a non-existent sub is called. Note that
if a non-existent variable is asked for then I<undef> is returned
without any warning. Caller should check the value received.

Remember to put an C<our($AUTOLOAD)> in the file scope to keep
C<use strict;> happy.

=cut

sub AUTOLOAD {
my $self=shift(@_);
my ($var) = $AUTOLOAD=~m/.*::(\w+)$/;
if(exists($self->{$var})) {
return($self->{$var});
}
else {
return(undef);
}
}
# }}}

return(1);

