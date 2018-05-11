package Sco::Fastq;
use 5.14.0;
use Bio::SeqIO;
use Bio::Seq;
our $AUTOLOAD;
use lib qw(/home/sco /home/sco/perllib);
use parent qw(Sco::Global);

sub new {
        my($class, $self);
        $class=shift(@_);
        $self={};
        bless($self, $class);
        return($self);
}

# {{{ sub gcFracn hash(file, n) returns hash(file, n, meanGC);
sub gcFracn {
my $self = shift(@_);
my %args = @_;
my $seqio = Bio::SeqIO->new(-file => $args{file});
my $lCnt = 0;
my $cumulative = 0;
while(my $seqobj = $seqio->next_seq()) {
my $seq = $seqobj->seq();
# print(STDERR "$seq\n"); # testing only.
my (undef, undef, $fraction) = $self->gc_frac($seq);
# print(STDERR "Fraction: $fraction\n"); # testing only.
$cumulative += $fraction;
$lCnt += 1;
if($lCnt >= $args{n}) {
  last;
}
}
my $mean = $cumulative / $lCnt;
return(meanGC => $mean, n => $args{n}, file => $args{file});
}
# }}}

# {{{ sub numReads
sub numReads {
my $self = shift(@_);
my %args = @_;
open(my $ifh, "<", $args{file});
my $lineCnt = 0;
while(<$ifh>) {
    $lineCnt += 1;
}
close($ifh);
my $numReads = $lineCnt / 4;
return($numReads);
}
# }}}

# {{{ sub seqContent hash(file) returns an integer;
sub seqContent {
my $self = shift(@_);
my %args = @_;
my $seqio = Bio::SeqIO->new(-file => $args{file});
my $cumlen = 0;
while(my $seqobj = $seqio->next_seq()) {
$cumlen += $seqobj->length();
}
return($cumlen);
}
# }}}


# {{{ ### sub AUTOLOAD ###
sub AUTOLOAD {
  my $self=shift(@_);
  my ($var) = $AUTOLOAD=~m/.*::(\w+)$/;
  if(exists $self->{$var}) {
    return($self->{$var});
  }
  else {
    return(undef);
  }
}
# }}}


return(1);


