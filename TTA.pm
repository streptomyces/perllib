package Sco::TTA;
use common::sense;
our $AUTOLOAD;
use lib qw(/home/sco /home/sco/perllib);
use Bio::SeqIO;
use Bio::Seq;

# {{{ new
sub new {
  my($class, $self);
  $class=shift(@_);
  my %argv = @_;
  $self={};
  foreach my $key (%argv) {
    $self->{$key} = $argv{$key}
  }
  bless($self, $class);
  return($self);
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

### more subs here ###

# {{{ inframeTTA (bioseq object or string of nt sequence) returns(@TTA positions)
sub inframeTTA {
my $self=shift(@_);
my $inseq=shift(@_);
my @retlist;
my $seq;
if(ref($inseq)) { $seq = $inseq->seq(); }
else { $seq = $inseq; }
#print(STDERR "$seq\n");
my $pos = 0;
while(my $triplet = substr($seq, $pos, 3)) {
  if($triplet=~m/[tu][tu]a/i) {
    push(@retlist, $pos + 1);
  }
$pos+=3;
}
return(@retlist); # change this
}
# }}}




# {{{ somesub (args) returns(retvals)
sub somesub {
my $self=shift(@_);
my %args=@_;
return('something'); # change this
}
# }}}








return(1);
