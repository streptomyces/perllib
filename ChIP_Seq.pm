package Sco::ChIP_Seq;   # change this to the file name
use 5.14.0;
our $AUTOLOAD;
use lib qw(/home/sco/perllib);

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

sub DESTROY {
my $self = shift(@_);
#$handle->disconnect();
return(1);
}


### more subs here ###

# {{{ featHash (featFile) returns(%annos)
sub featHash {
my $self=shift(@_);
my %args=@_;
my %annos;
open(my $feat, "<", $args{featFile}) or croak("Could not open $args{featFile}");
my $lineCnt = 0;
while(<$feat>) {
my $line=$_;
chomp($line);
my @llist=split(/\t/, $line);
my $lt = $llist[1];
my $anno = $llist[5];
my $geneLen = ($llist[3] - $llist[2]) + 1;
$annos{$lt} = {anno => $anno, len => $geneLen};
$lineCnt += 1;
}
close($feat);
print(STDERR "$lineCnt features read\n");
return(%annos);
}
# }}}

# {{{ sub bed2featList (bedfile)
sub bed2featList {
my $self = shift(@_);
my %args = @_;
my @retlist;
open(BED, "<", $args{bedfile});
my $trackLine = readline(BED);
while(<BED>) {
my $line = $_;
chomp($line);
my @llist = split(/\t/, $line);
my $strand;
if($llist[5] eq '+') { $strand = 1 }
elsif($llist[5] eq '-') { $strand = -1 }

push(@retlist, [ $llist[3], $llist[1], $llist[2], $strand ]);
}
close(BED);
return(@retlist);
}


# {{{ somesub (args) returns(retvals)
sub somesub {
my $self=shift(@_);
my %args=@_;
return('something'); # change this
}
# }}}






return(1);
