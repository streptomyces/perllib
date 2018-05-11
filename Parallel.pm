package Sco::Parallel;
use IO::Socket::INET;
use Bio::SeqIO;
use Getopt::Long;
use common::sense;
use Carp;
use File::Temp qw(tempfile tempdir);
my $tempdir = qw(/home/sco/volatile);
my $template="parallelXXXXX";
our($AUTOLOAD);

# {{{ new
sub new {
  my($class, $self);
  $class=shift(@_);
  my %argv = @_;
  $self={};
  foreach my $key (%argv) {
    $self->{$key} = $argv{$key}
  }
  open(CONF, "<$argv{conf}");
  while(<CONF>) {
    my $line = $_;
    chomp($line);
    if($line=~m/^\s*\#/ or $line=~m/^\s*$/) {next;}
    my @ll = split(/\t/, $line);
    $self->{$ll[0]} = $ll[1];
  }
  bless($self, $class);
  return($self);
}
# }}}

### The Query Server ###

# {{{ sub qServerStart
sub qServerStart {
  my $self = shift(@_);
  my $subref = shift(@_);
  my $infile = $self->input_file();
  if($self->{qstat} == 1) {
    print("Query server already running\n");
    return(0);
  }
#  my $seqio=Bio::SeqIO->new(-file => $infile);
  $self->{qstat} = 1;
my $retval = 1;
my $sock = IO::Socket::INET->new(Listen => 8,
LocalAddr => $self->query_host(),
LocalPort => $self->query_port(),
Proto     => 'tcp');
unless($sock) {carp "Socket Error"; return(0);}
my $client;

while($client=$sock->accept()) {
  my $seqout=Bio::SeqIO->new(-fh => $client, -format => 'fasta');
  while(<$client>) {
    my $cline=$_; chomp($cline);
    if($cline eq "close file and exit") {
      close($client);
      close($sock);
      exit;
    }
    elsif($cline eq 'next') {
      my $seqobj = $subref->();
      if($seqobj) {
      $seqout->write_seq($seqobj);
      }
      else { print($client undef); }
    }
    elsif($cline eq 'ping') {
      print($client "\npong\n");
    }
    print($client "\nEOT\n");
  }
}
return($retval);
}
# }}}


### The receiver ###

# {{{ sub receiverStart
sub receiverStart {
  my $self = shift(@_);
  my $subref = shift(@_);

my $sock = IO::Socket::INET->new(Listen => 8,
LocalAddr => $self->receiver_host(),
LocalPort => $self->receiver_port(),
Proto     => 'tcp');

unless($sock) {die "Socket Error";}
while(my $client=$sock->accept()) {
while(<$client>) {
  my $cline = $_; chomp($cline);
  if($cline eq 'close file and exit') {
    close($client);
    close($sock);
    exit;
  }
  else {
    $subref->($cline);
  }
}
}
}
# }}}

# {{{ sub _rstart_again
sub _rstart_again {
  my $fh = shift(@_);
  my $outfile = shift(@_);
  close($fh);
  open($fh, ">$outfile");
  return($fh);
}
# }}}


# {{{  sub AUTOLOAD 
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

### Functions for the worker ###

# {{{ sub get_next
sub get_next {
  my $self = shift(@_);
my $remote = IO::Socket::INET->new(
PeerAddr => $self->query_host(),
PeerPort => $self->query_port(),
Proto    => 'tcp');
print($remote "next\n");
my @rlines=();
while(<$remote>) {
  chomp();
if($_=~m/^EOT$/) {last;}
push(@rlines, $_);
}
close($remote);
return(join("\n", @rlines));
}
# }}}

# {{{ sub send_result
sub send_result {
  my $self = shift(@_);
my $result=shift(@_);
my $remote = IO::Socket::INET->new(
PeerAddr => $self->receiver_host(),
PeerPort => $self->receiver_port(),
Proto    => 'tcp');
print($remote "$result");
close($remote);
}
# }}}







return(1);


__END__


