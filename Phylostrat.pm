package Sco::Phylostrat;   # change this to the file name
use lib qw(/home/sco /home/sco/perllib);
use common::sense;
use parent qw(Sco::Orthologs);
our $AUTOLOAD;

my $phylostratHome = qw(/home/sco/customers/phylostratigraphy);
my $handle = Sco::Orthologs->get_handle();
#my $dbname = "bigortho";
#my $dbhost = "jiilin9.jic.bbsrc.ac.uk";
#my $handle=DBI->connect("DBI:Pg:dbname=$dbname;host=$dbhost", 'sco', 'tsinH4x');




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


### more subs here ###

sub taxId2RankName {
my $self = shift(@_);
my %argv = @_;
my $qstr = qq/select name from taxnames where taxid = $argv{taxid}/;
my ($name) = $handle->selectrow_array($qstr);
my $qstr = qq/select rank from taxnodes where taxid = $argv{taxid}/;
my ($rank) = $handle->selectrow_array($qstr);
return($rank, $name);
}


# sub dNdS hash(seq1, seq2);
sub dNdS {
my $self = shift(@_);
my %argv = @_;
my $ntA = $argv{seq1};

}


# {{{ ### sub charge ###
sub somesub {
my($self, @args);
my $self=shift(@_);
@args=@_;
# code here to do things #

return('something'); # change this
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
