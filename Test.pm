package Sco::Test;
use 5.14.0;
use diagnostics;
use Carp;

my $number = 37;

our $num1 = 137;


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

# {{{ ### sub AUTOLOAD ###
sub AUTOLOAD {
  our $AUTOLOAD;
  my $self=shift(@_);
  my ($var) = $AUTOLOAD=~m/.*::(\w+)$/;
  if(exists $self->{$var}) {
    return($self->{$var});
  }
  else {
    carp("$var not found.\n");
    return(undef);
  }
}
# }}}

sub DESTROY {
print(STDERR "From DESTROY in Test.pm\n");
}

return(1);
