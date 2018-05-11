package Sco::Conf;
use 5.14.0;
use lib qw(/home/sco/perllib);
use Sco::Common qw(tablist linelist tablistE linelistE tabhash tabhashE tabvals);
our $AUTOLOAD;

=head1 Name

Sco::Conf

=head2 Description

This module allows you to make an object to hold configuration
variables. The initial configuration may be read from a file which
should be a two column file where columns are space separated.
The first column should not have any spaces in it.

Lines beginning with '#' are ignored.

You can have comments at the end of lines and they are ignored
as well.

=head2 Configuration file

 # This is a comment and the whole line is ignored.
 key    The value can have spaces in it.    # This comment is ignored.
 key2 Another value.
 key3 Value3

=cut


sub new {
  my($class, $self);
  $class=shift(@_);
  my %argv = @_;
  $self={};
  foreach my $key (keys %argv) {
    $self->{$key} = $argv{$key}
  }
  if(exists($argv{file})) {
    open(CONF, "<", $argv{file}) or croak("Failed to open $argv{file}.");
    while(my $line = readline(CONF)) {
      chomp($line);
      if($line=~m/^\s*\#/ or $line=~m/^\s*$/) {next;}
      $line =~ s/\s*#.*$//;
      my @ll=split(' ', $line, 2);
      $self->{$ll[0]} = $ll[1];
    }
    close(CONF);
  }
  bless($self, $class);
  return($self);
}

sub dump {
my $self = shift(@_);
my @keys = keys(%{$self});
for my $key (sort @keys) {
tablist("key: $key", "value: $self->{$key}");
}
}



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


