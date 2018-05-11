package Sco::Dolist;   # change this to the file name
use 5.14.0;
use lib qw(/home/sco/perllib);
use Sco::Common qw(
tablist linelist tablistE linelistE tabhash tabhashE tabvals
tablistV tablistVE linelistV linelistVE tablistH linelistH
tablistER tablistVER linelistER linelistVER tabhashER tabhashVER
ymd ddmyhms
);
our $AUTOLOAD;

# {{{ sub new {
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


### more subs here ###

# {{{ sub mremove.
sub mremove {
  my $self = shift(@_);
  my $serial = shift(@_);
  my $count = shift(@_);
  my @do = $self->readFromFile();
    splice(@do, $serial, $count);
  $self->writeToFile(@do);
  return;
}
# }}}

# {{{ sub remove (also used for replacement).
# This sub is also used for replacement.
sub remove {
  my $self = shift(@_);
  my $serial = shift(@_);
  my $replacement = shift(@_);
  my @do = $self->readFromFile();
  if($replacement) {
    splice(@do, $serial, 1, $replacement);
  }
  else {
    splice(@do, $serial, 1);
  }
  $self->writeToFile(@do);
  return;
}
# }}}

# {{{ sub readFromDone
sub readFromDone {
my $self = shift(@_);
my @done;
open(my $ifh, "<", $self->{donefile});
while(<$ifh>) {
  chomp;
  push(@done, $_);
}
close($ifh);
return(@done);
}
# }}}

# {{{ sub readFromFile
sub readFromFile {
my $self = shift(@_);
my @do;
open(my $ifh, "<", $self->{file});
while(<$ifh>) {
  chomp;
  push(@do, $_);
}
close($ifh);
return(@do);
}
# }}}

# {{{ sub add
sub add {
  my $self = shift(@_);
  my @args = @_;
  open(my $ifn, "<", $self->{file});
  my $serial = 0;
  my @do;
  while(<$ifn>) {
    chomp();
    push(@do, $_);
  }
  if($args[0]) {
    if(defined($args[1])) {
      splice(@do, $args[1], 0, $args[0]);
    }
    else { push(@do, $args[0]); }
    $self->writeToFile(@do);
  }
  else { linelistE("No item provided"); }
  return;
}
# }}}

# {{{ sub move
sub move {
  my $self = shift(@_);
  my $opos = shift(@_);
  my $npos = shift(@_);
#  tablistE($opos, $npos);
  my @do;
  open(my $ifn, "<", $self->{file});
  while(<$ifn>) {
    chomp($_);
    push(@do, $_);
  }
  close($ifn);

  my $item = splice(@do, $opos, 1);
  splice(@do, $npos, 0, $item);

  open(my $ofn, ">", $self->{file});
  my $oldselect = select($ofn);
  linelist(@do);
  close($ofn);
  select($oldselect);
  return;
}
# }}}

# {{{ sub show
sub show {
  my $self = shift(@_);
  open(my $ifh, "<", $self->{file});
  my $serial = 0;
  print("\n");
  # linelist("File: $self->{file}");
  print("\n");
  while(<$ifh>) {
    chomp();
    tablist(sprintf("%02d", $serial), $_);
    $serial += 1;
  }
  print("\n\n");
  close($ifh);
}
# }}}

# {{{ sub done
sub done {
  my $self = shift(@_);
  my $serial = shift(@_);
  my $tcount = shift(@_);
  my $count = $tcount > 1 ? $tcount : 1;

  my @do = $self->readFromFile();
  my @dout = splice(@do, $serial, $count);
  $self->writeToDone(@dout);
  $self->writeToFile(@do);
  return;
}
# }}}


# {{{ sub showdone
sub showdone {
  my $self = shift(@_);
  open(my $ifh, "<", $self->{donefile});
  my $serial = 0;
  print("\n");
  # linelist("File: $self->{file}");
  print("\n");
  while(<$ifh>) {
    chomp();
    linelist($_);
    $serial += 1;
  }
  print("\n\n");
  close($ifh);
}
# }}}


# {{{ sub file
sub file {
  my $self = shift(@_);
  linelist($self->{file});
}
# }}}

# {{{ sub hostname
sub hostname {
  my $self = shift(@_);
  my $hostname = qx(hostname);
  chomp($hostname);
  linelist($hostname);
}
# }}}

# {{{ sub writeToDone
sub writeToDone {
  my $self = shift(@_);
  my @args = @_;
open(my $fh, ">>", $self->{donefile});
my $before = select($fh);
my @ds = split(/\s+/, localtime());
my $year = pop(@ds);
unshift(@ds, $year);

for my $arg (@args) {
tablist(join(" ", @ds), $arg);
}

select($before);
  close($fh);
}
# }}}

# {{{ sub writeToFile
sub writeToFile {
  my $self = shift(@_);
  my @args = @_;
open(my $fh, ">", $self->{file});
my $before = select($fh);
linelist(@args);
select($before);
  close($fh);
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

