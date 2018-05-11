package Sco::Common;
use 5.14.0;
use Exporter qw(import);
our @EXPORT_OK = qw(
tablist linelist tablistE linelistE tabhash tabhashE tabvals
tablistV tablistVE linelistV linelistVE tablistH linelistH
tablistER tablistVER linelistER linelistVER tabhashER tabhashVER
ymd ddmyhms
);
# use Sco::Common qw(tablist linelist tablistE linelistE tabhash tabhashE tabvals
#                    tablistV tablistVE linelistV linelistVE);

my @days = qw(Sun Mon Tue Wed Thu Fri Sat);

# {{{ subroutines tablist, linelist, tabhash and their *E versions.
# The E versions are for printing to STDERR.
# tabvals is also here.

sub tablistVER {
  if($main::verbose) {
  my @in = @_;
  print(STDERR join("\t", @in), "\r");
  }
}


sub tablistVE {
  if($main::verbose) {
  my @in = @_;
  print(STDERR join("\t", @in), "\n");
  }
}

sub tablistV {
  if($main::verbose) {
  my @in = @_;
  print(join("\t", @in), "\n");
  }
}

sub tablist {
  my @in = @_;
  print(join("\t", @in), "\n");
}

sub tablistER {
  my @in = @_;
  print(STDERR join("\t", @in), "\r");
}


sub tablistE {
  my @in = @_;
  print(STDERR join("\t", @in), "\n");
}


sub linelistV {
  if($main::verbose) {
  my @in = @_;
  print(join("\n", @in), "\n");
  }
}

sub linelistVE {
  if($main::verbose) {
  my @in = @_;
  print(STDERR join("\n", @in), "\n");
  }
}

sub linelistVER {
  if($main::verbose) {
  my @in = @_;
  print(STDERR join("\n", @in), "\r");
  }
}


sub linelist {
  my @in = @_;
  print(join("\n", @in), "\n");
}

sub linelistE {
  my @in = @_;
  print(STDERR join("\n", @in), "\n");
}

sub linelistER {
  my @in = @_;
  print(STDERR join("\n", @in), "\r");
}


sub tabhash {
  my %in = @_;
  for my $key (sort keys(%in)) {
    print(join("\t", $key, $in{$key}), "\n");
  }
}

sub tabhashE {
  my %in = @_;
  for my $key (sort keys(%in)) {
    print(STDERR join("\t", $key, $in{$key}), "\n");
  }
}

sub tabhashER {
  my %in = @_;
  for my $key (sort keys(%in)) {
    print(STDERR join("\t", $key, $in{$key}), "\r");
  }
}

sub tabvals {  # takes a hashref and a listref.
  my $hr = shift(@_);
  my $lr = shift(@_);
  my @temp;
  if(ref($lr)) {
    for my $key (@{$lr}) {
      push(@temp, $hr->{$key});
    }
  }
  else {
    for my $key (@{$hr}) {
      push(@temp, $hr->{$key});
    }
  }
  print(join("\t", @temp), "\n");
}

# }}}

# {{{ tablistH and linelistH
sub tablistH {
my $fh = shift(@_);
my @in = @_;
print($fh join("\t", @in), "\n");
}

sub linelistH {
  my $fh = shift(@_);
  my @in = @_;
  print($fh join("\n", @in), "\n");
}
# }}}

# {{{ sub ymd (year month date)
sub ymd {
my $time = shift(@_);
unless($time) { $time = time(); }
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime($time);
$year += 1900;
my $month = $mon + 1;
return($year, $month, $mday);
}
# }}}

# {{{ sub ddmyhms (weekday date month year hour min sec)
sub ddmyhms {
my $time = shift(@_);
unless($time) { $time = time(); }
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime($time);
$year += 1900;
my $xday = $days[$wday];
my $month = sprintf("%2d", $mon + 1);

my $rsec = sprintf("%2d", $sec);
my $rmin = sprintf("%2d", $min);
my $rhour = sprintf("%2d", $hour);

return($xday, $mday, $month, $year, $rhour, $rmin, $rsec);
}

# }}}

return(1);

