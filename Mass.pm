package Sco::Mass;   # change this to the file name
use 5.14.0;
use utf8;
use Carp;
use File::Basename;
our $AUTOLOAD;

# {{{ Weights.
my $water = 18.0106;
my $twohydro = 2.0157;
my $hydroion = 1.0073;
my $twohydroion = 2.0146;
my $sodium = 22.9892;
# 2.0157
my %modification = (
LossOf2Hydrogen => -($twohydro),
LossOfWater => -($water)
);
# }}}

# {{{ %lut LUT of aa weights.
my %lut = (
A => 71.037114,
R => 156.101111,
N => 114.042927,
D => 115.026943,
C => 103.009185,
E => 129.042593,
Q => 128.058578,
G => 57.021464,
H => 137.058912,
I => 113.084064,
L => 113.084064,
K => 128.094963,
M => 131.040485,
F => 147.068414,
P => 97.052764,
S => 87.032028,
T => 101.047679,
U => 150.95363,
W => 186.079313,
Y => 163.06332,
V => 99.068414
);
# }}}

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


# {{{ Internal sub _adduct. (mass, adductString);
# This is the mass after modification.
sub _adduct {
my $mass = shift(@_);
my $adduct = shift(@_);

if($adduct eq "None") {
return($mass);
}
elsif($adduct eq "+H") {
return($mass + $hydroion);
}
elsif($adduct eq "+2H") {
return(($mass + $twohydroion)/2);
}
elsif($adduct eq "+Na") {
return($mass + $sodium);
}
elsif($adduct eq "+Na+H") {
return(($mass + $sodium + $hydroion)/2);
}
}
# }}}

# {{{ sub adductedMasses3.
# peptideSequence
# lr of modifications
# lr of adducts
# lr of increases
# lr of decreases
# Returns a hash.
sub adductedMasses3 {
  my $self = shift(@_);
  my $temp = shift(@_);  # Sequence to be upcased later.
  my $modr = shift(@_);  # list ref to modifications.
  my $addr = shift(@_);  # list ref to adducts.
  my $incr = shift(@_);  # list ref to increases.
  my $decr = shift(@_);  # list ref to decreases.
  my $seq = uc($temp);
  my @mods;
  for my $temp (@{$modr}) {
    push(@mods, $modification{$temp});
  }

  for my $temp (@{$incr}) {
  if($temp != 0) { push(@mods, $temp); }
  }
  for my $temp (@{$decr}) {
  if($temp != 0) { push(@mods, -($temp)); }
  }


  print(STDERR join("\t", @mods), "\n");


  my @adducts = @{$addr};
  my $lumass = 0;
  for my $aa (split("", $seq)) {
    $lumass += $lut{$aa};
  }
  my $mass = $lumass + $water;
  my %mh;
  $mh{None} = $mass;

  my @mcr = _combinations(@mods); # mod combinations listrefs

  for my $mcr (@mcr) {
    my @tmod = @{$mcr};
    unless(@tmod) { next; }
    my $cycmass = $mass;
    my @prekey = ();
  for my $mod (@tmod) {
      $cycmass = $cycmass + $mod; 
      push(@prekey, $mod);
  }
  my $key = join(",", @prekey);
  $mh{$key} = $cycmass;
  }

  my %rh;
  for my $key (keys(%mh)) {
    for my $adduct (@adducts) {
    $rh{$key.":".$adduct} = _adduct($mh{$key}, $adduct);;
    }
  }
  # else { %rh = %mh; }
  return(%rh);
}
# }}}

# {{{ sub adductedMasses2. (peptideSequence, lr of modifications, 
# lr of adducts, lr of increases, lr of decreases).
# Returns a hash.
sub adductedMasses2 {
  my $self = shift(@_);
  my $temp = shift(@_);  # Sequence to be upcased later.
  my $modr = shift(@_);  # list ref to modifications.
  my $addr = shift(@_);  # list ref to adducts.
  my $incr = shift(@_);  # list ref to increases.
  my $decr = shift(@_);  # list ref to decreases.
  my $seq = uc($temp);
  my @mods = @{$modr};
  my @adducts = @{$addr};
  my $lumass = 0;
  for my $aa (split("", $seq)) {
    $lumass += $lut{$aa};
  }
  my $mass = $lumass + $water;
  my %mh;
  $mh{None} = $mass;
  for my $mod (@mods) {
      $mh{$mod} = $mass + $modification{$mod}; 
  }
  for my $inc (@{$incr}) {
    $mh{"+", $inc} = $mass + $inc;
  }
  for my $dec (@{$decr}) {
    $mh{"-", $dec} = $mass - $dec;
  }

  my %rh;
  for my $key (keys(%mh)) {
    for my $adduct (@adducts) {
    $rh{$key.":".$adduct} = _adduct($mh{$key}, $adduct);;
    }
  }
  # else { %rh = %mh; }
  return(%rh);
}
# }}}

# {{{ sub adductedMasses. (peptideSequence, lr of modifications, lr of adducts).
# Returns a hash.
sub adductedMasses {
  my $self = shift(@_);
  my $temp = shift(@_);  # Sequence to be upcased later.
  my $modr = shift(@_);  # list ref to modifications.
  my $addr = shift(@_);  # list ref to adducts.
  my $seq = uc($temp);
  my @mods = @{$modr};
  my @adducts = @{$addr};
  my $lumass = 0;
  for my $aa (split("", $seq)) {
    $lumass += $lut{$aa};
  }
  my $mass = $lumass + $water;
  my %mh;
  $mh{None} = $mass;
  for my $mod (@mods) {
    if( grep {$_ eq $mod} keys(%modification) ) {
      $mh{$mod} = $mass + $modification{$mod}; 
    }
    elsif ($mod !~ m/[^+\-\d\.]|[\d\.][+\-]|\..*\.|[+\-].*[+\-]/ and $mod != 0) {
      $mh{adjust} = $mass + $mod;
    }
  }

  my %rh;

  for my $key (keys(%mh)) {
    for my $adduct (@adducts) {
    $rh{$key.":".$adduct} = _adduct($mh{$key}, $adduct);;
    }
  }
  return(%rh);
}
# }}}

# {{{ sub masses. (peptideSequence, listofmodifications). returns a hash.
sub masses {
  my $self = shift(@_);
  my $temp = shift(@_);
  my @mods = @_;
  my $seq = uc($temp);
  my $lumass = 0;
  for my $aa (split("", $seq)) {
    $lumass += $lut{$aa};
  }
  my $mass = $lumass + $water;
  my %mh;
  $mh{None} = $mass;
  for my $mod (@mods) {
    if( grep {$_ eq $mod} keys(%modification) ) {
      $mh{$mod} = $mass + $modification{$mod}; 
    }
    elsif ($mod !~ m/[^+\-\d\.]|[\d\.][+\-]|\..*\.|[+\-].*[+\-]/ and $mod != 0) {
      $mh{adjust} = $mass + $mod;
    }
  }

  my %rh;

  for my $key (keys(%mh)) {
    $rh{$key.":plus1Hion"} = $mh{$key} + $hydroion;
    $rh{$key.":plus2Hion"} = ($mh{$key} + $twohydroion)/2;
    $rh{$key.":plusNa"} = $mh{$key} + $sodium;
    $rh{$key.":plusNaHion"} = ($mh{$key} + $hydroion + $sodium) / 2;
  }
  return(%rh);
}
# }}}

# {{{ ### sub charge ###
sub somesub {
my $self = shift(@_);
my @args = @_;
# code here to do things #

return('something'); # change this
}
# }}}

# {{{ Internal sub _combinations (list) (list of listrefs);
sub _combinations {
  return [] unless @_;
  my $first = shift;
  my @rest = _combinations(@_);
  return @rest, map { [$first, @$_] } @rest;
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


__END__

1. No need for “none” in the “Form" column in the results - the charged options
are sufficient.

2. Can the modifications be additive?  i.e. if I select loss of water, loss of
2H and addition of a given mass, the results list every possible combination of
these modifications (from none through to every modification made).

3. The “loss of hydrogen” checkbox is not necessary (but then isn’t a problem
if it isn’t checked so not essential to change).



Rules for searching for ribosomal peptide masses.

Calculate mass of any contiguous sequence of amino acids from a given sequence
within a peptide.
Test example: MDIRDEAGAEQATAGAGEPLPDLLELDLAQLGTVEHPVLREVLGELRARAAEPSEMLWGFDNSF

Use following table for masses:

A 71.037114
R 156.101111
N 114.042927
D 115.026943
C 103.009185
E 129.042593
Q 128.058578
G 57.021464
H 137.058912
I 113.084064
L 113.084064
K 128.094963
M 131.040485
F 147.068414
P 97.052764
S 87.032028
T 101.047679
U 150.95363
W 186.079313
Y 163.06332
V 99.068414

Mass of peptide = Sum of masses above + 18.0106
For each peptide, also calculate masses for peptides with the following losses:
H2 = - 2.0157
2 H2 = - 4.0313

Therefore, this gives 3 peptides per any given sequence:
1. Full
2. - H2
3. - 2H2

Constraint search for any peptides between 3 and 30 amino acids long.

Search for the following forms of each peptide:
[Peptide+H]+ = Peptide mass + 1.0073
[Peptide + 2H]2+ = (peptide mass + 2.0146)/2
[Peptide + Na]+ = Peptide mass + 22.9892
[Peptide + Na + H]2+ = (Peptide mass + 1.0073 + 22.9892)/2

Therefore, using the peptide sequence SEMLWGFDNS, the following 
masses should be obtained (if I've done my maths correctly):

591.2320795
592.2398795
593.2477295
602.2230295
603.2308295
604.2386795
1181.456859
1183.472459
1185.488159
1203.438759
1205.454359
1207.470059

An output list ordered by mass would be great (csv format or just different
lines in a text file?). Accuracy of 4 decimal places will be ideal.

