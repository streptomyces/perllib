package Sco::Alignment;
use common::sense;
our $AUTOLOAD;

# {{{ ### class data and methods ###
{

my $classdata1='classdata1';
my $classdata2='classdata2';

# {{{ ### A get or set method for class data ###
sub cd1 {
  my $self=shift(@_);
  my $temp=shift(@_);
  if (defined($temp)) {
    $classdata1=$temp;
    return(1);
  } else {
    return($classdata1);
  }
}
# }}} ##########################################

}
# }}} ### end of class data and methods ###

my $gapChar = '-';
my @aminoAcids = qw(A C D E F G H I K L M N P Q R S T V W Y);
push(@aminoAcids, $gapChar);

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

# {{{ clustal2hash (hash(file)) returns a hash
sub clustal2hash {
  my $self = shift(@_);
  my %args = @_;


open(ALN, "<$args{file}");

my %seqs; # Holds the aligned sequences
my @names; # holds all the sequence identifiers
my $clustHead = readline(ALN);
chomp($clustHead);
$clustHead = "CLUSTAL 2.0.9 multiple sequence alignment";

while(<ALN>) {
my $line=$_;
chomp($line);
if($line=~m/^\s*\#/ or $line=~m/^\s*$/ or $line=~m/^\s/) {next;}

my @llist=split(/\s+/, $line);   ### change the split char here.
unless(grep {$_ eq $llist[0]} @names) { push(@names, $llist[0]); } # names populated here
$seqs{$llist[0]} .= $llist[1];
}

close(ALN);

my $numSequences = scalar(@names);
my $alnLength;
return(%seqs);
}
# }}}

# {{{ removeGapCols (hash(alignment, gapthresh)) returns(alignment) original alignment is unaltered.
# alignment is a hash
sub removeGapCols {
my $self = shift(@_);
my %args = @_;
my %seqs = %{$args{alignment}};
my @names = keys %seqs;
my $gapThresh = $args{gapThresh};
my $numSequences = scalar(@names);

my @gapcols; # holds column numbers which have gaps in them.
my %gapsInCol; # number of gaps in different column numbers
foreach my $temp (0..length($seqs{$names[0]})) {
$gapsInCol{$temp} = 0;
}

my $pcnt = 0;
foreach my $name (@names) {
  $pcnt += 1;
#print(STDERR "$pcnt\t$name\n");
my $pos = -1;
while(($pos = index($seqs{$name}, '-', $pos)) > -1) {
my $temp = $pos - 1;
$gapsInCol{$pos} += 1;
unless(grep {$_ == $pos} @gapcols) { push(@gapcols, $pos); }
$pos += 1;
}
}

foreach my $key (sort {$a <=> $b} keys(%gapsInCol)) {
#print("$key\t$gapsInCol{$key}\t", $gapsInCol{$key}/$numSequences,"\n");
if($gapsInCol{$key}/$numSequences >= $gapThresh ) {
foreach my $name (@names) {
#  print(STDERR "$key\t$name\n");
substr($seqs{$name}, $key, 1, '!');
}
}
}

foreach my $name (@names) {
$seqs{$name}=~s/\!//g;
}
return(%seqs);
}
# }}}


# {{{ alignment2clustal (hash(alignment, outfile)) returns nothing. If outfile is given writes to
# outfile otherwise writes to STDOUT.
sub alignment2clustal {
my $self = shift(@_);
my %args = @_;
my %seqs = %{$args{alignment}};
my @names = keys %seqs;
my $numSequences = scalar(@names);

# below is for printing out the multiple sequence alignment in clustal format.
my %chunks;
my $num_chunks = 0;
foreach my $name (@names) {
  $num_chunks = 0;
my $pos = 0;
while(my $chunk = substr($seqs{$name}, $pos, 60) ) {
if($chunk=~m/\S+/) {
push(@{$chunks{$name}}, $chunk);
$pos += 60;
$num_chunks += 1;
}
}
}

if($args{outfile}) {
  open(OUT, ">$args{outfile}");
  select OUT;
}

print("CLUSTAL 2.0.9 multiple sequence alignment\n\n\n");
foreach my $cn (0..$num_chunks-1) {
foreach my $name (@names) {
print($name, " " x (16 - length($name)), $chunks{$name}->[$cn], "\n");
}
print("\n");
}
if($args{outfile}) {
  close(OUT);
}
select(STDOUT);
}
# }}}

# {{{ jaccard hash(alignment, c1, c2) returns(hash) jaccard indices for all aa pairs
sub jaccard {
  my $self = shift(@_);
  my %args = @_;
  my %seqs = %{$args{alignment}};
  my $c1 = $args{c1};
  my $c2 = $args{c2};
  my %jis;
  my (@aa1, @aa2);
  foreach my $name (keys %seqs) {
    unless(grep({$_ eq substr($seqs{$name}, $c1, 1)} @aa1)) {
      push(@aa1, substr($seqs{$name}, $c1, 1));
    }
    unless(grep({$_ eq substr($seqs{$name}, $c2, 1)} @aa2)) {
      push(@aa2, substr($seqs{$name}, $c2, 1));
    }
  }

my @pairs;
foreach my $a1 (@aa1) {
foreach my $a2 (@aa2) {
push(@pairs, [$a1, $a2]);
}
}

foreach my $pair (@pairs) {
  my (@names1, @names2);
  foreach my $name (keys %seqs) {
  if(substr($seqs{$name}, $c1, 1) eq $pair->[0]) {
    push(@names1, $name);
  }
  }
  foreach my $name (keys %seqs) {
  if(substr($seqs{$name}, $c2, 1) eq $pair->[1]) {
    push(@names2, $name);
  }
  }

my @intsect = _intersection(\@names1, \@names2);
my @union = _union(\@names1, \@names2);
my $n_intsect = scalar(@intsect);
my $n_union = scalar(@union);
my $ji = scalar(@intsect) / scalar(@union);
my $pairkey = $pair->[0] . $pair->[1];
$jis{$pairkey} = [$n_intsect, $n_union, $ji];
}

return(%jis);

}
# }}}

# {{{ internal subs _union and _intersection
sub _intersection {
my $lr1 = shift(@_);
my $lr2 = shift(@_);
my @retlist;
foreach my $el1 (@{$lr1}) {
foreach my $el2 (@{$lr2}) {
if($el1 eq $el2) {
unless(grep {$_ eq $el2} @retlist) { push(@retlist, $el2); }
}
}
}
return(@retlist);
}

sub _union {
my $lr1 = shift(@_);
my $lr2 = shift(@_);
my @retlist;
foreach my $el1 (@{$lr1}) {
unless(grep {$_ eq $el1} @retlist) { push(@retlist, $el1); }
}
foreach my $el2 (@{$lr2}) {
unless(grep {$_ eq $el2} @retlist) { push(@retlist, $el2); }
}
return(@retlist);
}
# }}}

# {{{ sectionsOfAlignment hash(alignment, positions) returns an alignment
sub sectionsOfAlignment {
my $self = shift(@_);
my %args = @_;
my %seqs = %{$args{alignment}};
my @names = keys %seqs;
my $numSequences = scalar(@names);
my @positions = @{$args{positions}};

print(join(" ", @positions), "\n");


my %retseqs;

while(@positions) {
my @pp;
$pp[0] = shift(@positions);
$pp[1] = shift(@positions);
foreach my $name (@names) {
$retseqs{$name} .= substr($seqs{$name}, $pp[0], $pp[1] - $pp[0] + 1);
}

}
return(%retseqs);
}
# }}}

return(1);
