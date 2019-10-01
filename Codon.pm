package Sco::Codon;
use 5.14.0;
our $AUTOLOAD;
use Bio::SeqIO;
use lib qw(/home/sco /home/sco/perllib);
use Sco::Common qw(tablist linelist tablistE linelistE tabhash tabhashE tabvals
    tablistV tablistVE linelistV linelistVE tablistH linelistH
    tablistER tablistVER linelistER linelistVER tabhashER tabhashVER csvsplit);

# # {{{ ### class data and methods ###
# {
# 
# my $classdata1='classdata1';
# my $classdata2='classdata2';
# 
# # {{{ ### A get or set method for class data ###
# sub cd1 {
#   my $self=shift(@_);
#   my $temp=shift(@_);
#   if (defined($temp)) {
#     $classdata1=$temp;
#     return(1);
#   } else {
#     return($classdata1);
#   }
# }
# # }}} ##########################################
# 
# }
# # }}} ### end of class data and methods ###

# {{{ %ct. Hash. Translation table.
my %ct = (
ttt => "F",
ttc => "F",
tta => "L",
ttg => "L",
tct => "S",
tcc => "S",
tca => "S",
tcg => "S",
tat => "Y",
tac => "Y",
taa => "*",
tag => "*",
tgt => "C",
tgc => "C",
tga => "*",
tgg => "W",
ctt => "L",
ctc => "L",
cta => "L",
ctg => "L",
cct => "P",
ccc => "P",
cca => "P",
ccg => "P",
cat => "H",
cac => "H",
caa => "Q",
cag => "Q",
cgt => "R",
cgc => "R",
cga => "R",
cgg => "R",
att => "I",
atc => "I",
ata => "I",
atg => "M",
act => "T",
acc => "T",
aca => "T",
acg => "T",
aat => "N",
aac => "N",
aaa => "K",
aag => "K",
agt => "S",
agc => "S",
aga => "R",
agg => "R",
gtt => "V",
gtc => "V",
gta => "V",
gtg => "V",
gct => "A",
gcc => "A",
gca => "A",
gcg => "A",
gat => "D",
gac => "D",
gaa => "E",
gag => "E",
ggt => "G",
ggc => "G",
gga => "G",
ggg => "G"
);
# }}}

=head1 Name

Codon.pm

=head2 Examples

 use Sco::Codon
 

=cut

my %a2c; # keys are aminoacids. values are lists of codons. Populated in the call to new.
my %naa; # number of codons encoding each amino acid. Populated in the call to new.

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
  %a2c = _aa2codons();
  %naa = _naa();
  $self->{a2c} = \%a2c;
  return($self);
}
# }}}

# {{{ sub codon_counts
sub codon_counts {
my $self = shift(@_);
my %args = @_;
my $in = $args{seq};
my $seq; # nt sequence string.
if(ref($in)) {
  $seq = $in->seq();
}
else {
  $seq = $in;
}

$seq = lc($seq);

my $pos = 0;
my $codon;
my %rethash = _codon_zero_hash();
while($codon = substr($seq, $pos, 3)) {
$rethash{$codon} += 1;
$pos += 3;
}
return(%rethash);
}
# }}}

# {{{ sub rscu
sub rscu {
my $self = shift(@_);
my %args = @_;
my $ccr = $args{codon_counts};
my %cc = %{$ccr};
my %rethash;
my %rscu;

for my $aa (sort keys(%naa)) {
my @codons = @{$a2c{$aa}};
my $ccnt = $naa{$aa};

my $aacnt;
for my $codon (@codons) {
$aacnt += $cc{$codon};
}

my $unbias = $aacnt / $ccnt;

for my $codon (@codons) {
$rscu{$codon} = sprintf("%0.3f", $cc{$codon} / $unbias);
}
}

return(%rscu);
}
# }}}

# {{{ sub write_rscu
sub write_rscu {
my $self = shift(@_);
my %args = @_;
my $rr = $args{rscu};
my %rscu = %{$rr};
my @retlist;
for my $aa (sort keys(%a2c)) {
my $temp = $a2c{$aa};
my @codons = @{$temp};
@codons = sort(@codons);
for my $codon (@codons) {
push(@retlist, [$aa, $codon, $rscu{$codon}]);
}
}
return(@retlist);
}
# }}}

# {{{ sub _aa2codons
# No arguments. Reads %ct. Returns
sub _aa2codons {
  my %a2c;
  for my $codon (sort keys(%ct)) {
    my $aa = $ct{$codon};
    push(@{$a2c{$aa}}, $codon)
  }
  return(%a2c);
}
# }}}

# {{{ sub _codon_zero_hash
sub _codon_zero_hash {
my %zero;
for my $codon (keys %ct) {
$zero{$codon} = 0;
}
return(%zero);
}
# }}}

# {{{ sub _naa. number of codons for given AA.
sub _naa {
my %naa;
for my $aa (values(%ct)) {
$naa{$aa} += 1;
}
return(%naa);
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
    return($AUTOLOAD);
  }
}
# }}}


return(1);

__END__

# Bio::Tools::CodonTable,
# Bio::WebAgent,
# Bio::CodonUsage::IO,
# Bio::CodonUsage::Table,
# Bio::DB::CUTG
