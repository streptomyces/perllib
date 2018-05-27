package Sco::Codon;
use 5.14.0;
our $AUTOLOAD;
use Bio::SeqIO;
use lib qw(/home/sco /home/sco/perllib);

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

# {{{ %ct
my %ct = (
TTT => "F",
TTC => "F",
TTA => "L",
TTG => "L",
TCT => "S",
TCC => "S",
TCA => "S",
TCG => "S",
TAT => "Y",
TAC => "Y",
TAA => "*",
TAG => "*",
TGT => "C",
TGC => "C",
TGA => "*",
TGG => "W",
CTT => "L",
CTC => "L",
CTA => "L",
CTG => "L",
CCT => "P",
CCC => "P",
CCA => "P",
CCG => "P",
CAT => "H",
CAC => "H",
CAA => "Q",
CAG => "Q",
CGT => "R",
CGC => "R",
CGA => "R",
CGG => "R",
ATT => "I",
ATC => "I",
ATA => "I",
ATG => "M",
ACT => "T",
ACC => "T",
ACA => "T",
ACG => "T",
AAT => "N",
AAC => "N",
AAA => "K",
AAG => "K",
AGT => "S",
AGC => "S",
AGA => "R",
AGG => "R",
GTT => "V",
GTC => "V",
GTA => "V",
GTG => "V",
GCT => "A",
GCC => "A",
GCA => "A",
GCG => "A",
GAT => "D",
GAC => "D",
GAA => "E",
GAG => "E",
GGT => "G",
GGC => "G",
GGA => "G",
GGG => "G"
);
# }}}

my %a2c;
_aa2codons();

# {{{ sub new {
# The returned instance will contain two hashrefs.
# $self->{c2a}
# $self->{a2c}
sub new {
  my($class, $self);
  $class=shift(@_);
  my %argv = @_;
  $self={};
  foreach my $key (%argv) {
    $self->{$key} = $argv{$key};
  }
    $self->{c2a} = \%ct;
#     _aa2codons();
    $self->{a2c} = \%a2c;
  bless($self, $class);
  return($self);
}
# }}}

sub cutable {
my $self = shift(@_);
my %args = @_;
my $infile = $args{file};
my $seqio=Bio::SeqIO->new(-file => $infile);
while(my $seqobj=$seqio->next_seq()) {
  foreach my $feature ($seqobj->all_SeqFeatures()) {
    if($feature->primary_tag() eq 'CDS') {
      my $aaobj = Sco::Genbank->feat_translate($feature);
}

}
}
}

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

# {{{ sub _aa2codons {
sub _aa2codons {
for my $codon (sort keys(%ct)) {
my $aa = $ct{$codon};
push(@{$a2c{$aa}}, $codon)
}
}
# }}}


return(1);
