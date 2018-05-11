package Sco::Bed;
use 5.14.0;
use Bio::SeqIO;
use Bio::Seq;
use Carp;
use DBI;
use File::Basename;
our $AUTOLOAD;
use lib qw(/home/sco /home/sco/perllib);
use parent qw(Sco::Fasta Sco::Global);
# use Scoglobal;
use File::Temp qw(tempfile tempdir);
my $tempdir = qw(/home/sco/volatile);
my $template="genbankpmXXXXX";


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
  return($self);
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

# {{{ sub head hash(name, description). Returns a list.
sub head {
my $self = shift(@_);
my %args = @_;
my @head = ("track", "name=$args{name}", qq/description="$args{description}"/,
qq/useScore=0/, qq/visibility=1/, qq/itemRgb=on/, qq/color=100/);
return(@head);
}
# }}}


# {{{ sub feature hash(molecule, id, strand, start, end, bpcoords, rgb|color).
# Returns a list.
# bpcoords is boolean. Pass true if coordinates are bioperl coordinates.
# If bpcoords is true strand should be specified as 1 or -1.
sub feature {
my $self = shift(@_);
my %args = @_;
my $mol = $args{molecule};
my $id = $args{id};
my $score = 0;
my $strand = $args{strand};
my $rgb = $args{rgb} ? $args{rgb} : $args{color};

unless($rgb) {
if($strand eq '+' or $strand == 1) {
$rgb = qq/63,255,63/;
}
else {
$rgb = qq/63,63,255/;
}
}

my $start;
if($args{bpcoords}) {
$strand = $args{strand} == 1 ? '+' : '-';
$start = $args{start} - 1;
}
my $end = $args{end};

my @retlist = (
$mol, $start, $end, $id, $score, $strand, $start, $end, $rgb
);

return(@retlist);
}
# }}}

# {{{ sub featHash hash(bedfile)
# The returned coordinates are bioperl coordinates. i.e. Start is
# incremented by 1.
sub featHash {
my $self = shift(@_);
my %args = @_;
my %rethash;
open(my $bedfh, "<", $args{bedfile});
my $header = readline($bedfh);

while(my $line = readline($bedfh)) {
chomp($line);
my @ll=split(/\t/, $line);
my $mol = $ll[0];
my $start = $ll[1];
my $end = $ll[2];
my $locus_tag = $ll[3];
my $lt = $locus_tag;
my $strand = $ll[5] eq '+' ? 1 : -1;

$rethash{$lt} = {
  molecule => $mol,
  start => $start + 1,
  end => $end,
  locus_tag => $lt,
  lt => $lt,
  strand => $strand
}
}
close($bedfh);
return(%rethash);
}
# }}}

# {{{ sub featList
# The returned coordinates are bioperl coordinates. i.e. Start is
# incremented by 1.
sub featList {
my $self = shift(@_);
my %args = @_;
my @retlist;
open(my $bedfh, "<", $args{bedfile});
my $header = readline($bedfh);

while(my $line = readline($bedfh)) {
chomp($line);
my @ll=split(/\t/, $line);
my $mol = $ll[0];
my $start = $ll[1];
my $end = $ll[2];
my $locus_tag = $ll[3];
my $lt = $locus_tag;
my $strand = $ll[5] eq '+' ? 1 : -1;

push(@retlist, {
  molecule => $mol,
  start => $start + 1,
  end => $end,
  locus_tag => $lt,
  lt => $lt,
  strand => $strand
});

}
close($bedfh);

my @sorted = sort {$a->{start} <=> $b->{start}} @retlist;


return(@sorted);
}
# }}}


return(1);

__END__
