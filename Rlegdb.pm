package Sco::Rlegdb;
use strict qw(vars subs);
use lib qw(/home/sco /home/sco/perllib);
use Sco::Global;
our @ISA = ("Sco::Global");
use DBI;
use Carp;
use Bio::SeqIO;
use Bio::Seq;
use Bio::SearchIO;
use File::Copy;
use File::Temp qw(tempfile tempdir);
use File::Basename;
my $tempdir = qw(/home/sco/volatile);
my $template="rlegdbXXXXX";

my $dbname = "rleg";
my $dbhost = "jiilin9.jic.bbsrc.ac.uk";
my $handle=DBI->connect("DBI:Pg:dbname=$dbname;host=$dbhost", 'sco', 'tsinH4x');


our $draft_dir = '/home/sco/seq/nt/draft_ftp';
our $finished_dir = '/home/sco/seq/nt/genbank_ftp';

# {{{ sub new
sub new {
  my($class, $self, %argv);
  $class=shift(@_);
  %argv = @_;
  $self={};
  foreach my $key (keys(%argv)) {
    $self->{$key} = $argv{$key};
  }
  bless($self, $class);
  return($self);
}
# }}}

# {{{ sub flankingFeats (hash(pre, start, end)) returns list of [lt, dist, product].
sub flankingFeats {
my $self = shift(@_);
my %args = @_;
my $pre = $args{pre};
my @retlist;

my $iqstr = qq/select accession from organisms where pre = '$pre'/;
my ($accession) = $handle->selectrow_array($iqstr);

my ($qstr, $lt, $start_pos, $end_pos, $product, $dist);

$qstr = qq/select id, start_pos, end_pos, product from features where accession = '$accession' and strand = -1 and end_pos <= $args{start} order by end_pos desc limit 1/;
if(($lt, $start_pos, $end_pos, $product) = $handle->selectrow_array($qstr)) {
$dist = $args{start} - $end_pos;
push(@retlist, [$lt, $dist, $product]);
}
else {
print(STDERR "$qstr\n");
}

$qstr = qq/select id, start_pos, end_pos, product from features where accession = '$accession' and strand = 1 and start_pos >= $args{end} order by start_pos asc limit 1/;
if(($lt, $start_pos, $end_pos, $product) = $handle->selectrow_array($qstr)) {
$dist = $start_pos - $args{end};
push(@retlist, [$lt, $dist, $product]);
}
else {
  print(STDERR "$qstr\n");
}
return(@retlist);
}
# }}}

# {{{ sub inCDS (hash(pre, start, end)) returns hash(lt, start_pos, end_pos, strand, product)
sub inCDS {
my $self = shift(@_);
my %args = @_;
my $pre = $args{pre};
my @retlist;

my $iqstr = qq/select accession from organisms where pre = '$pre'/;
my ($accession) = $handle->selectrow_array($iqstr);

my ($qstr, $lt, $start_pos, $end_pos, $product, $dist, $strand);

$qstr = qq/select id, start_pos, end_pos, strand, product from features where accession = '$accession' and ( (start_pos <= $args{end} and end_pos >= $args{end}) or (start_pos <= $args{start} and end_pos >= $args{start}) )/;
if(($lt, $start_pos, $end_pos, $strand, $product) = $handle->selectrow_array($qstr)) {
return(locus_tag => $lt, start_pos => $start_pos, end_pos => $end_pos, strand => $strand, product=> $product);
}
else {
return();
}
}
# }}}

# {{{ sub DESTROY
sub DESTROY {
  my $self = shift(@_);
  $handle->disconnect();
}
# }}}

