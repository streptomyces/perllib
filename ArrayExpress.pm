package Sco::ArrayExpress;
use common::sense;
use Carp;
use File::Copy;
#use Bio::SeqIO;
#use Bio::Seq;
#use XML::Simple;
use LWP::Simple;
use File::Temp qw(tempfile tempdir);
my $template='ScoAExpXXXXX';
my $dir='/home/sco/volatile';
#$ENV{http_proxy}='http://wwwcache.bbsrc.ac.uk:8080';

our $AUTOLOAD;
use lib qw(/home/sco /home/sco/perllib);
use Scoglobal;
our @ISA = qw(Scoglobal);

# {{{ Class variables
my $aeUrl='http://www.ebi.ac.uk/arrayexpress';
# }}}



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

# {{{ sub readADF hash(file). Returns a hash.
sub readADF {
my $self = shift(@_);
my %args = @_;
open(my $adfh, "<", $args{file});
my $relFlag = 0;
my %retHash;
while(<$adfh>) {
my $line = $_;
chomp($line);
if($line=~m/^Reporter Name/) { $relFlag = 1; next;}
if($relFlag) {
  my @llist = split(/\t/, $line);
  my $key = shift(@llist);
  $retHash{$key} = [ @llist ];
}
}
close($adfh);
return(%retHash);
}
# }}}

# {{{ sub readSDRF hash(file). Returns a hash.
sub readSDRF {
my $self = shift(@_);
my %args = @_;
open(my $sdrfh, "<", $args{file});
my $discard = readline($sdrfh);
#my $relFlag = 0;
my %retHash;
while(<$sdrfh>) {
my $line = $_;
chomp($line);
  my @llist = split(/\t/, $line);
  my $key = shift(@llist);
  $retHash{$key} = [ @llist ];
}
close($sdrfh);
return(%retHash);
}
# }}}

# {{{ sub id2names hash(filename file, hashref adf);
sub id2names {
my $self = shift(@_);
my %args = @_;
my $file = $args{file};
my $adf = $args{adf};
my %rethash;
open(my $infile, "<", $file);
while(<$infile>) {
  chomp();
  my @llist = split(/\t/, $_);
  if($adf->{$llist[0]}->[1]) {
    $rethash{$adf->{$llist[0]}->[1]} = [@llist];
#    splice(@llist, 0, 1, $adf->{$llist[0]}->[1]);
#    push(@retlist, join("\t", @llist));
  }
}
close($infile);
return(%rethash);
}
# }}}

return(1);

