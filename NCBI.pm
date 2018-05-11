package Sco::NCBI;
use 5.14.0;
use Carp;
use File::Copy;
use Bio::SeqIO;
use Bio::Seq;
use XML::Simple;
use LWP::Simple;
use Bio::SeqIO;
use Bio::Seq;
use Net::FTP;
use File::Temp qw(tempfile tempdir);
my $template='ScoNCBIXXXXX';
my $dir='/home/sco/volatile';
#$ENV{http_proxy}='http://wwwcache.bbsrc.ac.uk:8080';

our $AUTOLOAD;
use lib qw(/home/sco /home/sco/perllib);
use Scoglobal;
our @ISA = qw(Scoglobal);


# {{{ Class variables
my $esearchURL='http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?';
my $efetchURL='http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?';
my $esummaryURL='http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esummary.fcgi?';
my $ftp=Net::FTP->new(Host => "ftp.ncbi.nih.gov", Passive => 1);
my $draftDir = qq(genomes/ASSEMBLY_BACTERIA);
my $finishedDir = qq(genomes/Bacteria);
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
  $self->{draftDir} = $draftDir;
  $self->{finishedDir} = $finishedDir;
  bless($self, $class);
  return($self);
}
# }}}

# {{{ sub genus hash(type => finished or draft, genus). Returns a list.
# type is essential.
# genus is optional
sub genus {
  my $self = shift(@_);
  my %args = @_;
  my @retlist;
  my $genus = $args{genus};
  my $ford = $args{type}; # finished or draft
  $ftp->login('anonymous', 'govind.chandra@jic.ac.uk');
  if($ford =~ m/^f/i) {
  $ftp->cwd($finishedDir);
  }
  if($ford =~ m/^d/i) {
  $ftp->cwd($draftDir);
  }
  my @listing = $ftp->ls();
  if($genus) {
    for my $el (@listing) {
      if($el =~ m/^$genus/) {
        push(@retlist, $el);
      }
    }
  }
  else { @retlist = @listing; }
  $ftp->quit();
  return(@retlist);
}
# }}}


# {{{ sub genbankStats (Bio::Seq or a filename) returns(hash()) ### for nucleotide genbanks only
# keys in rethash(accession, binomial, taxonomy, n_cds, n_rrna, n_trna, n_psgene, n_pscds
# description, size, linear 
sub genbankStats {
  my $self = shift(@_);
  my $temp = shift(@_);
  my $seqobj;
  if(ref($temp)) { $seqobj = $temp; }
  elsif (-e $temp) {
    my $seqio = Bio::SeqIO->new(-file => $temp);
    $seqobj = $seqio->next_seq();
  }
  else {
    croak("Argument is neither a Bio::Seq nor an existing filename\n");
  }
  my %rh;  ### the returned hash
  $rh{accession} = $seqobj->accession();
  my ($species, @cl, $taxonomy, $binomial);
  $species=$seqobj->species();
  if($species) {
    @cl = $species->classification();
    $taxonomy = join(" : ", reverse(@cl));
    $binomial =$species->binomial();
  }
  else {
    $taxonomy = undef;
    $binomial = undef;
  }
  $rh{binomial} = $binomial;
  $rh{taxonomy} = $taxonomy;
  $rh{nt_len}=$seqobj->length();
  $rh{n_cds} = 0; 
  $rh{n_rrna} = 0; 
  $rh{n_trna} = 0;
  $rh{n_psgene} = 0;
  $rh{n_pscds} = 0;
  foreach my $feature ($seqobj->all_SeqFeatures()) {
    if($feature->primary_tag() eq 'CDS') {
      $rh{n_cds}+=1;
      if($feature->has_tag('pseudo')){
        $rh{n_pscds}+=1;
      }
    }
    elsif($feature->primary_tag() eq 'gene') {
      if($feature->has_tag('pseudo')){
        $rh{n_psgene}+=1;
      }
    }
    elsif($feature->primary_tag() eq 'tRNA') {
      $rh{n_trna}+=1;
    }
    elsif($feature->primary_tag() eq 'rRNA') {
      $rh{n_rrna}+=1;
    }
  }
  $rh{description}=$seqobj->description();
  $rh{size}=$seqobj->length();

  my @temp = Scoglobal->gc_frac($seqobj->seq());
  $rh{gc_frac} = sprintf("%.3f", $temp[2]);
  if($seqobj->is_circular()) {
    $rh{linear} = 0;
  }
  else {
    $rh{linear} = 1;
  }
#  $rh{description}=~s/[^\']\'[^\']/\'\'/g;
  return(%rh);
}

# }}}


# {{{ sub fetchGbk hash(uid, outfile) returns(Bio::Seq object) ### for nucleotide genbanks only
# if outfile is specified we get a file, otherwise we get a Bio::Seq object.
sub fetchGbk {
my $self = shift(@_);
my %args = @_;
my $uid = $args{uid};
my $outfile = $args{outfile};
my $efetch= $efetchURL;
$efetch .= 'db=nucleotide&';
$efetch .= "id=$uid\&";
$efetch .= "retmode=text\&";
$efetch .= "rettype=gbwithparts";
my ($fh, $filename)=tempfile($template, DIR => $dir, SUFFIX => '.gbk');
select($fh);
getprint($efetch);
select(STDOUT);
close($fh);
if($outfile) {
if(copy($filename, $outfile)) {
  unlink($filename);
  return(1);
}
else {
  carp("Copy of $filename to $outfile failed");
  unlink($filename);
  return(0);
}
}
else {
my $seqio = Bio::SeqIO->new(-file  => $filename);
my $seqobj = $seqio->next_seq();
unlink($filename);
return($seqobj);
}

}
# }}}


# {{{ sub fetchSummary hash(uid, outfile) returns(name of an xml file)
sub fetchSummary {
my $self = shift(@_);
my %args = @_;
my $uid = $args{uid};
my $outfile = $args{outfile};
my $efetch= $esummaryURL;
$efetch .= 'db=nucleotide&';
$efetch .= "id=$uid\&";
$efetch .= "version=2.0";
my ($fh, $filename)=tempfile($template, DIR => $dir, SUFFIX => '.xml');
select($fh);
getprint($efetch);
select(STDOUT);
close($fh);
if($outfile) {
if(copy($filename, $outfile)) {
  unlink($filename);
  return($outfile);
}
else {
  carp("Copy of $filename to $outfile failed");
  unlink($filename);
  return(0);
}
}
else {
  return($filename);
}
}
# }}}

# {{{ sub ntAccession2UID (nucleotide_accession) returns (UID) {
sub ntAccession2UID {
  my $self = shift(@_);
  my $accn = shift(@_);
my $esearch= $esearchURL;
my $term = "term\=$accn\[ACCN\]";
my $retmax='retmax=20';
my $retmode='retmode=xml';
my $database='db=nucleotide';
$esearch .= $term .'&' . $database . '&' . $retmax . '&' . $retmode;
#$esearch .= $term .'&' . $retmode;
#print(STDERR "$esearch\n");
my ($fh, $filename)=tempfile($template, DIR => $dir);
getstore($esearch, $filename);
my $xmlsim=XML::Simple->new();
my $tree=$xmlsim->XMLin($filename);
unlink($filename);
#my @uids=@{$tree->{IdList}->{Id}};
my $uid=$tree->{IdList}->{Id};
return($uid);
}
# }}}

# {{{ sub accession2UID hash(accession, database) returns (UID) {
sub accession2UID {
  my $self = shift(@_);
  my %args = @_;
  my $accn = $args{accession};
  my $database=q/db=nucleotide/;
  if(exists($args{database})) {
    $database=qq/db=$args{database}/;
  }
my $esearch= $esearchURL;
my $term = "term\=$accn\[ACCN\]";
my $retmax='retmax=20';
my $retmode='retmode=xml';
$esearch .= $term .'&' . $database . '&' . $retmax . '&' . $retmode;
#$esearch .= $term .'&' . $retmode;
#print(STDERR "$esearch\n");
my ($fh, $filename)=tempfile($template, DIR => $dir);
getstore($esearch, $filename);
my $xmlsim=XML::Simple->new();
my $tree=$xmlsim->XMLin($filename);
unlink($filename);
#my @uids=@{$tree->{IdList}->{Id}};
my $uid=$tree->{IdList}->{Id};
return($uid);
}
# }}}

return(1);

