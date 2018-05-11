package Sco::Pubmed;
use common::sense;
use Carp;
use XML::Simple;
use LWP::Simple;
use Data::Dumper;
use File::Temp qw(tempfile tempdir);
my $template='ScopbmdXXXXX';
my $tempdir='/home/sco/volatile';
#$ENV{http_proxy}='http://wwwcache.bbsrc.ac.uk:8080';
our $AUTOLOAD;
use lib qw(/home/sco /home/sco/perllib);
use Sco::Global;
our @ISA = qw(Sco::Global);

my $esearchURL='http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?';
my $efetchURL='http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?';
my $pubmedURL = "http://www.ncbi.nlm.nih.gov/pubmed";


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

### more subs here ###

# {{{ somesub (args) returns(retvals)
sub somesub {
my $self=shift(@_);
my %args=@_;
return('something'); # change this
}
# }}}

# {{{ sub getXML (list of UIDs) returns(a list of string refs);
sub getXML {
  my $self = shift(@_);
  my @uids = @_;
  my @retlist;
foreach my $uid (@uids) {
my $efetch .= $efetchURL . 'db=pubmed&';
$efetch .= "id=$uid\&";
$efetch .= "retmode=xml\&";
#$efetch .= "retmode=text\&";
$efetch .= "rettype=abstract";
#print(STDERR "$count\r");
my $abstract = get($efetch);
  my $xml = get($efetch);
  push(@retlist, \$xml);
}
return(@retlist);
}
# }}}

# {{{ sub getUIDs (hash(days, term)) returns(a listref of UIDs)
sub getUIDs {
  my $self = shift(@_);
  my %args = @_;
  my $days;
  if(exists($args{days})) {
    $days = $args{days};
  }
  else { $days = 7;}
  my $term = $args{term};
  if ($days == 0) {
    $term .= " ";
  }
  else {
    $term .= " AND \"last $days days\"[DP]";
  }
  my $retmax='retmax=40000';
  my $retmode='retmode=xml';
  my $reldate="reldate=$days";
  my $esearch .= $esearchURL . $term . '&' . $retmax . '&' . $reldate . '&' . $retmode;
#  print(STDERR "$esearch\n");
  my($fh, $fn)=tempfile($template, DIR => $tempdir, SUFFIX => '.xml');
#close($fh);
  getstore($esearch, $fn);
#getprint($esearch);
#print($fn);
  my $xs = XML::Simple->new(ForceArray => 1);
  my $ref = $xs->XMLin($fn);
  unlink($fn);
#print(Dumper($ref), "\n\n");
  my $idref = $ref->{IdList}->[0]->{Id};
  return($idref);
}
# }}}

# {{{ sub getAbstract (list of UIDs) returns(a list of hashrefs);
sub getAbstracts {
  my $self = shift(@_);
  my @uids = @_;
  my @retlist;
foreach my $uid (@uids) {
my $efetch .= $efetchURL . 'db=pubmed&';
$efetch .= "id=$uid\&";
$efetch .= "retmode=xml\&";
#$efetch .= "retmode=text\&";
$efetch .= "rettype=abstract";
#print(STDERR "$count\r");
my $abstract = get($efetch);
  my($fh, $fn)=tempfile($template, DIR => $tempdir, SUFFIX => '.xml');
#close($fh);
  getstore($efetch, $fn);
#getprint($esearch);
#print($fn);
  my $xs = XML::Simple->new(ForceArray => 1);
  my $tree = $xs->XMLin($fn);
#  print(Dumper($tree),"\n"); last;
  unlink($fn);
my $pmid=$tree->{PubmedArticle}->[0]->{MedlineCitation}->[0]->{PMID}->[0]->{content};

my $temp = $tree->{PubmedArticle}->[0]->{MedlineCitation}->[0]->{DateCreated};
# print(STDERR Dumper(\$temp), "\n");
my $cyear = $tree->{PubmedArticle}->[0]->{MedlineCitation}->[0]->{DateCreated}->[0]->{Year}->[0];
my $cmonth = $tree->{PubmedArticle}->[0]->{MedlineCitation}->[0]->{DateCreated}->[0]->{Month}->[0];
my $cday = $tree->{PubmedArticle}->[0]->{MedlineCitation}->[0]->{DateCreated}->[0]->{Day}->[0];

# print(STDERR "--- $cyear   $cmonth   $cday ---\n");
$cyear=~s/^0+//;
$cmonth=~s/^0+//;
$cday=~s/^0+//;

#print ("$pmid   $uid\n");
my $xarticle=$tree->{PubmedArticle}->[0]->{MedlineCitation}->[0]->{Article};
my $title = $xarticle->[0]->{ArticleTitle}->[0];
my $abstractref = $xarticle->[0]->{Abstract}->[0]->{AbstractText}; # always a listref
#my $abstract = $xarticle->[0]->{Abstract}->[0]->{AbstractText}->[0];
# print(STDERR Dumper(\$abstractref), "\n");
my $abstract = &_abstract($abstractref);
my $lastName = $xarticle->[0]->{AuthorList}->[0]->{Author}->[0]->{LastName}->[0];
my $foreName = $xarticle->[0]->{AuthorList}->[0]->{Author}->[0]->{ForeName}->[0];
my $author = join(" ", $foreName, $lastName);
my $journalTitle = $xarticle->[0]->{Journal}->[0]->{Title}->[0];
my $year = $xarticle->[0]->{Journal}->[0]->{JournalIssue}->[0]->{PubDate}->[0]->{Year}->[0];
my $month = $xarticle->[0]->{Journal}->[0]->{JournalIssue}->[0]->{PubDate}->[0]->{Month}->[0];
my $day = $xarticle->[0]->{Journal}->[0]->{JournalIssue}->[0]->{PubDate}->[0]->{Day}->[0];
my $journalDate = "$day $month $year";  
my $authors = &_authors($xarticle);
my $pubtype = &_pubtype($xarticle);
push(@retlist, { pmid => $pmid,
authors => $authors,
pubtype => $pubtype,
title => $title, 
abstract => $abstract,
journal => $journalTitle,
year => $cyear,
month => $cmonth,
day => $cday}
);
}
return(@retlist);
}
# }}}




# {{{ sub _abstract
sub _abstract {
my $refAt = shift(@_); # This is always a listref.
my @stringbuild;
foreach my $temp (@{$refAt}) {
if(ref($temp)) {
push(@stringbuild, "$temp->{Label}: $temp->{content}");
}
else {push(@stringbuild, $temp);}
}
#print(STDERR join(" ", @stringbuild), "\n");
return(join(" ", @stringbuild));
}
# }}}

# {{{ sub _authors {
sub _authors {
my $xa = shift(@_);
my @alist = @{$xa->[0]->{AuthorList}->[0]->{Author}};
#print(STDERR Dumper(@alist));
my @names;
foreach my $auth (@alist) {
#  print(STDERR ref($auth), "\n");
#  print(STDERR keys(%{$auth}), "\n");
my $lname = $auth->{LastName}->[0];
my $fname = $auth->{ForeName}->[0];
my $name = join(" ", ($fname, $lname));
#print(STDERR "--- $name ---");
push(@names, $name);
}
return(join(", ", @names));
}
# }}}

# {{{ sub _pubtype {
sub _pubtype {
my $xa = shift(@_);
my @plist = @{$xa->[0]->{PublicationTypeList}->[0]->{PublicationType}};
#print(STDERR Dumper(\@plist));
my @types;
foreach my $pt (@plist) {
push(@types, $pt);
}
return(join(", ", @types));
}
# }}}

return(1);

