use strict;
use XML::Simple;
use LWP::Simple;
use Bio::SeqIO;
use File::Temp qw(tempfile tempdir);
my $template='whiasearchXXXXX';
my $dir='/home/sco/volatile';

# {{{ Getopt::Long stuff
use Getopt::Long;
my $outfile;
GetOptions (
"outfile:s" => \$outfile
);

open(my $ofh, ">", $outfile);

my $joined_flag;


# }}}


my $esearch='http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?';
my @query;
push(@query, 'term=streptomyces[ORGN] AND cds[FKEY] AND 10000:10000000[SLEN]');
push(@query, 'retmax=200000');
push(@query, 'retmode=xml');
push(@query, 'db=nucleotide');
my $qstr = join('&', @query);
$esearch .= $qstr;
#$esearch .= $term .'&' . $retmode;
$esearch =~ s/\s/\+/g;
print(STDERR "$esearch\n");

my ($fh, $filename)=tempfile($template, DIR => $dir);
getstore($esearch, $filename);
my $xmlsim=XML::Simple->new();
my $tree=$xmlsim->XMLin($filename);
my @uids=@{$tree->{IdList}->{Id}};
print(STDERR "UID count is: ", scalar(@uids), "\n");

my $cnt=0;
foreach my $uid (@uids) {
$cnt+=1;
#print("$cnt $uid\n");
}

system("cp $filename temp.xml");

#print($filename);
unlink($filename);

# exit;

# {{{ fetch genpept for all uids in @uids
my $fetch_cnt = 0;
foreach my $uid (@uids) {
my $efetch='http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?';
$efetch .= 'db=nucleotide&';
$efetch .= "id=$uid\&";
$efetch .= "retmode=text\&";
$efetch .= "rettype=gbwithparts";

#my ($fh, $filename)=tempfile($template, DIR => $dir);
my $oldh = select($ofh);
getprint($efetch);
close($ofh);
select $oldh;

# my $protlen = &genpept_len($filename);

#unless(unlink($filename)) {
#  print(STDERR "Failed to unlink $filename");
#}

$fetch_cnt+=1;
print(STDERR "$fetch_cnt\r");
if($fetch_cnt > 1) {last;}
}

# }}}


exit;

### subs begin ###

# {{{ sub codedby
sub codedby {
my $infile = shift(@_);
my $uid = shift(@_);

my $seqio=Bio::SeqIO->new(-file => $infile);
my $seqobj=$seqio->next_seq();
my $coded_by;
my $has_coded_by = 0;
  foreach my $feature ($seqobj->all_SeqFeatures()) {
    if($feature->primary_tag() eq 'CDS' and $feature->has_tag('coded_by')) {
          $coded_by=join(" ", $feature->get_tag_values('coded_by'));
          $has_coded_by = 1;
        }
      }
unless($has_coded_by) {return();}

#print(STDERR "$uid\t$coded_by\n");
my $acc;
my $start;
my $end;
my $ntuid;
my $strand = 1;

  if($coded_by=~m/complement/) {
    $strand = 2;
    $coded_by=~s/complement\(//;
    $coded_by=~s/\)$//;
  }

if($coded_by=~m/join/i) {
  $coded_by=~s/join\(//;
  $coded_by=~s/\)//;
  $joined_flag = 1;
}
else {
  $joined_flag = 0;
}

  my @ranges = split(/\,\s*/, $coded_by);

  ($acc, undef) = split(/:/, $ranges[0]);

  my @ranges2;
  foreach my $range (@ranges) {
    $range=~s/$acc://;
    push(@ranges2, $range);
  }

  my @ends;
  foreach my $r2 (@ranges2) {
    $r2=~s/<|>//g;
    push(@ends, split(/\.\./, $r2));
  }
  

my @sortends = sort {$a <=> $b} @ends;
$ntuid = &accn_to_uid($acc);
$start = $sortends[0];
$end = $sortends[$#sortends];

#  print(STDERR "R2 ", join(" ", @ranges2), "\n");
#  print(STDERR "Ends ", join(" ", @ends), "\n");
#  print(STDERR "Sortends ", join(" ", @sortends), "\n");


#close(INFILE);
#print(STDERR "sub coded_by: $uid\t$acc\t$ntuid\t$start\t$end\t$strand\n");
#print(STDERR "\n\n");
return($ntuid, $start, $end, $strand);
#  last;
}
# }}}


# {{{ sub accn_to_uid {
sub accn_to_uid {
  my $accn = shift(@_);
my $esearch='http://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi?';

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


#{{{ sub subseq_ntuid {
sub subseq_ntuid {
  my $ntuid = shift(@_);
  my $start = shift(@_);
  my $stop = shift(@_);
  my $strand = shift(@_);
my $efetch='http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?';
$efetch .= 'db=nucleotide&';
$efetch .= "id=$ntuid\&";
$efetch .= "seq_start=$start\&";
$efetch .= "seq_stop=$stop\&";
$efetch .= "strand=$strand\&";
$efetch .= "complexity=1\&";
$efetch .= "retmode=text\&";
$efetch .= "rettype=fasta";


my $length = $stop - $start + 1;
#print(STDERR "$uid\n$efetch\n");
#print(STDERR "$length\t$efetch\n");
my $ntfas = get($efetch);
if($joined_flag) {
#print(FNA "$ntfas\n");
}
return($length);
}
# }}}

# {{{ sub genpept_len {
sub genpept_len {
my $infile = shift(@_);
my $seqio=Bio::SeqIO->new(-file => $infile);
my $seqobj=$seqio->next_seq();
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
if($taxonomy) {
#print(STDERR "$taxonomy\n");
if($taxonomy=~m/actinobacteria/i) {
#$actino_out->write_seq($seqobj);
}
else {
#$other_out->write_seq($seqobj);
}
}
else {
print(STDERR "No taxonomy\n");
}

return($seqobj->length());
}
# }}}

__END__

wget -O - 'http://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=protein&id=254593662&retmode=text&rettype=fasta' 


