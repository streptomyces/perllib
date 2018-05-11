package Sco::Orthologs2011;
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
my $template="orthoXXXXX";

my $dbname = "ortho2011";
my $dbhost = "jiilin9.jic.bbsrc.ac.uk";
my $handle=DBI->connect("DBI:Pg:dbname=$dbname;host=$dbhost", 'sco', 'tsinH4x');


our $draft_dir = '/home/sco/seq/nt/draft_ftp';
our $finished_dir = '/home/sco/seq/nt/genbank_ftp';

my $clustersTable = 'gbstats';
my $clustersTTATable = 'tta';
my $clustersView = 'vtta';
my $ttaPosnsTable = 'proteins';

my %strepDBaccessions;
$strepDBaccessions{sco} =  "AL645882";
$strepDBaccessions{scp1} = "AL589148";
$strepDBaccessions{scp2} = "AL645771";
 
my $rpsBlastDir = '/home/sco/blast_databases/rpsblast';
our $blastdbdir = '/home/sco/customers/orthologs/bldb';

# {{{

=head2 Tables

bigortho=> \d proteins
     Table "public.proteins"
  Column   |  Type   | Modifiers 
----------+---------+-----------
 filename  | text    | not null
 serial    | integer | not null
 nt_seq    | text    | not null
 aa_seq    | text    | not null
 identical | text    | 
 shahex    | text    | 
 id        | text    | 
Indexes:
    "proteins_filename_key" UNIQUE, btree (filename, serial)
    "protid" btree (id)

bigortho=> \d features
     Table "public.features"
  Column   |  Type   | Modifiers 
-----------+---------+-----------
 filename  | text    | not null
 serial    | integer | not null
 pri_tag   | text    | not null
 start_pos | integer | not null
 end_pos   | integer | not null
 strand    | integer | 
 locus_tag | text    | 
 product   | text    | 
Indexes:
    "features_filename_key" UNIQUE, btree (filename, serial)

        Table "public.recihits"
 Column  |       Type       | Modifiers 
---------+------------------+-----------
 query   | text             | not null
 hit     | text             | 
 frac_id | double precision | 
 qcov    | double precision | 
 hcov    | double precision | 
 signif  | double precision | 
 bits    | double precision | 
Indexes:
    "recihits_query_key" UNIQUE, btree (query, hit)

=cut

# }}}

# {{{ bghash (hash(fnafile)) returns(a hash of mono and dinucleotide frequencies)
sub bghash {
my $self = shift(@_);
my %args = @_;
my $seqio = Bio::SeqIO->new(-file => $args{fnafile});
my %counts;
my $total;
while(my $seqobj = $seqio->next_seq()) {
  my $seq = lc($seqobj->seq());
  my $pos = 0;
  while(my $dinuc = substr($seq, $pos, 2)) {
    if($dinuc=~m/[^acgt]/) { $pos += 1; next; }
    my ($na, $nb) = split("", $dinuc);
    $counts{$dinuc} += 1;
    $counts{$na} += 1;
    $total += 1;
    $pos += 1;
  }
}
my %freqs;
foreach my $key (sort keys(%counts)) {
$freqs{$key} = $counts{$key} / $total
}
return(%freqs);
}
# }}}

# {{{ sub dbstats no arguments. Prints to STDOUT.
sub dbstats {
my $self = shift(@_);
foreach my $table (qw(orgtab features proteins recihits rnas contigs)) {
my $qstr = qq/select count(*) from $table/;
my ($count) = $handle->selectrow_array($qstr);
  print("$table\t$count\n");
}

my $qstr1 = qq/select count(distinct organism) from orgtab/;
my ($count) = $handle->selectrow_array($qstr1);
print("Orgs in orgtab\t$count\n");

$qstr1 = qq/select count(distinct organism) from orgtab where filename in
(select filename from features)/;
($count) = $handle->selectrow_array($qstr1);
print("Orgs in features\t$count\n");
}
# }}}

# {{{ sub flankingFeats (hash(filename, start, end)) returns list of [lt, dist].
sub flankingFeats {
my $self = shift(@_);
my %args = @_;
my $filename = $args{filename};
my @retlist;
my ($qstr, $lt, $start_pos, $end_pos, $product, $dist);

$qstr = qq/select locus_tag, start_pos, end_pos, product from features where filename = '$filename' and strand = -1 and end_pos <= $args{start} order by end_pos desc limit 1/;
if(($lt, $start_pos, $end_pos, $product) = $handle->selectrow_array($qstr)) {
$dist = $args{start} - $end_pos;
push(@retlist, [$lt, $dist, $product]);
}
else {
print(STDERR "$qstr\n");
}

$qstr = qq/select locus_tag, start_pos, end_pos, product from features where filename = '$filename' and strand = 1 and start_pos >= $args{end} order by start_pos asc limit 1/;
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

# {{{ sub organism_fna
sub organism_fna {
  my $self = shift(@_);
  my %args = @_;
  my $org = $args{organism};
  my @files = $self->organism2files($org);
  my $fh;
  open($fh, ">$args{outfile}");
    my $seqout = Bio::SeqIO->new(-fh => $fh, -format => 'fasta');
  foreach my $file (@files) {
    my ($tfn, $tdir, $ext)= fileparse($file, qr/\.[^.]*/);
    my $fn = $tfn . $ext;
    my $qstr = qq/select nt_seq from contigs where filename = '$fn'/;
    my ($nt_seq) = $handle->selectrow_array($qstr);
    my $outobj = Bio::Seq->new(-seq => $nt_seq);
    $outobj->display_name($fn);
    $seqout->write_seq($outobj);
  }
  close($fh);
  return(1);
}
# }}}

# {{{ sub ttaposns (id) returns a list.
sub ttaposns {
  my $self = shift(@_);
  my $id = shift(@_);
  my $qstr = "select ttapos from $ttaPosnsTable where id = '$id'";
  my $stmt = $handle->prepare($qstr);
  $stmt->execute();
  my $hr = $stmt->fetchrow_hashref();
  return(@{$hr->{ttapos}});
}
# }}}

# {{{ listorganisms return(list of distinct organism names)
sub listorganisms {
  my $self = shift(@_);
  my $organism = shift(@_);
  my $qstr = "select distinct organism from orgtab where n_cds > 0 order by organism";
  my $stmt = $handle->prepare($qstr);
  $stmt->execute();
  my @orgnames;
  while(my $hr = $stmt->fetchrow_hashref()) {
    my $org = $hr->{organism};
    push(@orgnames, $org);
  }
  $stmt->finish();
  return(@orgnames);
}
# }}}

# {{{ organism2files_nodir (organism) return(list of filenames)
sub organism2files_nodir {
  my $self = shift(@_);
  my $organism = shift(@_);
  my $qstr = "select filename from orgtab where organism = '$organism' order by nt_len DESC";
  my $stmt = $handle->prepare($qstr);
  $stmt->execute();
  my @filenames;
  while(my $hr = $stmt->fetchrow_hashref()) {
    my $fn = $hr->{filename};
    push(@filenames, $fn);
  }
  $stmt->finish();
  return(@filenames);
}
# }}}

# {{{ genus2files_nodir (genus) return(list of filenames)
sub genus2files_nodir {
  my $self = shift(@_);
  my $genus = shift(@_);
  my $qstr = "select filename from orgtab where organism ~* '^$genus' order by nt_len DESC";
  my $stmt = $handle->prepare($qstr);
  $stmt->execute();
  my @filenames;
  while(my $hr = $stmt->fetchrow_hashref()) {
    my $fn = $hr->{filename};
    push(@filenames, $fn);
  }
  $stmt->finish();
  return(@filenames);
}
# }}}

# {{{ organism2files (organism) return(list of filenames)
sub organism2files {
  my $self = shift(@_);
  my $organism = shift(@_);
  my $qstr = "select dorf, filename from orgtab where organism like '%$organism%' order by nt_len DESC";
  my $stmt = $handle->prepare($qstr);
  $stmt->execute();
  my @filenames;
  while(my $hr = $stmt->fetchrow_hashref()) {
    my $fn = $hr->{filename};
    my $dorf = $hr->{dorf};
    my $seqdir = $dorf=~m/^f/i ? $finished_dir : $draft_dir;
    my $filename;
    $filename = $seqdir . '/' . $fn;
    push(@filenames, $filename);
  }
  $stmt->finish();
  return(@filenames);
}
# }}}

# {{{ taxclass (organism) returns a string of taxonomic class. Used for CSS class.
sub taxclass {
  my $self = shift(@_);
  my $org = shift(@_);
  my @taxo = $self->organism2taxo($org);
  if (grep {$_=~m/^corynebacterineae/i} @taxo) {
    return("coryne");
  }
  elsif (grep {$_=~m/^streptosporangineae/i} @taxo) {
    return("strsp");
  }
  elsif (grep {$_=~m/^propionibacterineae/i} @taxo) {
    return("propi");
  }
  elsif (grep {$_=~m/^pseudonocardineae/i} @taxo) {
    return("pnoca");
  }
  elsif (grep {$_=~m/^frankineae/i} @taxo) {
    return("frank");
  }
  elsif (grep {$_=~m/^catenulisporineae/i} @taxo) {
    return("catenuli");
  }
  elsif (grep {$_=~m/^streptomycineae/i} @taxo) {
    return("strepto");
  }
  elsif (grep {$_=~m/^micrococcineae/i} @taxo) {
    return("microco");
  }
  elsif (grep {$_=~m/^micromonosporineae/i} @taxo) {
    return("micromo");
  }
  elsif (grep {$_=~m/^kineosporiineae/i} @taxo) {
    return("kineo");
  }
  elsif (grep {$_=~m/^actinomycineae/i} @taxo) {
    return("actino");
  }
  elsif (grep {$_=~m/^glycomycineae/i} @taxo) {
    return("glyco");
  }
  elsif (grep {$_=~m/^bacillaceae/i} @taxo) {
    return("bsu");
  }
  elsif (grep {$_=~m/^enterobacteriaceae/i} @taxo) {
    return("eco");
  }
  else { return("hit"); }
}
# }}}

# {{{ sub TTAuidHTML

{
my $org_cnt = 0;
sub TTAuidHTML {
my $self = shift(@_);
my $uid = shift(@_);
my $qstr = "select * from $clustersView where uid = $uid";
my $stmt = $handle->prepare($qstr);
$stmt->execute();
my $last_uid = "This cannot be a UID";
my $row_cnt = 0;
my $retstr;

while(my $hr = $stmt->fetchrow_hashref()) {

my $temp = $hr->{tta_pos};
my @ttapos = @{$temp};
my @relpos = map { sprintf("%.3f", $_/$hr->{gene_len}) } @ttapos;
my $ttaPosStr = join(", ", @ttapos);
my $relPosStr = join(", ", @relpos);

if($hr->{uid} != $last_uid) {
$org_cnt += 1;
$retstr .= <<"UID";
<tr><td>$org_cnt</td><th colspan = 5>$hr->{definition}</th></tr>
<tr><td colspan = 6>Taxonomy: $hr->{taxonomy}</td></tr>
<tr><td><a href="http://www.ncbi.nlm.nih.gov/nucleotide/$hr->{uid}" target = "_blank">$hr->{uid}</a></td><td colspan = 2>$hr->{organism}</td><td>GC: $hr->{gc_frac}</td>
<td>Genes: $hr->{n_cds}</td><td>Length: $hr->{nt_len}</td></tr>
UID
$last_uid = $hr->{uid};
$row_cnt = 0;
}
$row_cnt += 1;

$retstr .= <<"ROW";
<tr><td>$row_cnt</td><td>$hr->{id}</td><td>$hr->{product}</td><td>$hr->{gene_len}</td>
<td>$ttaPosStr</td><td>$relPosStr</td>
</tr>
ROW
}

return($retstr);
}
}
# }}}

# {{{ sub clusterUIDsInOrder (columnName) returns(@uids)
sub clusterUIDsInOrder {
my $self = shift(@_);
my $ordcol = shift(@_);
my $qstr = "select uid from $clustersTable order by $ordcol";
my @retlist;
my $stmt = $handle->prepare($qstr);
$stmt->execute();
while(my $hr = $stmt->fetchrow_hashref()) {
push(@retlist, $hr->{uid});
}
return(@retlist);
}
# }}}

# {{{ sub rpsblast (%(id, db)) returns(Filename)
sub rpsblast {
  my $self = shift(@_);
  my %args = @_;
  my $id = $args{id};
#  print("$id\n");
  my $db = $args{db};
  my $rpsblastdb = $rpsBlastDir . '/' . $db;
  my($fh, $fn)=tempfile($template, DIR => $tempdir, SUFFIX => '.faa');
  close($fh);
  my $ofn = $self->id2faa($id, $fn);
  my($fh1, $fn1)=tempfile($template, DIR => $tempdir, SUFFIX => '.blast');
  close($fh1);
  `/usr/local/blastplus/bin/rpsblast -query $fn -db $rpsblastdb -evalue 10 -out $fn1`;
  unlink($fn);
#  print("$fn $ofn $fn1\n");
  return($fn1);
}
# }}}

# {{{ sub strepDBlink (locus_tag) returns(aURL) # currently for sco, scp1, and scp2 only
sub strepDBlink {
  my $self = shift(@_);
  my $locusTag = shift(@_);
  my $accession;
  if($locusTag=~m/^SCO/) {
    $accession = $strepDBaccessions{sco};
  }
  elsif($locusTag=~m/^SCP1/) {
    $accession = $strepDBaccessions{scp1};
  }
  elsif($locusTag=~m/^SCP2/) {
    $accession = $strepDBaccessions{scp2};
  }
  else { return(undef); }
  my $url = "http://strepdb.streptomyces.org.uk/cgi-bin/dc3.pl?name=$locusTag&accession=$accession";
  return($url);
}
# }}}

# {{{ sub query2List (SQLQueryString) returns(ListOfSomething)
sub query2List {
my $self = shift(@_);
my $qstr = shift(@_);
my $stmt = &_executedStmt($qstr);
my @retlist;
while(my $ar = $stmt->fetchrow_arrayref()) {
push(@retlist, [@{$ar}]);
}
return(@retlist);
}
# }}}

# {{{ internal sub  _executedStmt
sub _executedStmt {
my $qstr = shift(@_);
my $stmt = $handle->prepare($qstr);
$stmt->execute();
return($stmt);
}
# }}}

# {{{ sub organism2preList (organismName) returns(listOfPres)
sub organism2preList {
  my $self=shift(@_);
  my $org = shift(@_);
  my $qstr = "select pre from orgtab where organism = '$org'";
  my $stmt = &_executedStmt($qstr);
  my @pres;
  while(my $hr = $stmt->fetchrow_hashref()) {
    push(@pres, $hr->{pre});
  }
  return(@pres);
}
# }}}

# {{{ queriesByHitCount (list of qids) returns(list of qids)
sub queriesByHitCount {
  my $self = shift(@_);
  my @qids = @_;
  my $sorter = sub {
    my $qstra = "select count(*) from recihits where query = '$a'";
    my $qstrb = "select count(*) from recihits where query = '$b'";
    my $counta = $handle->selectrow_array($qstra);
    my $countb = $handle->selectrow_array($qstrb);
    return($countb <=> $counta);
  };
  my @oqids = sort $sorter @qids;
  return(@oqids);
}
# }}}

# {{{ sub colFromFile (colnum, filename) returns a list
sub colFromFile {
my $self = shift(@_);
my $colnum = shift(@_); # zero based
my $file = shift(@_);
my @retlist;
open(IN, "<$file");
while(<IN>) {
my $line=$_;
chomp($line);
if($line=~m/^\s*\#/ or $line=~m/^\s*$/) {next;}
my @llist=split(/\t/, $line);
push(@retlist, $llist[$colnum]);
}
close(IN);
return(@retlist);
}

# }}}

# {{{ query2hits (query_id) returns a list of hashrefs {id, frac_id, qcov, hcov} 
sub query2hits {
my $self = shift(@_);
my $qid = shift(@_);
my $qstr = "select hit, frac_id, qcov, hcov from recihits where query = '$qid'";
my $stmt = $handle->prepare($qstr);
$stmt->execute();
my @hits;
while(my $hr = $stmt->fetchrow_hashref()) {
push(@hits, { id => $hr->{hit}, frac_id => $hr->{frac_id}, qcov => $hr->{qcov}, hcov => $hr->{hcov} });
}
my @sorted = sort _sortHits @hits;
return(@sorted);
}
# }}}

# {{{ query2orghit2 (query, organism)
# returns a list of hashrefs {query, hit, frac_id, qcov, hcov} 
sub query2orghit2 {
my $self = shift(@_);
my %args = @_;
my $qid = $args{query};
my $horg = $args{organism};
my @horgfiles = $self->genus2files_nodir($horg);
# print(STDERR join("\t", @horgfiles), "\n");
my @hrl;

my $qstr = "select query, hit, frac_id, qcov, hcov from recihits where query = '$qid'";
my $stmt = $handle->prepare($qstr);
$stmt->execute();

while(my $hr = $stmt->fetchrow_hashref()) {
my (undef, $hfile, undef, $hser) = $self->id2fileser(id => $hr->{hit});
#print(STDERR "$hfile\n");
#push(@hits, $hr->{hit});
if(grep {$_ eq $hfile} @horgfiles) {
push(@hrl, $hr);
}
}
return(@hrl);
}
# }}}

# {{{ internal sub _sortHits
sub _sortHits {
if($b->{frac_id} == $a->{frac_id}) {
return($b->{qcov} <=> $a->{qcov});
}
else {
return($b->{frac_id} <=> $a->{frac_id});
}
}
# }}}

# {{{ query2orghit (query, organism) returns a  hashref {query, hit, frac_id, qcov, hcov} 
sub query2orghit {
my $self = shift(@_);
my %args = @_;
my $qid = $args{query};
my $horg = $args{organism};
my @horgfiles = $self->organism2files_nodir($horg);

my $qstr = "select hit from recihits where query = '$qid'";
my $stmt = $handle->prepare($qstr);
$stmt->execute();
my @hits;
while(my $hr = $stmt->fetchrow_hashref()) {
push(@hits, $hr->{hit});
}
foreach my $hit (@hits) {
my (undef, $hfile, undef, $hser) = $self->id2fileser(id => $hit);
if(grep {$_ eq $hfile} @horgfiles) {
my $qu = qq/select query, hit, frac_id, qcov, hcov from recihits where query = '$qid' and hit = '$hit'/;
my $st = $handle->prepare($qu); $st->execute();
my $hr = $st->fetchrow_hashref();
return($hr);
}
}
return();
}
# }}}

# {{{ id2organism (protein_id) returns a string of organism name
sub id2organism {
my $self = shift(@_);
my $id = shift(@_);
#my $qstr = "select organism from idbldb where id = '$id'";
my $qstr = "select filename, serial from proteins where id = '$id'";
my ($filename, $serial) = $handle->selectrow_array($qstr);
my $qstr2 = qq/select organism from orgtab where filename = '$filename'/;
my ($organism) = $handle->selectrow_array($qstr2);
return($organism);
}
# }}}

# {{{ organism2taxo (organism) return (a list of taxonomy);
sub organism2taxo {
my $self=shift(@_);
my $organism = shift(@_);

my $qstr="select taxonomy from orgtab where organism = '$organism'";
my ($taxonomy)=$handle->selectrow_array($qstr);
my @taxlist = split(/\;\s*/, $taxonomy);
return(@taxlist);
}
# }}}

# {{{ hitSpeciesCount (list containing references) returns(integer)
sub hitSpeciesCount {
my $self = shift(@_);
my @list = @_;
my ($totCnt, $refCnt);
foreach my $el (@list) {
$totCnt += 1;
if(ref($el)) {
$refCnt += 1;
}
}
return($refCnt);
}
# }}}

# {{{ sub lt2id (locus_tag) returns a list of ids.
sub lt2id {
  my $self = shift(@_);
  my $lt = shift(@_);
  my @retlist;
  my $qstr = "select id from proteins where locus_tag = '$lt'";
  my $stmt = $handle->prepare($qstr);
  $stmt->execute();
  while(my $hr = $stmt->fetchrow_hashref()) {
    my $id = $hr->{id};
    push(@retlist, $id);
  }
return(@retlist);
}
# }}}

# {{{ sub organismProteinIds (organism) returns a list of ids
sub organismProteinIds {
  my $self = shift(@_);
  my $organism = shift(@_);
  my @filenames = $self->organism_fn($organism);
  my @retlist;
  foreach my $filename (@filenames) {
    my $qstr1 = "select id from proteins where filename = '$filename' order by serial";
    my $stmt1 = $handle->prepare($qstr1);
    $stmt1->execute();
    while(my $hr1 = $stmt1->fetchrow_hashref()) {
      my $id = $hr1->{id};
      push(@retlist, $id);
    }
  }
  return(@retlist);
}
# }}}

# {{{ sub orgSelecProtIds hash(organism, locus_tags) returns a list of ids
sub orgSelectProtIds {
  my $self = shift(@_);
  my %args = @_;
  my $organism = $args{organism};
#  print(STDERR "$organism\n");
  my @lts = @{$args{locus_tags}};
  my @filenames = $self->organism_fn($organism);
  my @temp = map { "'" . $_ . "'"} @filenames;
  my $qfn = join(", ", @temp);
  my @retlist;
  foreach my $lt (@lts) {
    my $qstr1 = "select filename || '_' || serial as id from features where filename in ($qfn) and locus_tag = '$lt'";
#    print(STDERR "$qstr1\n");
    my $stmt1 = $handle->prepare($qstr1);
    $stmt1->execute();
    while(my $hr1 = $stmt1->fetchrow_hashref()) {
      my $id = $hr1->{id};
      push(@retlist, $id);
    }
  }
  return(@retlist);
}
# }}}

# {{{ id2identicals
sub id2identicals {
my $self = shift(@_);
my $id = shift(@_);
my $temp = $handle->selectrow_array("select identical from proteins where id = '$id'");
  my @identicals = split(/\,\s*/, $temp);
  return(@identicals);
}
# }}}

# {{{ sub protein (hash(filename, serial)) returns(hash(aa_seq, locus_tag, organism, product);
sub protein {
  my $self = shift(@_);
  my %args = @_;
  my $qstr1 = "select locus_tag, product from features where filename = '$args{filename}' and serial = $args{serial}";
  my ($locus_tag, $product) = $handle->selectrow_array($qstr1);
  my $qstr2 = "select aa_seq from proteins where filename = '$args{filename}' and serial = $args{serial}";
  my ($aa_seq) = $handle->selectrow_array($qstr2);
  my $qstr3 = "select organism from orgtab where filename = '$args{filename}'";
  my ($organism) = $handle->selectrow_array($qstr3);
  return(aa_seq => $aa_seq, locus_tag => $locus_tag, organism => $organism, product => $product);
}
# }}}

# {{{ sub lt2faa (locus_tag, filename) returns a list of ids.
sub lt2faa {
  my $self = shift(@_);
  my $lt = shift(@_);
  my $fn = shift(@_);
  my @idlist;
  my $qstr = "select filename, serial from features where locus_tag = '$lt'";
  my $stmt = $handle->prepare($qstr);
  $stmt->execute();
  while(my $hr = $stmt->fetchrow_hashref()) {
    my $id = $hr->{filename} . '_' . $hr->{serial};
    push(@idlist, $id);
  }
  if(scalar(@idlist) > 1) {
    print(STDERR "More than one identifier for $lt\n");
    return(undef);
  }
  else {
  my $id = shift(@idlist);
  open(OF, ">$fn");
  my $qstr1 = "select aa_seq from proteins where id = '$id'";
  my ($aa_seq) = $handle->selectrow_array($qstr1);
  $aa_seq =~ s/\*$//;
  #my $aa_seq;
  #print(STDERR "--- $filename\t$serial\n$aa_seq\n");
  print(OF ">$lt\n$aa_seq\n");
  close(OF);
  }
  return();
}
# }}}

# {{{ sub id2faa (id, filename) returns();
sub id2faa {
  my $self = shift(@_);
  my $id = shift(@_);
  my $fn = shift(@_);
  open(OF, ">$fn");
  my $qstr1 = "select aa_seq from proteins where id = '$id'";
  my ($aa_seq) = $handle->selectrow_array($qstr1);
  #my $aa_seq;
  #print(STDERR "--- $filename\t$serial\n$aa_seq\n");
  print(OF ">$id\n$aa_seq\n");
  close(OF);
  return();
}
# }}}

# {{{ sub id2aa_seq (id) returns(string);
sub id2aa_seq {
  my $self = shift(@_);
  my $id = shift(@_);
  my $qstr1 = "select aa_seq from proteins where id = '$id'";
  my ($aa_seq) = $handle->selectrow_array($qstr1);
  return($aa_seq);
}
# }}}

# {{{ sub id2nt_seq (id) returns(string);
sub id2nt_seq {
  my $self = shift(@_);
  my $id = shift(@_);
  my $qstr1 = "select nt_seq from proteins where id = '$id'";
  my ($nt_seq) = $handle->selectrow_array($qstr1);
  return($nt_seq);
}
# }}}

# {{{ id2lt (protein_id) return(locus_tag)
sub id2lt {
my $self = shift(@_);
my $id = shift(@_);
my %ret = $self->id2fileser(id => $id);
my $filename = $ret{filename};
my $serial = $ret{serial};
my $qstr = "select locus_tag from features where filename = '$filename' and serial = $serial";
my ($lt) = $handle->selectrow_array($qstr);
return($lt);
}
# }}}

# {{{ pre2org (pre) return(organism_name)
sub pre2org {
my $self = shift(@_);
my $pre = shift(@_);
my $qstr = "select organism from orgtab where pre = '$pre'";
my ($org) = $handle->selectrow_array($qstr);
return($org);
}
# }}}

# {{{ id2gbkFile (protein_id) return(a hash with locus_tag, description, organism, tta, porg);
sub id2gbkFile {
my $self = shift(@_);
my $id = shift(@_);
my $qstr = "select filename from proteins where id = '$id'";
my ($filename) = $handle->selectrow_array($qstr);
$qstr = "select dorf from orgtab where filename = '$filename'";
my ($dorf) = $handle->selectrow_array($qstr);

my $dir;
if($dorf=~m/^d/i) {
$dir = $draft_dir;
}
else {
$dir = $finished_dir;
}
my $fullPath = $dir . '/' . $filename;
return($fullPath);
}
# }}}

# {{{ id2detail (protein_id) return(a hash with organism, locus_tag, aa_seq, product, dorf, porg);
sub id2detail {
my $self = shift(@_);
my $id = shift(@_);
  my $qstr1 = "select locus_tag, pdesc, aa_seq, organism, porg, filename from vprot where id = '$id'";
  my ($locus_tag, $product, $aa_seq, $organism, $porg, $filename) = $handle->selectrow_array($qstr1);
  return(aa_seq => $aa_seq, locus_tag => $locus_tag,
  organism => $organism, product => $product, porg => $porg,
  filename => $filename 
  );
}
# }}}

# {{{ organism_fn (organism) return(list of filenames)
sub organism_fn {
  my $self = shift(@_);
  my $organism = shift(@_);
  my $qorg = $handle->quote($organism);
  my $qstr = "select filename from orgtab where organism = $qorg order by nt_len DESC";
  my $stmt = $handle->prepare($qstr);
  $stmt->execute();
  my @filenames;
  while(my $hr = $stmt->fetchrow_hashref()) {
    push(@filenames, $hr->{filename});
  }
  $stmt->finish();
  return(@filenames);
}
# }}}

# {{{ sub organism_faa (organism, outfile) faa gets written to outfile.
sub organism_faa {
  my $self = shift(@_);
  my $organism = shift(@_);
  my $outfile = shift(@_);
  my @filenames = $self->organism_fn($organism);
### print(">>>>", join("\t", @filenames), "\n"); ### debugging
  my $seqout=Bio::SeqIO->new(-file => ">$outfile", -format => 'fasta');
  foreach my $filename (@filenames) {
    my $qstr1 = "select id, aa_seq from proteins where filename = '$filename' order by ordvec";
    my $stmt1 = $handle->prepare($qstr1);
    $stmt1->execute();
    while(my $hr1 = $stmt1->fetchrow_hashref()) {
      my $id = $hr1->{id};
      my $aa_seq = $hr1->{aa_seq};
      my $outobj = Bio::Seq->new(-seq => $aa_seq);
      $outobj->display_name($id);
      $seqout->write_seq($outobj);
    }
  }
}
# }}}

# {{{ sub makeBlastDb (organism, directory) Makes a faa file and a blast database for the organism in directory. 
# This is not the one used by bigortho.
sub makeBlastDb {
  my $self = shift(@_);
  my $organism = shift(@_);
  my $bldbdir = shift(@_);
  if($bldbdir) {
  $bldbdir=~s/\/+$//; # get rid of any trailing slashes
  $bldbdir .= '/'; # and add just one slash 
  }
  my ($temp) = $handle->selectrow_array("select blastdb from orgbldb where organism = '$organism'");
  my $bldbname = $bldbdir . $temp;
  my $faa_name = $bldbdir . $temp . '.faa';
  $self->prog("$organism\t$bldbname\t$faa_name", "n");
  $self->organism_faa($organism, $faa_name);
  my @syslist = ("/usr/local/blastplus/bin/makeblastdb", '-in', $faa_name,
  '-title', $organism, '-out', $bldbname, '-dbtype' , 'prot');
  system(@syslist);
}

# }}}

# {{{ sub query2blastpdb (hash(query, title, out, exact)) Makes a blast database for the query organism(s). 
sub query2blastpdb {
  my $self = shift(@_);
  my %args = @_;
  my $qstr0;
  if($args{exact}) {
  $qstr0 = qq/select filename from orgtab where organism = '$args{query}'/;
  }
  else {
  $qstr0 = qq/select filename from orgtab where organism ~* '$args{query}'/;
  }
  my $stmt0 = $handle->prepare($qstr0);
  $stmt0->execute();
  my @fns;
  while(my $hr = $stmt0->fetchrow_hashref()) {
    push(@fns, $hr->{filename});
  }
  my @temp = map { "'" . $_ . "'"} @fns;
  my $fnstr = join(" ", @temp);
  my $qstr = qq/select id, aa_seq from proteins where filename in ($fnstr)/;
  my $stmt = $handle->prepare($qstr);
  $stmt->execute();
  my($fh, $fn)=tempfile($template, DIR => $tempdir, SUFFIX => '.faa');
  while(my $hr = $stmt->fetchrow_hashref()) {
    my $id = $hr->{id};
    print($fh ">$id\n$hr->{aa_seq}\n");
  }
  close($fh);
  qx(/usr/local/blastplus/bin/makeblastdb -in $fn -title "$args{title}" -out $args{out} -dbtype prot);
#  copy($fn, "lastfaa.faa");
  unlink($fn);
  return();
}

# }}}

# {{{ sub organism2blastpdb (hash(organism, outdir)). Makes a blast database for the query organism. 
# This is the one used by bigortho.
sub organism2blastpdb {
  my $self = shift(@_);
  my %args = @_;
  my $qorg = $handle->quote($args{organism});
  my $checkq = qq/select sum(n_cds) from orgtab where organism = $qorg/;
  my ($org_n_cds) = $handle->selectrow_array($checkq);
  if($org_n_cds == 0) {
    carp "$args{organism} has no CDSs $org_n_cds";
    return();
  }
  my $qstr0 = qq/select shaorg from orgtab where organism = $qorg/;
  my ($shaorg) = $handle->selectrow_array($qstr0);
  my $dbout = $args{outdir} . '/' . $shaorg;
  my $bltitle = $args{organism};
  $bltitle=~s/\'|\"//g;
  
  my $qstr1 = qq/select filename from orgtab where organism = $qorg/;
  my $stmt1 = $handle->prepare($qstr1);
  $stmt1->execute();
  my @fns;
  while(my $hr = $stmt1->fetchrow_hashref()) {
    push(@fns, $hr->{filename});
  }
  if(scalar(@fns) > 10) {
    print(STDERR "More than 10 files for $qorg. Skipping\n");
    return(0);
  }
  my @temp = map { "'" . $_ . "'"} @fns;
  my $fnstr = join(", ", @temp);

  my $tstr = qq/select count(*) from proteins where filename in ($fnstr) and id is NULL/;
  my ($noIdCount) = $handle->selectrow_array($tstr);
  if($noIdCount > 0) {
    print(STDERR "Some proteins without any id.\n");
    print(STDERR "Not making any BLAST database for $fnstr.\n");
    return(0);
  }

  my $qstr = qq/select id, aa_seq from proteins where filename in ($fnstr)/;
  my $stmt = $handle->prepare($qstr);
  $stmt->execute();
  my($fh, $fn)=tempfile($template, DIR => $tempdir, SUFFIX => '.faa');
  while(my $hr = $stmt->fetchrow_hashref()) {
    my $id = $hr->{id};
    print($fh ">$id $hr->{product}\n$hr->{aa_seq}\n");
  }
  close($fh);
  qx(/usr/local/blastplus/bin/makeblastdb -in $fn -title "$bltitle" -out $dbout -dbtype prot);
#  copy($fn, "lastfaa.faa");
  unlink($fn);
  return(1);
}

# }}}

# {{{ sub get_handle
sub get_handle {
my $self = shift(@_);
return($handle);
}
# }}}

# {{{ sub get_newhandle # you have to disconnect this one in the calling routine.
sub get_newhandle {
my $self = shift(@_);
my $newhandle=DBI->connect("DBI:Pg:dbname=$dbname;host=$dbhost", 'sco', 'tsinH4x');
return($newhandle);
}
# }}}

# {{{ sub get_gbdir
sub get_gbdir {
my $self = shift(@_);
my $type = shift(@_);
if($type=~m/^f/i) { return($finished_dir); }
elsif($type=~m/^d/i) { return($draft_dir); }
else { croak "Provide a type (f or d) to get the gbdir\n"; }
}
# }}}

# {{{ taxonomy (accession) return (organism, taxonomy, porg);
sub taxonomy {
my $self=shift(@_);
my $accession = shift(@_);

my $qstr="select organism, taxonomy, definition from orgtab where accession like '\%$accession%'";
my ($organism, $taxo, $definition)=$handle->selectrow_array($qstr);
my $porg;
if($definition=~m/plasmid/) {
  $porg = 'P'
}
else {$porg = 'G';}


return($organism, $taxo, $porg);
}
# }}}

# {{{ sub ltord
sub ltord {
  my $self = shift(@_);
  my $accession = shift(@_);
  my $gbdir = '/home/sco/seq/nt/genbank_ftp';
  my %order_hash;
  my $id = shift(@_);
  my $gbfile = $gbdir .'/'. $accession . '.gbk';
#  print(STDERR "$gbfile\n");
# sometimes there is no genbank file for an accession.
  unless(-e $gbfile) {print(STDERR "$gbfile does not exist\n"); return();}
  my $seqio=Bio::SeqIO->new(-file => $gbfile);
  my $seqobj=$seqio->next_seq();
  my $count=0;
  my @unordered = $seqobj->all_SeqFeatures();
  my @ordered = sort feature_sorter @unordered;
  foreach my $feature (@ordered) {
    if($feature->primary_tag() eq 'CDS') {
      my $locus_tag;
      if($feature->has_tag('locus_tag')) {
      $locus_tag = join(" ", $feature->get_tag_values('locus_tag'));
      }
      elsif($feature->has_tag('systematic_id')) {
      my @temp = $feature->get_tag_values('systematic_id');
      $locus_tag = $temp[0];
      }
      else {$locus_tag = $accession . '_' . $count;}
      $order_hash{$locus_tag} = $count;
      $count+=1;
    }
  }
  return(\%order_hash);
}
# }}}

# {{{ internal sub feature_sorter used by ltord
sub feature_sorter {
return($a->start() <=> $b->start());
}
# }}}

# {{{ sub get_anno
sub get_anno {
  my $self = shift(@_);
  my $accession = shift(@_);
  my $wanted_lt = shift(@_);
  my $gbdir = '/home/sco/seq/nt/genbank_ftp';
  my $gbfile = $gbdir .'/'. $accession . '.gbk';
#  print(STDERR "$gbfile\n");
# sometimes there is no genbank file for an accession.
  unless(-e $gbfile) {print(STDERR "$gbfile does not exist\n"); return();}
  my $seqio=Bio::SeqIO->new(-file => $gbfile);
  my $seqobj=$seqio->next_seq();
  my $count=0;
  my @unordered = $seqobj->all_SeqFeatures();
      my $product;
  foreach my $feature (@unordered) {
    if($feature->primary_tag() eq 'CDS') {
      my $locus_tag;
      if($feature->has_tag('product')) {
      $product = join(" ", $feature->get_tag_values('product'));
      }
      if($feature->has_tag('locus_tag')) {
      $locus_tag = join(" ", $feature->get_tag_values('locus_tag'));
      }
      elsif($feature->has_tag('systematic_id')) {
      my @temp = $feature->get_tag_values('systematic_id');
      $locus_tag = $temp[0];
      }
      else {$locus_tag = $accession . '_' . $count;}
      $count+=1;
      if($locus_tag eq $wanted_lt) {
        return($product);
      }
    }
  }
  return();
}
# }}}

# {{{ sub DESTROY
sub DESTROY {
  my $self = shift(@_);
  $handle->disconnect();
}
# }}}

# {{{ sub get_orgdet
sub get_orgdet {
my $self = shift(@_);
my $filename = shift(@_);
my $type = shift(@_);
my $qstr = "select accession, locus, filename, definition, organism, taxonomy, description, nt_len,
gcratio, n_cds, n_rrna, n_trna, n_psgenes, n_pscds, linear
from orgtab where filename = ?";
my $stmt = $handle->prepare($qstr);
$stmt->execute($filename);
my $hr = $stmt->fetchrow_hashref();
return($hr);
#my($accession, $locus, $filename, $definition, $organism, $taxonomy, $description, $nt_len,
#$gcratio, $n_cds, $n_rrna, $n_trna, $n_psgenes, $n_pscds, $linear);
#$stmt->bind_columns(\$accession, \$locus, \$filename, \$definition, \$organism, \$taxonomy, \$description, \$nt_len,
#\$gcratio, \$n_cds, \$n_rrna, \$n_trna, \$n_psgenes, \$n_pscds, \$linear);
#$stmt->fetch();
}
# }}}

# {{{ sub proteome (organism, type) populates proteins from genbank file.
sub proteome {
  my $self = shift(@_);
  my $organism = shift(@_);
  my $type = shift(@_);
  my $gbdir = $self->get_gbdir($type);
  my @filepre;
  my $qstr = "select filename, pre from orgtab where organism = '$organism'";
  my $stmt = $handle->prepare($qstr);
  $stmt->execute();
  while(my $hr=$stmt->fetchrow_hashref()) {
    push(@filepre, [$hr->{filename}, $hr->{pre}]);
  }
  $stmt->finish();
  foreach my $fpr (@filepre) {
    my $fullpath=$gbdir . '/' . $fpr->[0];
    my $pre = $fpr->[1];
    my $seqio=Bio::SeqIO->new(-file => $fullpath);
    my $seq=$seqio->next_seq();
    my $accession = $seq->accession_number();
### get all features ####
    my @featurelist=($seq->all_SeqFeatures());
    my @sorted_featurelist=sort _feature_sorter @featurelist;
    my $count=1;
    foreach my $feature (@sorted_featurelist) {
      if($feature->primary_tag() eq 'CDS') {
	my $locus_tag;
	if($feature->has_tag('locus_tag')) {
	  $locus_tag = join(" ", $feature->get_tag_values('locus_tag'));
	}
	elsif($feature->has_tag('systematic_id')) {
	  my @temp = $feature->get_tag_values('systematic_id');
	  $locus_tag = $temp[0]; 
	}
	else {$locus_tag = $accession . '_' . $count;}
	my $product;
	if($feature->has_tag('product')) {
	  $product=join(" ", $feature->get_tag_values('product'));
	}
	my $gene;
	if($feature->has_tag('gene')) {
	  $gene=join(" ", $feature->get_tag_values('gene'));
	}
	my $feat_aaobj=&_feat_trans($feature);
	my $id = $pre . '_' . $count;
	$feat_aaobj->display_name($id);
	my $aa_seq = $feat_aaobj->seq();
	my $description = join(" ", ($gene, $product));
	my $qdesc = $handle->quote($description);
	$feat_aaobj->description("$gene $product");
#$seqout->write_seq($feat_aaobj);
	$count+=1;
	my $instr = "insert into proteins (filename, id, locus_tag, aa_seq, description) values 
	  ('$fpr->[0]', '$id', '$locus_tag', '$aa_seq', $qdesc)"; 
	  if($handle->do($instr)) {
	    print("Inserted $fpr->[0]\t$id\n");
	  }
	  else {
	    print(STDERR "$instr\n");
	  }
      }
    }
  }
}
# }}}

# {{{ sub _feature_sorter {
sub _feature_sorter {
return($a->start() <=> $b->start());
}
# }}}

# {{{ sub _feat_trans

=head2 sub _feat_trans

Given a Bio::SeqFeature::Generic returns a Bio::Seq for the
translation of the spliced seq


=cut

sub _feat_trans {
  my $feature=shift(@_);
  my $codon_start=1;
  if($feature->has_tag('codon_start')) {
      ($codon_start) = $feature->get_tag_values('codon_start');
      }
      my $offset=1;
      if($codon_start > 1) { $offset = $codon_start;}
      my $featobj=$feature->spliced_seq(-nosort => '1');
      my $aaobj=$featobj->translate(-offset => $offset);
  return($aaobj);
}
# for translate() see Bio::PrimarySeqI documentation.
# }}}

# {{{ sub blastp (hash(query, database)) returns (filename) # unlink returned filename;
sub blastp {
my $self = shift(@_);
my %args = @_;
  my($fh, $fn)=tempfile($template, DIR => $tempdir, SUFFIX => '.blast');
  close($fh);
  `/usr/local/blastplus/bin/blastp -query $args{query} -db $args{database} -out $fn -evalue $args{ethresh}`;
return($fn);
}
# }}}


# {{{ sub hit2query (hash(id)) returns (queryname)
sub hit2query {
my $self = shift(@_);
my %args = @_;
my $qstr = "select query from recihits where hit = '$args{id}'";
my ($hit) = $handle->selectrow_array($qstr);
return($hit);
}
# }}}


# {{{ sub genbank2seqstr_nobp (genbankFileName) returns(nucleotide sequence string)
sub genbank2seqstr_nobp {
my $self = shift(@_);
my $filename = shift(@_);
open(IN, "<$filename");
my $seqflag = 0;
my $seq;
while(<IN>) {
my $line = $_;
chomp($line);
if($line=~m/^ORIGIN/) {
$seqflag = 1; next;
}
if($line=~m/^\/\//) { $seqflag = 0; last; }
if($seqflag) {
$line=~s/\d+|\s+//g;
$seq .= $line;
}
}
return($seq);
}
# }}}

# {{{ sub id2fileser (hash(id)) returns(hash(filename, serial));
sub id2fileser {
my $self = shift(@_);
my %args = @_;
my $id = $args{id};
my ($filename, $serial) = split(/_(?!.*_.*$)/, $id);
return(filename => $filename, serial => $serial);
}
# }}}


# {{{ sub upstream_n (hash(filename, serial, n)) returns(string of nt sequence) 
sub upstream_n {
  my $self = shift(@_);
  my %args = @_;
  my $filename = $args{filename};
  my $serial = $args{serial};
  my $qstr = qq/select start_pos, end_pos, strand from features where filename = '$filename' and serial = $serial/;
  my($start_pos, $end_pos, $strand) = $handle->selectrow_array($qstr);
  my $nt_seq;
  if($strand == 1) {
  my $sectionEnd = $start_pos - 1;
  my $sectionStart = $sectionEnd - ($args{n} - 1) >= 1 ? $sectionEnd - ($args{n} -1) : 1;
  $nt_seq = nt_subseq($self, filename => $filename,
                      end => $sectionEnd,
                      start => $sectionStart,
                      revcom => 0);
  }
  if($strand == -1) {
  my $sectionStart = $end_pos + 1;
  my $sectionEnd = $sectionStart + ($args{n} -1);
  $nt_seq = nt_subseq($self, filename => $filename,
                     start => $sectionStart,
                     end => $sectionEnd,
                     revcom => 1);
  }
  return($nt_seq);
}
# }}}

# {{{ sub upstream_nm (hash(filename, serial, n, m)) returns(string of nt sequence)
# upstream n and then m nucleotides into the gene
sub upstream_nm {
  my $self = shift(@_);
  my %args = @_;
  my $filename = $args{filename};
  my $serial = $args{serial};
  my $qstr = qq/select start_pos, end_pos, strand from features where filename = '$filename' and serial = $serial/;
  my($start_pos, $end_pos, $strand) = $handle->selectrow_array($qstr);
  my $nt_seq;
  if($strand == 1) {
  my $sectionEnd = $start_pos + $args{m};
  my $sectionStart = $start_pos - ($args{n} - 1) >= 1 ? $start_pos - ($args{n} -1) : 1;
  $nt_seq = nt_subseq($self, filename => $filename,
                      end => $sectionEnd,
                      start => $sectionStart,
                      revcom => 0);
#  print(STDERR "$sectionStart\t$sectionEnd\n");
  }
  if($strand == -1) {
  my $sectionStart = $end_pos - $args{m};
  my $sectionEnd = $end_pos + ($args{n} -1);
  $nt_seq = nt_subseq($self, filename => $filename,
                     start => $sectionStart,
                     end => $sectionEnd,
                     revcom => 1);
#  print(STDERR "$sectionStart\t$sectionEnd\n");
  }
  return($nt_seq);
}
# }}}


# {{{ sub upstream (hash(filename, serial)) returns(string of nt sequence) 
sub upstream {
  my $self = shift(@_);
  my %args = @_;
  my $filename = $args{filename};
  my $serial = $args{serial};
  my $qstr = qq/select start_pos, end_pos, strand from features where filename = '$filename' and serial = $serial/;
  my($start_pos, $end_pos, $strand) = $handle->selectrow_array($qstr);
  my $nt_seq;
  if($strand == 1) {
  $qstr = qq/select max(end_pos) from features where end_pos < $start_pos and filename = '$filename'/;
  my ($up_end) = $handle->selectrow_array($qstr);
  $nt_seq = nt_subseq($self, filename => $filename,
                      start => $up_end + 1,
                      end => $start_pos - 1,
                      revcom => 0);
  }
  if($strand == -1) {
  $qstr = qq/select min(start_pos) from features where start_pos > $end_pos and filename = '$filename'/;
  my ($up_start) = $handle->selectrow_array($qstr);
  $nt_seq = nt_subseq($self, filename => $filename,
                     start => $end_pos + 1,
                     end => $up_start - 1,
                     revcom => 1);
  }
  return($nt_seq);
}
# }}}

# {{{ sub nt_subseq (hash(filename, start, end, revcom)) returns(nt seq string)
sub nt_subseq {
  my $self = shift(@_);
  my %args = @_;
  if($args{start} >= $args{end}) { return(); }
#  print(STDERR "$args{start}\t$args{end}\n");
  my $filename = $args{filename};
  my $wantlen = ($args{end} - $args{start}) + 1;
  my $qstr = qq/select substring(nt_seq from $args{start} for $wantlen) as nt_seq from contigs where filename = '$filename'/;
  my $nt_seq = $handle->selectrow_array($qstr);
  if($args{revcom}) {
    my $seqobj = Bio::Seq->new(-seq => $nt_seq);
    my $revobj = $seqobj->revcom();
    my $revseq = $revobj->seq();
    return($revseq);
  }
  else {
  return($nt_seq);
  }
}
# }}} 


return(1);

__END__

perl ../emblcds2faa.pl -file Afu_chr1.genbank -fo genbank -pre afu > afu.faa
perl ../emblcds2faa.pl -file Afu_chr2.genbank -fo genbank -pre afu >> afu.faa
perl ../emblcds2faa.pl -file Afu_chr3.genbank -fo genbank -pre afu >> afu.faa
perl ../emblcds2faa.pl -file Afu_chr4.genbank -fo genbank -pre afu >> afu.faa
perl ../emblcds2faa.pl -file Afu_chr5.genbank -fo genbank -pre afu >> afu.faa
perl ../emblcds2faa.pl -file Afu_chr6.genbank -fo genbank -pre afu >> afu.faa
perl ../emblcds2faa.pl -file Afu_chr7.genbank -fo genbank -pre afu >> afu.faa
perl ../emblcds2faa.pl -file Afu_chr8.genbank -fo genbank -pre afu >> afu.faa



}

}


1;

__END__




