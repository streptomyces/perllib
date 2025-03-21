package Sco::Igbk;
use 5.14.0;
use Bio::SeqIO;
use Bio::Seq;
use Carp;
use DBI;
use File::Basename;
our $AUTOLOAD;
use lib qw(/home/sco /home/sco/perllib);

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

# {{{ sub fna %(seqobj, ofn) returns true.
sub fna {
  my $self = shift(@_);
  my %args = @_;
  my $filename;
  if($args{localname}) {
    $filename = $args{file};
  }
  else {
  $filename = $gbkDir . '/' . $args{file};
  }
  if(not -s $filename) {
    linelistE("$filename does not exist or is zero in size");
    return(0);
  }
  my $seqio = Bio::SeqIO->new(-file => $filename);
  my $seqout = Bio::SeqIO->new(-file => ">$args{outfilename}", -format => 'fasta');
  my $seqobj = $seqio->next_seq();
  $seqobj->description(undef);
  if($args{name}) {
    $seqobj->display_name($args{name});
  }
  $seqout->write_seq($seqobj);
  return(1);
}
# }}}

# {{{ sub gbkVitals %(file or seqobj) returns a hash

sub gbkVitals {
  my $self = shift(@_);
  my %args = @_;
  my $dorf;
  my ($seqio, $seqobj, $gbkfile);
  if($args{file}) {
  $gbkfile = $args{file};
  $seqio=Bio::SeqIO->new(-file => $gbkfile);
  $seqobj=$seqio->next_seq();
  }
  elsif($args{seqobj}) {
    $seqobj = $args{seqobj};
  }

  my $acc=$seqobj->accession();
  my ($species, @cl, $taxonomy, $binomial);
  $species=$seqobj->species();
  if($species) {
    @cl = $species->classification();
    $taxonomy = join(" : ", reverse(@cl));
    $binomial =$species->binomial('FULL');
  }
  else {
    $taxonomy = undef;
    $binomial = undef;
  }
  my $nt_len=$seqobj->length();
  my $n_cds = 0; 
  my $n_rrna = 0; 
  my $n_trna = 0;
  my $n_psgene = 0;
  my $n_pscds = 0;
  foreach my $feature ($seqobj->all_SeqFeatures()) {
    if($feature->primary_tag() eq 'CDS') {
      $n_cds+=1;
      if($feature->has_tag('pseudo')){
        $n_pscds+=1;
      }
    }
    elsif($feature->primary_tag() eq 'gene') {
      if($feature->has_tag('pseudo')){
        $n_psgene+=1;
      }
    }
    elsif($feature->primary_tag() eq 'tRNA') {
      $n_trna+=1;
    }
    elsif($feature->primary_tag() eq 'rRNA') {
      $n_rrna+=1;
    }
  }
  my $desc=$seqobj->description();
  my $porg;
  if($desc=~m/plasmid/i) { $porg = 'p'; }  else { $porg = 'g'; }
  my (undef, undef, $gcratio) = Sco::Global->gc_frac($seqobj->seq());
  my $linear;
  if($seqobj->is_circular()) {
    $linear = 0;
  }
  else {
    $linear = 1;
  }
  $desc=~s/\'/\'\'/g;
  $taxonomy=~s/\'/\'\'/g;
  $binomial=~s/\'/\'\'/g;
  my %rethash = (accession => $acc, filename => $gbkfile, porg => $porg,
  description => $desc, binomial => $binomial, organism => $binomial,
  taxonomy => $taxonomy, nt_len => $nt_len, gc_ratio => $gcratio, n_cds => $n_cds,
  n_rrna => $n_rrna, n_trna => $n_trna, n_psgene => $n_psgene, n_pscds => $n_pscds,
  linear => $linear, definition => $desc);
#  my $instr = qq/insert into orgtab (accession, filename, definition, organism, taxonomy, nt_len, gcratio, n_cds, n_rrna, n_trna, n_psgenes, n_pscds, linear, dorf, porg) values ('$acc', '$filext', '$desc', '$binomial', '$taxonomy', $size, $gcratio, $n_cds, $n_rrna, $n_trna, $n_psgene, $n_pscds, $linear, '$dorf', '$porg')/;
  return(%rethash);
}

# }}}

# {{{ sub bedfile (hash(accession|file, molecule, outfilename,
# optional:tagasid, features[]))
# Writes bed file for the molecule to the outfilename.
sub bedfile {
my $self = shift(@_);
my %args = @_;
my $molecule = $args{molecule};
my $gbkfile;

if($args{accession}) {
$gbkfile = $gbkDir . '/' . $args{accession} . '.gbk';
}
elsif($args{file}) {
$gbkfile = $args{file};
}

my $tagasid = qw(locus_tag);
if(exists($args{tagasid})) {
$tagasid = $args{tagasid};
}

my %featcol = (
rRNA => q/255,255,0/,
tRNA => q/200,200,200/,
misc_RNA => q/0,128,128/
);

my @do_feats = qw(tRNA rRNA CDS misc_RNA nc_RNA);
if($args{features}) {
push(@do_feats, @{$args{features}});
}
my $score=0;
my $track = qw(features);

my $seqio=Bio::SeqIO->new(-file => $gbkfile);

my $ofh;
open($ofh, ">$args{outfilename}");


print $ofh <<"TRACK";
track\tname=$track\tdescription="Features"\tuseScore=0\tvisibility=1\titemRgb=on\tcolor=100,20,20
TRACK

# print("track\tname=Features\tdescription="Features"\tuseScore=0\tvisibility=1\tcolor=255,255,255");
while(my $seqobj=$seqio->next_seq()) {
  foreach my $feature ($seqobj->all_SeqFeatures()) {
    my $pritag = $feature->primary_tag();

    if(grep {$_ eq $pritag} @do_feats) {
      my $start = $feature->start();
      my $end = $feature->end();
      my $strand = $feature->strand() == 1 ? '+' : '-';
      my $itemrgb = $feature->strand() == 1 ? '0,255,0' : '0,0,255';
      my $annotation;
      my $id;
      my @tags = $feature->get_all_tags();
      foreach my $tag (@tags) {
        if($tag=~m/product|gene|note/i) {
          my $annostr=join(" ", $feature->get_tag_values($tag));
          $annotation.=$annostr;
        }
        if($tag=~m/^$tagasid/) {
          my $lt=join("|", $feature->get_tag_values($tag));
          $id .= $lt;
        }
      }
      #if($pritag=~m/rRNA/) {$itemrgb = '255,255,0';}
      #if($pritag=~m/tRNA/) {$itemrgb = '200,200,200';}
      unless($pritag eq 'CDS') {
        $itemrgb = $featcol{$pritag};
      }
      unless($id) {$id = $pritag;}
      my $bedstart = $start - 1;
      print($ofh "$molecule\t$bedstart\t$end\t$id\t$score\t$strand\t$bedstart\t$end\t$itemrgb\n");
    }
  }
#$emblout->write_seq($seqobj);
}
close($ofh);
}
# }}}

# {{{ sub seqobjFeatPrexSQLite %(seqobj, sqlitefn, organism) returns true.
sub seqobjFeatPrexSQLite {
  my $self = shift(@_);
  my %args = @_;
  my $seqobj = $args{seqobj};
  my $orgname = $args{organism};
  my @temp = $seqobj->all_SeqFeatures();
  my $accession = $seqobj->accession();
  my @features = sort _feat_sorter(@temp);

  my($fh, $fn);
  if($args{filename}) { $fn = $args{filename}; }
  else {
  ($fh, $fn)=tempfile($template, DIR => $tempdir, SUFFIX => '.sqlite');
  close($fh);
  }
  
  my $handlite=DBI->connect("DBI:SQLite:dbname=$args{sqlitefn}", '', '');
  my $qacc = $handlite->quote($accession);

# create table features (organism text, accession text, locus_tag text,
# pritag text, start_pos integer, end_pos integer, strand integer,
# product text, olt text,
# unique (accession, locus_tag)
# );

  my @retlist;
  my $featCnt = 0;
  my $insertCnt = 0;
  foreach my $feat (@features) {
    my $pritag = $feat->primary_tag();
    my $start_pos = $feat->start();
    my $end_pos = $feat->end();
    my $strand = $feat->strand();
    if($pritag eq 'CDS' or $pritag=~m/[tr]RNA/) {
      $featCnt += 1;
      my $locus_tag;
      if($feat->has_tag('locus_tag')) {
        my @lts = $feat->get_tag_values('locus_tag');
        $locus_tag = $lts[0];
      }
      my $old_locus_tag;
      if($feat->has_tag('old_locus_tag')) {
        my @lts = $feat->get_tag_values('old_locus_tag');
        $old_locus_tag = $lts[0];
      }
      my $product;
      if($feat->has_tag('product')) {
        my @prods = $feat->get_tag_values('product');
        $product = join(" ", @prods);
      }
      my $dbxref;
      if($feat->has_tag('db_xref')) {
        my @temp = $feat->get_tag_values('db_xref');
        $dbxref = join(" ", @temp);
      }
      my $proteinid;
      if($feat->has_tag('protein_id')) {
        my @temp = $feat->get_tag_values('protein_id');
        $proteinid = $temp[0];
      }
      my $id = defined($locus_tag) ? $locus_tag : $pritag . "_" . $featCnt;
      unless($old_locus_tag) { $old_locus_tag = "NULL"; }
      #my %featrec = (id => $id, pritag => $pritag, start => $start_pos,
      #end => $end_pos, strand => $strand, product => $product);
      #if($old_locus_tag) { $featrec{olt} = $old_locus_tag; }
      #if($proteinid) { $featrec{proteinid} = $proteinid; }
      #if($dbxref) { $featrec{dbxref} = $dbxref; }
      my $qid = $handlite->quote($id);
      my $qptag = $handlite->quote($pritag);
      my $qprod = $handlite->quote($product);
      my $qorg = $handlite->quote($orgname);
      my $qolt;
      if($old_locus_tag eq 'NULL') { $qolt = $old_locus_tag; }
      else {
      $qolt = $handlite->quote($old_locus_tag);
      }
      my $instr = qq/insert into features (organism, accession, locus_tag, pritag, start_pos, end_pos, strand, product, olt)
      values ($qorg, $qacc, $qid, $qptag, $start_pos, $end_pos, $strand, $qprod, $qolt)/;
      if($handlite->do($instr)) {
        $insertCnt += 1;
      }
      else {
        print(STDERR "$instr\n");
        return(0);
        }
        
#      push(@retlist, {%featrec});
    }
  }
  return($insertCnt);
}
# }}}

# {{{ sub seqObj2featuresSQLite %(seqobj) returns SQLite DBIhandle, SQLite filename
sub seqObj2featuresSQLite {
  my $self = shift(@_);
  my %args = @_;
  my $seqobj = $args{seqobj};
  my @temp = $seqobj->all_SeqFeatures();
  my @features = sort _feat_sorter(@temp);

  my($fh, $fn);
  if($args{filename}) { $fn = $args{filename}; }
  else {
  ($fh, $fn)=tempfile($template, DIR => $tempdir, SUFFIX => '.sqlite');
  close($fh);
  }
  
  my $handlite=DBI->connect("DBI:SQLite:dbname=$fn", '', '');
  $handlite->do(qq/create table features (id text, pritag text, start_pos integer,
  end_pos integer, strand integer, product text, olt text)/);

  my @retlist;
  my $featCnt = 0;
  foreach my $feat (@features) {
    my $pritag = $feat->primary_tag();
    my $start_pos = $feat->start();
    my $end_pos = $feat->end();
    my $strand = $feat->strand();
    if($pritag eq 'CDS' or $pritag=~m/[tr]RNA/) {
      $featCnt += 1;
      my $locus_tag;
      if($feat->has_tag('locus_tag')) {
        my @lts = $feat->get_tag_values('locus_tag');
        $locus_tag = $lts[0];
      }
      my $old_locus_tag;
      if($feat->has_tag('old_locus_tag')) {
        my @lts = $feat->get_tag_values('old_locus_tag');
        $old_locus_tag = $lts[0];
      }
      my $product;
      if($feat->has_tag('product')) {
        my @prods = $feat->get_tag_values('product');
        $product = join(" ", @prods);
      }
      my $dbxref;
      if($feat->has_tag('db_xref')) {
        my @temp = $feat->get_tag_values('db_xref');
        $dbxref = join(" ", @temp);
      }
      my $proteinid;
      if($feat->has_tag('protein_id')) {
        my @temp = $feat->get_tag_values('protein_id');
        $proteinid = $temp[0];
      }
      my $id = defined($locus_tag) ? $locus_tag : $pritag . "_" . $featCnt;
      my %featrec = (id => $id, pritag => $pritag, start => $start_pos, end => $end_pos,
      strand => $strand, product => $product);
      if($old_locus_tag) { $featrec{olt} = $old_locus_tag; }
      if($proteinid) { $featrec{proteinid} = $proteinid; }
      if($dbxref) { $featrec{dbxref} = $dbxref; }
      my $qid = $handlite->quote($id);
      my $qptag = $handlite->quote($pritag);
      my $qprod = $handlite->quote($product);
      my $instr = qq/insert into features (id, pritag, start_pos, end_pos, strand, product)
      values ($qid, $qptag, $start_pos, $end_pos, $strand, $qprod)/;
      unless($handlite->do($instr)) {
        print(STDERR "$instr\n");
        return();
        }
#      push(@retlist, {%featrec});
    }
  }
  return($handlite, $fn);
}
# }}}

# {{{ sub organism2organism %(organism) returns a list of organism names
# This is to get the full name of an organism given its partial name.
sub organism2organism {
my $self = shift(@_);
my %args = @_;
my @retlist;
my $qstr = qq/select organism from orgtab where organism like '%$args{organism}%'/;
my $stmt = $handle->prepare($qstr);
$stmt->execute();
while(my $hr = $stmt->fetchrow_hashref()) {
unless(grep {$_ eq $hr->{organism}} @retlist) {
  push(@retlist, $hr->{organism});
}
}
return(@retlist);
}
# }}}

# {{{ sub seqObj2featuresTable %(seqobj) returns a list of hashrefs;
sub seqObj2featuresTable {
  my $self = shift(@_);
  my %args = @_;
  my $seqobj = $args{seqobj};
  my @temp = $seqobj->all_SeqFeatures();
  my @features = sort _feat_sorter(@temp);
  my @retlist;
  my $featCnt = 0;
  foreach my $feat (@features) {
    my $pritag = $feat->primary_tag();
    my $start_pos = $feat->start();
    my $end_pos = $feat->end();
    my $strand = $feat->strand();
    
    if($pritag eq 'CDS' or $pritag=~m/[tr]RNA/) {
      $featCnt += 1;
      my $locus_tag;
      if($feat->has_tag('locus_tag')) {
        my @lts = $feat->get_tag_values('locus_tag');
        $locus_tag = $lts[0];
      }
      my $old_locus_tag;
      if($feat->has_tag('old_locus_tag')) {
        my @lts = $feat->get_tag_values('old_locus_tag');
        $old_locus_tag = $lts[0];
      }
      my $product;
      if($feat->has_tag('product')) {
        my @prods = $feat->get_tag_values('product');
        $product = join(" ", @prods);
      }
      my $gene;
      if($feat->has_tag('gene')) {
        my @temp = $feat->get_tag_values('gene');
        $gene = join(" ", @temp);
      }
      my $note;
      if($feat->has_tag('note')) {
        my @temp = $feat->get_tag_values('note');
        $note = join(" ", @temp);
      }
      my $dbxref;
      if($feat->has_tag('db_xref')) {
        my @temp = $feat->get_tag_values('db_xref');
        $dbxref = join(" ", @temp);
      }
      my $proteinid;
      if($feat->has_tag('protein_id')) {
        my @temp = $feat->get_tag_values('protein_id');
        $proteinid = $temp[0];
      }

      my $id = defined($locus_tag) ? $locus_tag : $pritag . "_" . $featCnt;
      my %featrec = (id => $id, pritag => $pritag, start => $start_pos, end => $end_pos,
      strand => $strand, product => $product, note => $note, gene => $gene);
      if($old_locus_tag) { $featrec{olt} = $old_locus_tag; }
      if($proteinid) { $featrec{proteinid} = $proteinid; }
      if($dbxref) { $featrec{dbxref} = $dbxref; }
      push(@retlist, {%featrec});
    }
  }
  return(@retlist);
}
# }}}

# {{{ sub gbkfile2featuresTable %(file) returns a list of hashrefs;
sub gbkfile2featuresTable {
  my $self = shift(@_);
  my %args = @_;
  my $filename;
  if(-r $args{file}) { $filename = $args{file}; } 
  else { $filename = $gbkDir . '/' . $args{file}; }
  my $seqio = Bio::SeqIO->new(-file => $filename);
  my $seqobj = $seqio->next_seq();
  my @temp = $seqobj->all_SeqFeatures();
  my @features = sort _feat_sorter(@temp);
  my @retlist;
  my $featCnt = 0;
  foreach my $feat (@features) {
    my $pritag = $feat->primary_tag();
    my $start_pos = $feat->start();
    my $end_pos = $feat->end();
    my $strand = $feat->strand();
    
    if($pritag eq 'CDS' or $pritag=~m/[tr]RNA/) {
      $featCnt += 1;
      my $locus_tag;
      if($feat->has_tag('locus_tag')) {
        my @lts = $feat->get_tag_values('locus_tag');
        $locus_tag = $lts[0];
      }
      my $old_locus_tag;
      if($feat->has_tag('old_locus_tag')) {
        my @lts = $feat->get_tag_values('old_locus_tag');
        $old_locus_tag = $lts[0];
      }
      my $product;
      if($feat->has_tag('product')) {
        my @prods = $feat->get_tag_values('product');
        $product = join(" ", @prods);
      }
      my $gene;
      if($feat->has_tag('gene')) {
        my @temp = $feat->get_tag_values('gene');
        $gene = join(" ", @temp);
      }
      my $note;
      if($feat->has_tag('note')) {
        my @temp = $feat->get_tag_values('note');
        $note = join(" ", @temp);
      }
      my $dbxref;
      if($feat->has_tag('db_xref')) {
        my @temp = $feat->get_tag_values('db_xref');
        $dbxref = join(" ", @temp);
      }
      my $proteinid;
      if($feat->has_tag('protein_id')) {
        my @temp = $feat->get_tag_values('protein_id');
        $proteinid = $temp[0];
      }

      my $id = defined($locus_tag) ? $locus_tag : $pritag . "_" . $featCnt;
      my %featrec = (id => $id, pritag => $pritag, start => $start_pos, end => $end_pos,
      strand => $strand, product => $product, gene => $gene, note => $note);
      if($old_locus_tag) { $featrec{olt} = $old_locus_tag; }
      if($locus_tag) { $featrec{lt} = $locus_tag; }
      if($proteinid) { $featrec{proteinid} = $proteinid; }
      if($dbxref) { $featrec{dbxref} = $dbxref; }
      push(@retlist, {%featrec});
    }
  }
  return(@retlist);
}
# }}}

# {{{ sub gbkfile2featuresHash %(file) returns a hash of hashrefs;
sub gbkfile2featuresHash {
  my $self = shift(@_);
  my %args = @_;
  my $filename;
  if(-r $args{file}) { $filename = $args{file}; } 
  else { $filename = $gbkDir . '/' . $args{file}; }
  my $seqio = Bio::SeqIO->new(-file => $filename);
  my $seqobj = $seqio->next_seq();
  my @temp = $seqobj->all_SeqFeatures();
  my @features = sort _feat_sorter(@temp);
  my %rethash;
  my $featCnt = 0;
  foreach my $feat (@features) {
    my $pritag = $feat->primary_tag();
    my $start_pos = $feat->start();
    my $end_pos = $feat->end();
    my $strand = $feat->strand();
    
    if($pritag eq 'CDS' or $pritag=~m/[tr]RNA/) {
      $featCnt += 1;
      my $locus_tag;
      if($feat->has_tag('locus_tag')) {
        my @lts = $feat->get_tag_values('locus_tag');
        $locus_tag = $lts[0];
      }
      my $old_locus_tag;
      if($feat->has_tag('old_locus_tag')) {
        my @lts = $feat->get_tag_values('old_locus_tag');
        $old_locus_tag = $lts[0];
      }
      my $product;
      if($feat->has_tag('product')) {
        my @prods = $feat->get_tag_values('product');
        $product = join(" ", @prods);
      }
      my $gene;
      if($feat->has_tag('gene')) {
        my @temp = $feat->get_tag_values('gene');
        $gene = join(" ", @temp);
      }
      my $note;
      if($feat->has_tag('note')) {
        my @temp = $feat->get_tag_values('note');
        $note = join(" ", @temp);
      }
      my $dbxref;
      if($feat->has_tag('db_xref')) {
        my @temp = $feat->get_tag_values('db_xref');
        $dbxref = join(" ", @temp);
      }
      my $proteinid;
      if($feat->has_tag('protein_id')) {
        my @temp = $feat->get_tag_values('protein_id');
        $proteinid = $temp[0];
      }

      my $id = defined($locus_tag) ? $locus_tag : $pritag . "_" . $featCnt;
      my %featrec = (id => $id, pritag => $pritag, start => $start_pos, end => $end_pos,
      strand => $strand, product => $product, gene => $gene, note => $note);
      if($old_locus_tag) { $featrec{olt} = $old_locus_tag; }
      if($proteinid) { $featrec{proteinid} = $proteinid; }
      if($dbxref) { $featrec{dbxref} = $dbxref; }
      $rethash{$id} = {%featrec};
    }
  }
  return(%rethash);
}
# }}}

# {{{ sub gbkfile2seqobj %(file) returns a Bio::Seq object
sub gbkfile2seqobj {
  my $self = shift(@_);
  my %args = @_;
  my $filename = $gbkDir . '/' . $args{file};
  unless(-e $filename) { $filename = $args{file}; }
  my $seqio = Bio::SeqIO->new(-file => $filename);
  my $seqobj = $seqio->next_seq();
  return($seqobj);
}
# }}}

# {{{ sub organism2taxonomy %(organism, partial) returns %(organism, taxonomy)
sub organism2taxonomy {
my $self = shift(@_);
my %args = @_;
my $qstr;
if($args{partial}) {
$qstr = qq/select organism, taxonomy from orgtab where organism like '$args{organism}%' limit 1/;
} else {
$qstr = qq/select organism, taxonomy from orgtab where organism = '$args{organism}' limit 1/;
}
my ($organism, $taxonomy) = $handle->selectrow_array($qstr);
if($handle->state()) { print(STDERR "$qstr\n");
return();
}
else {
return(organism => $organism, taxonomy => $taxonomy);
}
}
# }}}

# {{{ sub gbkfile2molecule %(file) return scalar 'P' or 'G'
sub gbkfile2molecule {
my $self = shift(@_);
my %args = @_;
my $filename = $args{file};
my $qfn = $handle->quote($filename);
my $qstr = qq/select definition from orgtab where filename = $qfn/;
my ($definition) = $handle->selectrow_array($qstr);
my $porg;
if($definition=~m/plasmid/i) {
  $porg = 'P';
}
else {
  $porg = 'G';
}
return($porg);
}
# }}}

# {{{ sub organism2molecules {
sub organism2molecules {
my $self = shift(@_);
my %args = @_;
my $organism = $args{organism};
my $qorg = $handle->quote($organism);
my $qstr = qq/select definition from orgtab where organism = $qorg/;
my $stmt = $handle->prepare($qstr);
$stmt->execute();
my @molecules;
while(my $hr = $stmt->fetchrow_hashref()) {
  my $definition = $hr->{definition};
  if($definition=~m/plasmid/i) {
    push(@molecules, 'P');
  }
  else { push(@molecules, 'G'); }
}
return(@molecules);
}
# }}}

# {{{ sub orgFileMol {
sub orgFileMol {
my $self = shift(@_);
my %args = @_;
my $organism = $args{organism};
my $qorg = $handle->quote($organism);
my $qstr = qq/select definition, filename, nt_len from orgtab where organism = $qorg/;
$qstr .= qq/ and masked != 1 order by nt_len desc/;
my $stmt = $handle->prepare($qstr);
unless($stmt) {
croak("Failed to prepare: $qstr\n");
}
$stmt->execute();
my @molecules;
while(my $hr = $stmt->fetchrow_hashref()) {
  my $definition = $hr->{definition};
  my $file = $hr->{filename};
  my $nt_len = $hr->{nt_len};
  if($definition=~m/plasmid/i) {
    push(@molecules, { file => $file,  molecule => 'P', nt_len => $nt_len });
  }
  else { push(@molecules, { file => $file,  molecule => 'G', nt_len => $nt_len });}
}
$stmt->finish();
return(@molecules);
}
# }}}

# {{{ sub gbkfile2organism %(file) returns a hash;
sub gbkfile2organism {
my $self = shift(@_);
my %args = @_;
my $qstr = qq/select * from orgtab where filename = '$args{file}'/;
my $stmt = $handle->prepare($qstr);
$stmt->execute();
my $hr = $stmt->fetchrow_hashref();
my %rethash = %{$hr};
#my ($organism, $taxonomy) = $handle->selectrow_array($qstr);
#return(organism => $organism, taxonomy => $taxonomy);
return(%rethash);
}
# }}}

# {{{ sub feature2product %(feature) %(locus_tag, product, gene)
sub feature2product {
my $self = shift(@_);
my %args = @_;
my $feat = $args{feature};
my $locus_tag;
if($feat->has_tag('locus_tag')) {
my @lts = $feat->get_tag_values('locus_tag');
$locus_tag = $lts[0];
}
my $product;
if($feat->has_tag('product')) {
my @prods = $feat->get_tag_values('product');
$product = join(" ", @prods);
}
my $gene;
if($feat->has_tag('gene')) {
my @genes = $feat->get_tag_values('gene');
$gene = join(" ", @genes);
}
return(locus_tag => $locus_tag, product => $product, gene => $gene);
}
# }}}

# {{{ sub locusTag2aaObj %(file or seqobj, locus_tag) returns a single Bio::Seq object
sub locusTag2aaObj {
  my $self = shift(@_);
  my %args = @_;
  my $seqobj;
  if($args{file}) {
  my $filename = $args{file};
  unless (-e $filename) {
  $filename = $gbkDir . '/' . $args{file};
  }
  my $seqio = Bio::SeqIO->new(-file => $filename);
  $seqobj = $seqio->next_seq();
  }
  elsif($args{seqobj}) {
  $seqobj = $args{seqobj};
  }
  else {
    return(undef);
  }

  my @features = $seqobj->all_SeqFeatures();
  foreach my $feat (@features) {
    unless ($feat->primary_tag() eq 'CDS') { next ; }
    if($feat->has_tag('locus_tag')) {
      my @lts = $feat->get_tag_values('locus_tag');
      if(grep {$_ eq $args{locus_tag}} @lts) {
        my $aaobj = _feat_translate($feat);
        my $product;
        if($feat->has_tag('product')) {
          $product = join(" ", $feat->get_tag_values('product'));
        }
        #print("=== $product\n");
        $aaobj->display_name($args{locus_tag});
        $aaobj->description($product);
        return($aaobj);
      }
    }
  }
  return(undef);
}
# }}}

# {{{ sub neighbours %(file or seqobj, locus_tag, leftn, rightn) returns @(feature objects)
sub neighbours {
  my $self = shift(@_);
  my %args = @_;
  my $seqobj;
  if($args{file}) {
  my $filename = $gbkDir . '/' . $args{file};
  my $seqio = Bio::SeqIO->new(-file => $filename);
  $seqobj = $seqio->next_seq();
  }
  elsif($args{seqobj}) {
  $seqobj = $args{seqobj};
  }

  my @temp = $seqobj->all_SeqFeatures();
  my @temp1;
  foreach my $feat (@temp) {
    my $pritag = $feat->primary_tag();
    if($pritag eq 'CDS' or $pritag=~m/RNA/) { 
      push(@temp1, $feat);
      }
  }
  my @features = sort _feat_sorter @temp1;
  my $featCnt = 0;
  my $ltSerial;
  foreach my $feat (@features) {
    if($feat->has_tag('locus_tag')) {
    my @lts = $feat->get_tag_values('locus_tag');
    if(grep {$_ eq $args{locus_tag}} @lts) {
      $ltSerial = $featCnt;
    }
    }
    $featCnt += 1;
  }
  my $leftSerial = $ltSerial - $args{leftn};
  my $rightSerial = $ltSerial + $args{rightn};
  if($leftSerial < 1) { $leftSerial = 1; }
  if($rightSerial > $featCnt) { $rightSerial = $featCnt; }
  my @wantedFeatures = @features[$leftSerial..$rightSerial];
  #print(STDERR "--- $leftSerial\t$ltSerial\t$rightSerial\n");
  return(@wantedFeatures);
}
# }}}

# {{{ sub genbank2blastnDB %([files], name, title) returns %(name, files);
sub genbank2blastnDB {
my $self = shift(@_);
my %args = @_;
# print(STDERR "$args{organism}\t$args{name}\n");

my($fh, $fn)=tempfile($template, DIR => $tempdir, SUFFIX => '.fna');
my $seqout = Bio::SeqIO->new(-fh => $fh, -format => 'fasta');

my @gbkNames = @{$args{files}};

foreach my $temp (@gbkNames)  {
my $filename;
if($temp=~m/\//) { $filename = $temp; }
else {
$filename = $gbkDir . '/' . $temp;
}

my $dispname = $temp; $dispname=~s/\.gbk$//;
my $seqio = Bio::SeqIO->new(-file => $filename);
my $seqobj = $seqio->next_seq();
$seqobj->display_name($dispname);
$seqout->write_seq($seqobj);
}
close($fh);
qx(/usr/local/bin/makeblastdb -in $fn -title "$args{title}" -dbtype nucl -out $args{name});
#print(STDERR "$fn\n");
unlink($fn);
return(name => $args{name}, files => [@gbkNames]);
}
# }}}

# {{{ sub filenames2blastnDB %([files], name, title) returns %(name, files);
sub filenames2blastnDB {
my $self = shift(@_);
my %args = @_;
# print(STDERR "$args{organism}\t$args{name}\n");

my($fh, $fn)=tempfile($template, DIR => $tempdir, SUFFIX => '.fna');
my $seqout = Bio::SeqIO->new(-fh => $fh, -format => 'fasta');

my @gbkNames = @{$args{files}};

foreach my $temp (@gbkNames)  {
my $filename;
if($temp=~m/\//) { $filename = $temp; }
else {
$filename = $gbkDir . '/' . $temp;
}

my $dispname = $temp; $dispname=~s/\.gbk$//;
my $seqio = Bio::SeqIO->new(-file => $filename);
my $seqobj = $seqio->next_seq();
$seqobj->display_name($dispname);
$seqout->write_seq($seqobj);
}
close($fh);
qx(/usr/local/bin/makeblastdb -in $fn -title "$args{title}" -dbtype nucl -out $args{name});
#print(STDERR "$fn\n");
unlink($fn);
return(name => $args{name}, files => [@gbkNames]);
}
# }}}

# {{{ sub organism2blastnDB %(organism, name) returns %(name, files);
sub organism2blastnDB {
my $self = shift(@_);
my %args = @_;
# print(STDERR "$args{organism}\t$args{name}\n");
my $qstr = qq/select filename from orgtab where organism = '$args{organism}'/;
my $stmt = $handle->prepare($qstr);
$stmt->execute();
my($fh, $fn)=tempfile($template, DIR => $tempdir, SUFFIX => '.fna');
my $seqout = Bio::SeqIO->new(-fh => $fh, -format => 'fasta');
my @gbkNames;
while(my $hr = $stmt->fetchrow_hashref()) {
my $filename = $gbkDir . '/' . $hr->{filename};
push(@gbkNames, $filename);
}
foreach my $temp (@gbkNames)  {
my $filename;
if($temp=~m/\//) { $filename = $temp; }
else {
$filename = $gbkDir . '/' . $temp;
}
my $dispname = $temp; $dispname=~s/\.gbk$//;
my $seqio = Bio::SeqIO->new(-file => $filename);
my $seqobj = $seqio->next_seq();
$seqobj->display_name($dispname);
$seqout->write_seq($seqobj);
}
close($fh);
qx(/usr/local/bin/makeblastdb -in $fn -title "$args{organism}" -dbtype nucl -out $args{name});
#print(STDERR "$fn\n");
unlink($fn);
return(name => $args{name}, files => [@gbkNames]);
}
# }}}

# {{{ sub genbank2blastpDB %([files], name, title, faafn)
# returns %(name, [files], faafn);
# If you supply a faafn then it is your responsibility to unlink it.
sub genbank2blastpDB {
my $self = shift(@_);
my %args = @_;

my($fh, $fn);
my $keepfaa = 0;
if($args{faafn} =~ m/\w+/) {
$fn = $args{faafn};
open($fh, ">", $fn);
$keepfaa = 1;
}
else {
($fh, $fn)=tempfile($template, DIR => $tempdir, SUFFIX => '.faa');
}

my $seqout = Bio::SeqIO->new(-fh => $fh, -format => 'fasta');

my @gbkNames = @{$args{files}};

my $cdsCnt = 0;
foreach my $temp (@gbkNames)  {
  my $draftName = $draftDir . '/' . $temp;
  my $finName = $gbkDir . '/' . $temp;
  my $filename;
  if(-e $temp) { $filename = $temp; }
  elsif (-e $finName) { $filename = $finName; }
  elsif (-e $draftName) { $filename = $draftName; }
  else { croak("Could not find $temp genbank file\n"); }
my $seqio = Bio::SeqIO->new(-file => $filename);
my $seqobj = $seqio->next_seq();
  foreach my $feature ($seqobj->all_SeqFeatures()) {
    if($feature->primary_tag() eq 'CDS') {
      $cdsCnt += 1;
      my $product;
      my $id;
      my $gene;
      my @tags = $feature->get_all_tags();
      foreach my $tag (@tags) {
        if($tag=~m/product/i) {
          $product=join(" ", $feature->get_tag_values($tag));
        }
        if($tag=~m/^locus_tag/) {
          my $lt=join("|", $feature->get_tag_values($tag));
          $id = $lt;
        }
        if($tag=~m/gene/) {
          my $temp=join("|", $feature->get_tag_values($tag));
          $gene = $temp;
        }
      }
      my $fr = $feature->strand() == 1 ? 'F' : 'R';
      unless($id) { $id = $fr . "_CDS_at_" . $feature->start(); }
      my $aaobj = _feat_translate($feature);
      $aaobj->display_name($id);
      if($gene) { $product .= " gene: $gene"; }
      $aaobj->description($product);
      $seqout->write_seq($aaobj);
    }
  }
#$emblout->write_seq($seqobj);
}
close($fh);
if($cdsCnt) {
  my $runbin = $blastbindir ."/makeblastdb";
qx($runbin -in $fn -title "$args{title}" -dbtype prot -out $args{name});
# print(STDERR qq($runbin -in $fn -title "$args{title}" -dbtype prot -out $args{name}), "\n");
}
#print(STDERR "$fn\n");
unless($keepfaa) { unlink($fn); }
if($cdsCnt) {
  if($keepfaa) {
return(name => $args{name}, files => [@gbkNames], faafn => $fn);
  }
  else {
return(name => $args{name}, files => [@gbkNames]);
  }
}
else {
  return();
}
}
# }}}

# {{{ sub genbank2blDBFTFaa %(file, name, title) returns %(name, faa, feats, seqobj);
sub genbank2blDBFTFaa {
my $self = shift(@_);
my %args = @_;
my($fh, $fn)=tempfile($template, DIR => $tempdir, SUFFIX => '.faa');
my $seqout = Bio::SeqIO->new(-fh => $fh, -format => 'fasta');

my $cdsCnt = 0;
my $temp = $args{file};
  my $draftName = $draftDir . '/' . $temp;
  my $finName = $gbkDir . '/' . $temp;
  my $filename;
  if(-e $temp) { $filename = $temp; }
  elsif (-e $finName) { $filename = $finName; }
  elsif (-e $draftName) { $filename = $draftName; }
  else { croak("Could not find $temp genbank file\n"); }
my $seqio = Bio::SeqIO->new(-file => $filename);
my $seqobj = $seqio->next_seq();

# {{{ FT
my %feathash; # for the FT block
FT: {
  my @temp = $seqobj->all_SeqFeatures();
  my @features = sort _feat_sorter(@temp);
  my $featCnt = 0;
  foreach my $feat (@features) {
    my $pritag = $feat->primary_tag();
    my $start_pos = $feat->start();
    my $end_pos = $feat->end();
    my $strand = $feat->strand();
    
    if($pritag eq 'CDS') {
      $featCnt += 1;
      my $locus_tag;
      if($feat->has_tag('locus_tag')) {
        my @lts = $feat->get_tag_values('locus_tag');
        $locus_tag = $lts[0];
      }
      my $id = defined($locus_tag) ? $locus_tag : $pritag . "_" . $featCnt;
      $feathash{$id} = $feat;
    }
  }
}
# }}}


  foreach my $feature ($seqobj->all_SeqFeatures()) {
    if($feature->primary_tag() eq 'CDS') {
      $cdsCnt += 1;
      my $product;
      my $id;
      my $gene;
      my @tags = $feature->get_all_tags();
      foreach my $tag (@tags) {
        if($tag=~m/product/i) {
          $product=join(" ", $feature->get_tag_values($tag));
        }
        if($tag=~m/^locus_tag/) {
          my $lt=join("|", $feature->get_tag_values($tag));
          $id = $lt;
        }
        if($tag=~m/gene/) {
          my $temp=join("|", $feature->get_tag_values($tag));
          $gene = $temp;
        }
      }
      my $fr = $feature->strand() == 1 ? 'F' : 'R';
      unless($id) { $id = $fr . "_CDS_at_" . $feature->start(); }
      my $aaobj = _feat_translate($feature);
      $aaobj->display_name($id);
      if($gene) { $product .= " gene: $gene"; }
      $aaobj->description($product);
      $seqout->write_seq($aaobj);
    }
  }
#$emblout->write_seq($seqobj);


close($fh);
if($cdsCnt) {
qx(/usr/local/bin/makeblastdb -in $fn -title "$args{title}" -dbtype prot -out $args{name});
print(STDERR qq(/usr/local/bin/makeblastdb -in $fn -title "$args{title}" -dbtype prot -out $args{name}), "\n");
}
#print(STDERR "$fn\n");
#unlink($fn);
if($cdsCnt) {
return(name => $args{name}, faa => $fn, feats => \%feathash, seqobj => $seqobj);
}
else {
  return();
}
}
# }}}

# {{{ sub genbank2faa %([files], old_locus_tag, ofh, tfh) returns %(name, [files]);
# tfh is optional. If a filehandle is given then a table of CDS is written
# to that filehandle and the filehandle closed.
sub genbank2faa {
my $self = shift(@_);
my %args = @_;

my $ofh = $args{ofh};
my $seqout = Bio::SeqIO->new(-fh => $ofh, -format => 'fasta');

my $tfh;
my $tableWanted = 0;
if(exists($args{tfh})) {
$tableWanted = 1;
$tfh = $args{tfh};
}

my @gbkNames = @{$args{files}};

my $cdsCnt = 0;
foreach my $temp (@gbkNames)  {
  my $filename;
  if(-e $temp or -l $temp) { $filename = $temp; }
  else { $filename = $gbkDir . '/' . $temp; }
  unless(-r $filename) {
    carp("Could not read $filename\n");
  }
  if(-z $filename) {
    carp("Zero size of $filename\n");
  }
my $seqio = Bio::SeqIO->new(-file => $filename);
my ($gbkBn, $directory, $ext) = fileparse($filename, qr/\.[^.]*/);
while(my $seqobj = $seqio->next_seq()) {
  my $seqid = $seqobj->display_id();
my $species = $seqobj->species();
my $binomial;
if($species) {
$binomial = $species->binomial('FULL');
}
  my @temp = $seqobj->all_SeqFeatures();
  foreach my $feature (sort _feat_sorter @temp) {
    if($feature->primary_tag() eq 'CDS') {
      $cdsCnt += 1;
      my $product;
      my $id;
      my $gene;

# If old_locus_tag is specified then attempt to get the old_locus_tag
# failing which, get the locus_tag.
          if($args{old_locus_tag}) {
            if($feature->has_tag("old_locus_tag")) {
              my $lt=join("|", $feature->get_tag_values("old_locus_tag"));
              $id = $lt;
            }
            elsif($feature->has_tag("locus_tag")) {
              my $lt=join("|", $feature->get_tag_values("locus_tag"));
              $id = $lt;
            }
          }
          elsif($feature->has_tag("locus_tag")) {
            my $lt=join("|", $feature->get_tag_values("locus_tag"));
            $id = $lt;
          }

# Get product and gene.
          if($feature->has_tag("product")) {
            $product=join(" ", $feature->get_tag_values("product"));
          }
          if($feature->has_tag("gene")) {
            $gene=join(" ", $feature->get_tag_values("gene"));
          }

      my $lig = $feature->start() . ":" . $feature->end(); # location in genbank
      $lig .= ":" . $feature->strand();
      my $fr = $feature->strand() == 1 ? 'F' : 'R';
      unless($id) { $id = $gbkBn . "_" . sprintf("%05d", $cdsCnt); }
      my $aaobj = _feat_translate($feature);
      if($aaobj) {
      $aaobj->display_name($id);
      if($gene) { $product .= " gene: $gene"; }
      if($binomial) {
      my $desc = $seqid . " " . $product . " : " . $binomial;
      $aaobj->description($desc);
      }
      else {
      my $desc = $seqid . " " . $product;
      $aaobj->description($desc);
      }
      $seqout->write_seq($aaobj);
      if($tableWanted) {
        print($tfh join("\t", $id, $aaobj->length(), $feature->start(),
        $feature->end(), $feature->strand(), $product), "\n");
      }
      }
      else { carp("Failed to translate $id from $filename\n"); }
    }
  }
}
}
close($ofh);
if($tableWanted) {
close($tfh);
}
}
# }}}

# {{{ sub genbank16sfna %([files], ofh, [ids] ) returns %(name, [files]);
# ids, if given should have to same number of elements as the number of files.
# to that filehandle and the filehandle closed.
sub genbank16sfna {
  my $self = shift(@_);
  my $errh = $self->{errh};
  my %args = @_;
  my @ids = @{$args{ids}};

  my $ofh = $args{ofh};
  my $seqout = Bio::SeqIO->new(-fh => $ofh, -format => 'fasta');

  my @gbkNames = @{$args{files}};

my $ndx = 0;

  foreach my $temp (@gbkNames)  {
    my $filename;
    if(-e $temp or -l $temp) { $filename = $temp; }
    else { $filename = $gbkDir . '/' . $temp; }
    unless(-r $filename) {
      carp("Could not read $filename\n");
    }
    if(-z $filename) {
      carp("Zero size of $filename\n");
    }
    my $seqio = Bio::SeqIO->new(-file => $filename);
    my ($gbkBn, $directory, $ext) = fileparse($filename, qr/\.[^.]*/);
    my $seqobj = $seqio->next_seq();
    my @temp = $seqobj->all_SeqFeatures();
    my $rrnaCnt = 0;
    foreach my $feature (sort _feat_sorter @temp) {
      if($feature->primary_tag() eq 'rRNA'
        and $feature->length() > 1200
        and $feature->length() < 1800
        ) {
        $rrnaCnt += 1;
        my $product;
        if($feature->has_tag("product")) {
          my @temp = $feature->get_tag_values("product");
          $product = join(", ", @temp);
        }
        my $id;
        if($ids[$ndx]) {
        $id = $ids[$ndx] . "_" . $rrnaCnt;
        }
        else {
        $id = $gbkBn . "_" . $rrnaCnt;
        }
        my $featobj=$feature->spliced_seq(-nosort => '1');
        if($featobj) {
          $featobj->display_name($id);
          if($product) {
            $featobj->description($product);
          }
          $seqout->write_seq($featobj);
        }
        else {
          print($errh "Failed to get a feature object for $id\n");
        }
      }
    }
#$emblout->write_seq($seqobj);
    $ndx += 1;
  }
}
# }}}

# {{{ sub genbank2protfna %([files], old_locus_tag, ofh) returns %(name, [files]);
# tfh is optional. If a filehandle is given then a table of CDS is written
# to that filehandle and the filehandle closed.
sub genbank2protfna {
  my $self = shift(@_);
  my $errh = $self->{errh};
  my %args = @_;

  my $ofh = $args{ofh};
  my $seqout = Bio::SeqIO->new(-fh => $ofh, -format => 'fasta');

  my $tfh;
  my $tableWanted = 0;
  if(exists($args{tfh})) {
    $tableWanted = 1;
    $tfh = $args{tfh};
  }

  my @gbkNames = @{$args{files}};

  my $cdsCnt = 0;
  foreach my $temp (@gbkNames)  {
    my $filename;
    if(-e $temp or -l $temp) { $filename = $temp; }
    else { $filename = $gbkDir . '/' . $temp; }
    unless(-r $filename) {
      carp("Could not read $filename\n");
    }
    if(-z $filename) {
      carp("Zero size of $filename\n");
    }
    my $seqio = Bio::SeqIO->new(-file => $filename);
    my ($gbkBn, $directory, $ext) = fileparse($filename, qr/\.[^.]*/);
    while(my $seqobj = $seqio->next_seq()) {
      my $seqid = $seqobj->display_id();
      my $species = $seqobj->species();
      my $binomial;
      if($species) {
        $binomial = $species->binomial('FULL');
      }
      my @temp = $seqobj->all_SeqFeatures();
      foreach my $feature (sort _feat_sorter @temp) {
        if($feature->primary_tag() eq 'CDS') {
          $cdsCnt += 1;
          my $product;
          my $id;
          my $gene;

# If old_locus_tag is specified then attempt to get the old_locus_tag
# failing which, get the locus_tag.
          if($args{old_locus_tag}) {
            if($feature->has_tag("old_locus_tag")) {
              my $lt=join("|", $feature->get_tag_values("old_locus_tag"));
              $id = $lt;
            }
            elsif($feature->has_tag("locus_tag")) {
              my $lt=join("|", $feature->get_tag_values("locus_tag"));
              $id = $lt;
            }
          }
          elsif($feature->has_tag("locus_tag")) {
            my $lt=join("|", $feature->get_tag_values("locus_tag"));
            $id = $lt;
          }

# Get product and gene.
          if($feature->has_tag("product")) {
            $product=join(" ", $feature->get_tag_values("product"));
          }
          if($feature->has_tag("gene")) {
            $gene=join(" ", $feature->get_tag_values("gene"));
          }


          my $fr = $feature->strand() == 1 ? 'F' : 'R';
          unless($id) { $id = $gbkBn . "_" . $cdsCnt; }
          my $featobj=$feature->spliced_seq(-nosort => '1');
          if($featobj) {
            $featobj->display_name($id);
            if($gene) { $product .= " gene: $gene"; }
            if($binomial) {
              $featobj->description($seqid . " " . $product . " : " . $binomial);
            }
            else {
              $featobj->description($seqid . " " . $product);
            }
            $seqout->write_seq($featobj);
          }
          else {
            print($errh "Failed to get a feature object for $id\n");
          }
        }
      }
    }
  }
  close($ofh);
}
# }}}

# {{{ sub organism2blastpDB %(organism, name, keepfasta) returns %(name, files, fasta);
sub organism2blastpDB {
my $self = shift(@_);
my %args = @_;
# print(STDERR "$args{organism}\t$args{name}\n");
my $qstr = qq/select filename from orgtab where organism = '$args{organism}' and masked = 0/;
my $stmt = $handle->prepare($qstr);
$stmt->execute();
my($fh, $fn)=tempfile($template, DIR => $tempdir, SUFFIX => '.faa');
my $seqout = Bio::SeqIO->new(-fh => $fh, -format => 'fasta');

my @gbkNames;

while(my $hr = $stmt->fetchrow_hashref()) {
my $filename = $gbkDir . '/' . $hr->{filename};
push(@gbkNames, $filename);
my $seqio = Bio::SeqIO->new(-file => $filename);
my $seqobj = $seqio->next_seq();
my $cdsCnt = 0;
  foreach my $feature ($seqobj->all_SeqFeatures()) {
    if($feature->primary_tag() eq 'CDS') {
      $cdsCnt += 1;
      my $product;
      my $id;
      my $gene;
      my @tags = $feature->get_all_tags();
      foreach my $tag (@tags) {
        if($tag=~m/product/i) {
          $product=join(" ", $feature->get_tag_values($tag));
        }
        if(lc($tag) eq "locus_tag" or lc($tag) eq "systematic_id") {
        # if($tag=~m/locus_tag|systematic_id/) {
          my $lt=join("|", $feature->get_tag_values($tag));
          $id = $lt;
        }
        if($tag=~m/gene/) {
          my $temp=join("|", $feature->get_tag_values($tag));
          $gene = $temp;
        }
      }
      my $aaobj = _feat_translate($feature);
      $aaobj->display_name($id);
      if($gene) { $product .= " gene: $gene"; }
      $aaobj->description($product);
      $seqout->write_seq($aaobj);
    }
  }
#$emblout->write_seq($seqobj);
}
close($fh);
qx(/usr/local/bin/makeblastdb -in $fn -title "$args{organism}" -dbtype prot -out $args{name});
#print(STDERR "$fn\n");
if($args{keepfasta}) {
return(name => $args{name}, files => [@gbkNames], fasta => $fn);
}
else {
unlink($fn);
return(name => $args{name}, files => [@gbkNames]);
}

}
# }}}

# {{{ sub getHandle
sub getHandle {
my $self = shift(@_);
return($handle);
}
# }}}

# {{{ sub getDraftsHandle
sub getDraftsHandle {
my $self = shift(@_);
return($draftsHandle);
}
# }}}

# {{{ seq2obj (hash(seq, id, description)) returns(Bio::Seq); 
sub seq2obj {
my $self = shift(@_);
my %args = @_;
my $retobj = Bio::Seq->new(-seq => $args{seq});
$retobj->display_name($args{id});
if($args{description}) {
$retobj->description($args{description});
}
return($retobj);
}
# }}}

# {{{ sub upstreamOfLocusTag (hash(file|seqobj, locus_tag, length)) returns(Bio::Seq);
# A fixed length of sequence upstream of a given feature.
sub upstreamOfLocusTag {
  my $self = shift(@_);
  my %args = @_;
  my $locusTag = $args{locus_tag};
  my $wantLen = $args{length};
  my $seqobj;
  if($args{seqobj}) {
    $seqobj = $args{seqobj};
  }
  elsif($args{file}) {
    my $gbkfile = $args{file};
    my $seqio = Bio::SeqIO->new(-file => $gbkfile);
    $seqobj=$seqio->next_seq();
  }

  foreach my $feature ($seqobj->all_SeqFeatures()) {
    if($feature->primary_tag() eq 'CDS') {
      if($feature->has_tag('locus_tag')) {
        my @lts = $feature->get_tag_values('locus_tag');
        if(grep {$_ eq $locusTag} @lts) {
          my $start = $feature->start();
          my $end = $feature->end();
          my $strand = $feature->strand();
          my $subseq;
          if($strand == 1 ) {
            $subseq = $seqobj->trunc($start-$wantLen, $start-1);
          }
          elsif($strand == -1 ) {
            $subseq = $seqobj->trunc($end+1, $end+$wantLen)->revcom();
          }
          $subseq->display_name("$locusTag");
          $subseq->description("$wantLen upstream nucleotides");
          return($subseq);
        }
      }
    }
  }
return();
}
# }}}

# {{{ sub freeUpstreamOfLocusTag (hash(file, locus_tag, minlen)) returns(Bio::Seq);
# Region upstream of feature which is devoid if any CDS or tRNA or rRNA features.
sub freeUpstreamOfLocusTag {
  my $self = shift(@_);
  my %args = @_;
  my $gbkfile = $args{file};
  my $locusTag = $args{locus_tag};
  my $minLen = 20;
  if(exists($args{minlen})) {$minLen = $args{minlen};}
  my $seqio = Bio::SeqIO->new(-file => $gbkfile);
  my $seqobj=$seqio->next_seq();
  my @two = (undef, undef);
  foreach my $feature ($seqobj->all_SeqFeatures()) {
    unless($feature->primary_tag()=~m/CDS|[rt]RNA/) { next; }
    shift(@two); push(@two, $feature);
    unless(defined($two[0]) and defined($two[1])) { next; }
    my $subseq;
    if($two[1]->has_tag('locus_tag') and
      grep {$_ eq $locusTag} $two[1]->get_tag_values('locus_tag') and
      $two[1]->strand() == 1
      ) {
      if(($two[1]->start() - 1) - ($two[0]->end() + 1) >= $minLen ) {
      $subseq = $seqobj->trunc($two[0]->end() + 1, $two[1]->start() - 1);
      }
    }
    elsif($two[0]->has_tag('locus_tag') and
      grep {$_ eq $locusTag} $two[0]->get_tag_values('locus_tag') and
      $two[0]->strand() == -1 ) {
      if(($two[1]->start() - 1) - ($two[0]->end() + 1) >= $minLen ) {
      $subseq = $seqobj->trunc($two[0]->end() + 1, $two[1]->start() - 1)->revcom();
      }
    }
    if($subseq) {
    if($subseq->length() >= $minLen) {
          $subseq->display_name("$locusTag");
          $subseq->description("Upstream nucleotides");
          return($subseq);
    }
    else { return(); }
    }
    else { next; }
  }
return();
}
# }}}

# {{{ sub freeUpstreamOfLocusTagMulti (hash(file, [locus_tags], minlen)) returns(Bio::Seq);
# Region upstream of feature which is devoid if any CDS or tRNA or rRNA features.
sub freeUpstreamOfLocusTagMulti {
  my $self = shift(@_);
  my %args = @_;
  my $gbkfile = $args{file};
  my $ltref = $args{locus_tags};
  my @lts = @{$ltref};
  my $minLen = 20;
  if(exists($args{minlen})) {$minLen = $args{minlen};}
  my $seqio = Bio::SeqIO->new(-file => $gbkfile);
  my $seqobj=$seqio->next_seq();
  my @features = $seqobj->all_SeqFeatures();
  my @retlist;
  LT: for my $locusTag (@lts) {
#    print(STDERR $locusTag, "         \r");
  my @two = (undef, undef);
  my $ltenc = 0;
  foreach my $feature (@features) {
    my $pritag = $feature->primary_tag();
    unless($pritag=~m/CDS|[rt]RNA/) { next; }
    shift(@two); push(@two, $feature);
    unless(defined($two[0]) and defined($two[1])) { next; }
    my $subseq;
    if($two[1]->has_tag('locus_tag') and
      grep {$_ eq $locusTag} $two[1]->get_tag_values('locus_tag') and
      $two[1]->strand() == 1
      ) {
      $ltenc = 1;
      if(($two[1]->start() - 1) - ($two[0]->end() + 1) >= $minLen ) {
      $subseq = $seqobj->trunc($two[0]->end() + 1, $two[1]->start() - 1);
      }
    }
    elsif($two[0]->has_tag('locus_tag') and
      grep {$_ eq $locusTag} $two[0]->get_tag_values('locus_tag') and
      $two[0]->strand() == -1 ) {
      $ltenc = 1;
      if(($two[1]->start() - 1) - ($two[0]->end() + 1) >= $minLen ) {
      $subseq = $seqobj->trunc($two[0]->end() + 1, $two[1]->start() - 1)->revcom();
      }
    }
    if($subseq) {
    if($subseq->length() >= $minLen) {
      my $seqstr = $subseq->seq();
      my $pushobj = Bio::Seq->new(-seq => $seqstr);
      $pushobj->display_name($locusTag);
      $pushobj->description("Upstream nucleotides");
      push(@retlist, $pushobj);
    }
    }
    if($ltenc) {
      next LT;
    }
    else { next; }
  }
  }
return(@retlist);
}
# }}}

# {{{ sub locusTagCodonPos %(file or seqobj, locus_tag or \@locus_tag, codon) returns a list of integers
sub locusTagCodonPos {
  my $self = shift(@_);
  my %args = @_;
  my $seqobj;
  if($args{file}) {
  my $temp = $gbkDir . '/' . $args{file};
  my $filename;
  if(-e $temp) { $filename = $temp; } else { $filename = $args{file}; }
  my $seqio = Bio::SeqIO->new(-file => $filename);
  $seqobj = $seqio->next_seq();
  }
  elsif($args{seqobj}) {
  $seqobj = $args{seqobj};
  }
  my $codon = 'tta';
  if($args{codon}) { $codon = $args{codon}; }

  my @features = $seqobj->all_SeqFeatures();
  foreach my $feat (@features) {
    unless ($feat->primary_tag() eq 'CDS') { next ; }
    if($feat->has_tag('locus_tag')) {
      my @lts = $feat->get_tag_values('locus_tag');
      if(grep {$_ eq $args{locus_tag}} @lts) {
        my $featObj = $feat->spliced_seq(-nosort => 1);
        my $featSeq = $featObj->seq();
        my $pos = 0;
        my @retlist;
        while(my $triplet = lc(substr($featSeq, $pos, 3))) {
          if($triplet eq lc($codon)) {
            push(@retlist, $pos + 1);
          }
          $pos += 3;
        }
        return(@retlist);
      }
    }
  }
  return(undef);
}
# }}}

# {{{ sub CDSstartCodonPos %(file or seqobj, start, codon) returns a list of integers
# Often CDSs do not have locus_tags. Here we use the start position of the
# CDS feature to identify it. Reports all the positions at which a given
# codon exists in the CDS starting at the given start position.
sub CDSstartCodonPos {
  my $self = shift(@_);
  my %args = @_;
  my $seqobj;
  if($args{file}) {
  my $temp = $gbkDir . '/' . $args{file};
  my $filename;
  if(-e $temp) { $filename = $temp; } else { $filename = $args{file}; }
  my $seqio = Bio::SeqIO->new(-file => $filename);
  $seqobj = $seqio->next_seq();
  }
  elsif($args{seqobj}) {
  $seqobj = $args{seqobj};
  }
  my $codon = 'tta';
  if($args{codon}) { $codon = $args{codon}; }

  my @features = $seqobj->all_SeqFeatures();
  foreach my $feat (@features) {
    unless ($feat->primary_tag() eq 'CDS') { next ; }
    if($feat->start() == $args{start}) {
        my $featObj = $feat->spliced_seq(-nosort => 1);
        my $featSeq = $featObj->seq();
        my $pos = 0;
        my @retlist;
        while(my $triplet = lc(substr($featSeq, $pos, 3))) {
          if($triplet eq lc($codon)) {
            push(@retlist, $pos + 1);
          }
          $pos += 3;
        }
        return(@retlist);
      }
  }
  return(undef);
}
# }}}

# {{{ sub CDSTags %(file or seqobj, locus_tag or start) returns a hash;
sub CDSTags {
  my $self = shift(@_);
  my %args = @_;
  my $seqobj;
  if($args{file}) {
    my $temp = $gbkDir . '/' . $args{file};
    my $filename;
    if(-e $temp) { $filename = $temp; } else { $filename = $args{file}; }
    my $seqio = Bio::SeqIO->new(-file => $filename);
    $seqobj = $seqio->next_seq();
  }
  elsif($args{seqobj}) {
    $seqobj = $args{seqobj};
  }
  if(exists $args{start}) {
    my @features = $seqobj->all_SeqFeatures();
    foreach my $feat (@features) {
      unless ($feat->primary_tag() eq 'CDS') { next ; }
      if($feat->start() == $args{start}) {
        my %rethash;
        my $featObj = $feat->spliced_seq(-nosort => 1);
        my $featSeq = $featObj->seq();
        $rethash{seq} = $featSeq;
        my @tags = $feat->get_all_tags();
        foreach my $tag (@tags) {
          my $valstr = join("; ", $feat->get_tag_values($tag));
          $rethash{$tag} = $valstr;
        }
        return(%rethash);
      }
    }
  }

  if(exists $args{locus_tag}) {
    my @features = $seqobj->all_SeqFeatures();
    foreach my $feat (@features) {
      unless ($feat->primary_tag() eq 'CDS' and $feat->has_tag('locus_tag')) { next ; }
      my @lts = $feat->get_tag_values('locus_tag');
      if(grep {$_ eq $args{locus_tag} } @ lts) {
        my %rethash;
        my $featObj = $feat->spliced_seq(-nosort => 1);
        my $featSeq = $featObj->seq();
        $rethash{seq} = $featSeq;
        my @tags = $feat->get_all_tags();
        foreach my $tag (@tags) {
          my $valstr = join("; ", $feat->get_tag_values($tag));
          $rethash{$tag} = $valstr;
        }
        return(%rethash);
      }
    }
  }

  return(undef);
}
# }}}

# {{{ sub _feat_translate {
sub _feat_translate {
  my $feature=shift(@_);
  my $codon_start=1;
  if($feature->has_tag('codon_start')) {
      ($codon_start) = $feature->get_tag_values('codon_start');
      }
      my $aaobj;
      eval {
      my $offset=1;
      if($codon_start > 1) { $offset = $codon_start;}
      my $featobj=$feature->spliced_seq(-nosort => '1');
      $aaobj = $featobj->translate(-offset => $offset, -complete => 1);
      };
  return($aaobj);
}
# }}}

# {{{ sub _feat_translate_in {
sub _feat_translate_in {
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
# }}}

# {{{ sub _feat_sorter # internal sub
sub _feat_sorter {
return($a->start() <=> $b->start());
}
# }}}

# {{{ sub startPos2aaObj %(file or seqobj, start) returns a single Bio::Seq object
sub startPos2aaObj {
  my $self = shift(@_);
  my %args = @_;
  my $seqobj;
  if($args{file}) {
  my $filename;
  if(-e $args{file}) { $filename = $args{file}; }
  else {
  $filename = $gbkDir . '/' . $args{file};
  }
  my $seqio = Bio::SeqIO->new(-file => $filename);
  $seqobj = $seqio->next_seq();
  }
  elsif($args{seqobj}) {
  $seqobj = $args{seqobj};
  }

  my @features = $seqobj->all_SeqFeatures();
  foreach my $feat (@features) {
    unless ($feat->primary_tag() eq 'CDS' and $feat->start() == $args{start})
    { next; }
    my $alphaStrand = $feat->strand == 1 ? 'F' : 'R';
    my $dn = join("_", $alphaStrand, "CDS", "at", $args{start});
        my $aaobj = _feat_translate($feat);
        my $product;
        if($feat->has_tag('product')) {
          $product = join(" ", $feat->get_tag_values('product'));
        }
        $aaobj->display_name($dn);
        $aaobj->description($product);
        return($aaobj);
  }
  return(undef);
}
# }}}

# {{{ sub stringSearch hash(filename or seqobj, regex). Returns a list of integers.
# This sub does not search both strands. This is by design.
sub stringSearch {
my $self = shift(@_);
my %args = @_;
my $filename;
my $seqObj;
if($args{file}) {
  if(-e $args{file}) { $filename = $args{file}; }
  else { $filename = $gbkDir . '/' . $args{file}; }
  $seqObj = $self->gbkfile2seqobj(file => $args{file});
}
elsif($args{seqobj}) {
$seqObj = $args{seqobj};
}
my $seq = $seqObj->seq();
my $regex = $args{regex};
my @foundAt;
while($seq=~m/$regex/gi) {
push(@foundAt, pos($seq));
}
return(@foundAt);
}
# }}}

# {{{ sub lt2faa %(file, locus_tag, ofh) returns 1
sub lt2faa {
  my $self = shift(@_);
  my %args = @_;
  my ($seqio, $seqobj);
  if($args{file}) {
  my $temp = $args{file};
  my $draftName = $draftDir . '/' . $temp;
  my $finName = $gbkDir . '/' . $temp;
  my $filename;
  if(-e $temp) { $filename = $temp; }
  elsif (-e $finName) { $filename = $finName; }
  elsif (-e $draftName) { $filename = $draftName; }
  else { croak("Could not find $temp genbank file\n"); }
  $seqio = Bio::SeqIO->new(-file => $filename);
  $seqobj = $seqio->next_seq();
  }
  elsif($args{seqobj}) {
  $seqobj = $args{seqobj};
  }

  my $seqout = Bio::SeqIO->new(-fh => $args{ofh}, -format => 'fasta');
  my @features = $seqobj->all_SeqFeatures();
  foreach my $feat (@features) {
    unless ($feat->primary_tag() eq 'CDS') { next ; }
    if($feat->has_tag('locus_tag')) {
      my @lts = $feat->get_tag_values('locus_tag');
      if(grep {$_ eq $args{locus_tag}} @lts) {
        my $aaobj = _feat_translate($feat);
        my $product;
        if($feat->has_tag('product')) {
          $product = join(" ", $feat->get_tag_values('product'));
        }
        #print("=== $product\n");
        $aaobj->display_name($args{locus_tag});
        $aaobj->description($product);
        $seqout->write_seq($aaobj);
      }
    }
  }
  close($args{ofh});
  return(1);
}
# }}}

# {{{ sub lt2fna %(file, locus_tag, ofh) returns 1
sub lt2fna {
  my $self = shift(@_);
  my %args = @_;
  my ($seqio, $seqobj);
  if($args{file}) {
  my $temp = $args{file};
  my $draftName = $draftDir . '/' . $temp;
  my $finName = $gbkDir . '/' . $temp;
  my $filename;
  if(-e $temp) { $filename = $temp; }
  elsif (-e $finName) { $filename = $finName; }
  elsif (-e $draftName) { $filename = $draftName; }
  else { croak("Could not find $temp genbank file\n"); }
  $seqio = Bio::SeqIO->new(-file => $filename);
  $seqobj = $seqio->next_seq();
  }
  elsif($args{seqobj}) {
  $seqobj = $args{seqobj};
  }

  my $seqout = Bio::SeqIO->new(-fh => $args{ofh}, -format => 'fasta');
  my @features = $seqobj->all_SeqFeatures();
  foreach my $feat (@features) {
    unless ($feat->primary_tag() eq 'CDS') { next ; }
    if($feat->has_tag('locus_tag')) {
      my @lts = $feat->get_tag_values('locus_tag');
      if(grep {$_ eq $args{locus_tag}} @lts) {
        my $ntobj = _feat_ntseq($feat);
        my $product;
        if($feat->has_tag('product')) {
          $product = join(" ", $feat->get_tag_values('product'));
        }
        #print("=== $product\n");
        $ntobj->display_name($args{locus_tag});
        $ntobj->description($product);
        $seqout->write_seq($ntobj);
      }
    }
  }
  close($args{ofh});
  return(1);
}
# }}}

# {{{ sub subseq %(file (or seqobj), start, end, revcom) returns Bio::Seq.
sub subseq {
  my $self = shift(@_);
  my %args = @_;
  my ($seqio, $seqobj);
  if($args{file}) {
  my $temp = $args{file};
  my $draftName = $draftDir . '/' . $temp;
  my $finName = $gbkDir . '/' . $temp;
  my $filename;
  if(-e $temp) { $filename = $temp; }
  elsif (-e $finName) { $filename = $finName; }
  elsif (-e $draftName) { $filename = $draftName; }
  else { croak("Could not find $temp genbank file\n"); }
  $seqio = Bio::SeqIO->new(-file => $filename);
  $seqobj = $seqio->next_seq();
  }
  elsif($args{seqobj}) {
  $seqobj = $args{seqobj};
  }

  my $seqstr;
  $seqstr = $seqobj->subseq($args{start}, $args{end});
  my $tempobj = Bio::Seq->new(-seq => $seqstr);
  my $retobj;
  if($args{revcom}) {
  $retobj = $tempobj->revcom();
  }
  else { $retobj = $tempobj; }
  $retobj->display_id("subsequence");
  return($retobj);
}
# }}}

# {{{ sub _feat_ntseq {
sub _feat_ntseq {
  my $feature=shift(@_);
  my $codon_start=1;
  if($feature->has_tag('codon_start')) {
      ($codon_start) = $feature->get_tag_values('codon_start');
      }
      my $offset=1;
      if($codon_start > 1) { $offset = $codon_start;}
      my $featobj=$feature->spliced_seq(-nosort => '1');
  return($featobj);
}
# }}}

# {{{ sub feat_translate {
sub feat_translate {
  my $self = shift(@_);
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
# }}}

# {{{ subroutines tablist, linelist, tabhash meant for internal use.

sub tablist {
my @in = @_;
print(join("\t", @in), "\n");
}
sub tablistE {
my @in = @_;
print(STDERR join("\t", @in), "\n");
}

sub linelist {
my @in = @_;
print(join("\n", @in), "\n");
}
sub linelistE {
my @in = @_;
print(STDERR join("\n", @in), "\n");
}

sub tabhash {
my %in = @_;
for my $key (sort keys(%in)) {
print(join("\t", $key, $in{$key}), "\n");
}
}
sub tabhashE {
my %in = @_;
for my $key (sort keys(%in)) {
print(STDERR join("\t", $key, $in{$key}), "\n");
}
}

# }}}

# {{{ sub genbank2faa2 %([files], ofh, tfh) returns %(name, [files]);
# tfh is optional. If a filehandle is given then a table of CDS is written
# to that filehandle and the filehandle closed.
sub genbank2faa2 {
my $self = shift(@_);
my $errh = $self->{errh};
my %args = @_;

my $ofh = $args{ofh};
my $seqout = Bio::SeqIO->new(-fh => $ofh, -format => 'fasta');

my $tfh;
my $tableWanted = 0;
if(exists($args{tfh})) {
$tableWanted = 1;
$tfh = $args{tfh};
}

my @gbkNames = @{$args{files}};

my $cdsCnt = 0;
foreach my $temp (@gbkNames)  {
  my $filename;
  if(-e $temp or -l $temp) { $filename = $temp; }
  else { $filename = $gbkDir . '/' . $temp; }
  unless(-r $filename) {
    carp("Could not read $filename\n");
  }
  if(-z $filename) {
    carp("Zero size of $filename\n");
  }
my $seqio = Bio::SeqIO->new(-file => $filename);
my ($gbkBn, $directory, $ext) = fileparse($filename, qr/\.[^.]*/);
while(my $seqobj = $seqio->next_seq()) {
  my $seqid = $seqobj->display_id();
my $species = $seqobj->species();
my $binomial;
if($species) {
$binomial = $species->binomial('FULL');
}
  my @temp = $seqobj->all_SeqFeatures();
  foreach my $feature (sort _feat_sorter @temp) {
    if($feature->primary_tag() eq 'CDS') {
      $cdsCnt += 1;
      my $product;
      my $id;
      my $gene;
      my @tags = $feature->get_all_tags();
      foreach my $tag (@tags) {
        if($tag=~m/product/i) {
          $product=join(" ", $feature->get_tag_values($tag));
        }
        if($tag=~m/^locus_tag/) {
          my $lt=join("|", $feature->get_tag_values($tag));
          $id = $lt;
        }
        if($tag=~m/gene/) {
          my $temp=join("|", $feature->get_tag_values($tag));
          $gene = $temp;
        }
      }
      my $lig = $feature->start() . ":" . $feature->end(); # location in genbank
      $lig .= ":" . $feature->strand();
      my $fr = $feature->strand() == 1 ? 'F' : 'R';
      # unless($id) { $id = $gbkBn . "_" . sprintf("%05d", $cdsCnt); }
      $id = $seqid . "_" . sprintf("%05d", $cdsCnt);
      my $aaobj = _feat_translate($feature);
      if($aaobj) {
      $aaobj->display_name($id);
      if($gene) { $product .= " gene: $gene"; }
      if($binomial) {
      my $desc = $lig . " " . $seqid . " " . $product . " : " . $binomial;
      $aaobj->description($desc);
      }
      else {
      my $desc = $lig . " " . $seqid . " " . $product;
      $aaobj->description($desc);
      }
      $seqout->write_seq($aaobj);
      if($tableWanted) {
        print($tfh join("\t", $id, $aaobj->length(), $feature->start(),
        $feature->end(), $feature->strand(), $product), "\n");
      }
      }
      else { print($errh "Failed to translate $id from $filename\n");
      }
    }
  }
}
}
close($ofh);
if($tableWanted) {
close($tfh);
}
}
# }}}

return(1);


__END__
