package Sco::Blast;
use 5.16.0;
use Carp;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::Seq;
use File::Copy;
use File::Temp qw(tempfile tempdir);
my $tempdir = qw(/home/sco/volatile);
my $template="BlastXXXXX";
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;
# $XML::SAX::ParserPackage = 'XML::SAX::PurePerl';

our($AUTOLOAD);
my $blastbindir = qq(/usr/local/bin);

sub new {
        my($class, $self);
        $class=shift(@_);
        $self={};
        bless($self, $class);
        return($self);
}

### more subs go below ###

# {{{ hitHashesFH (blastOutputFH, format(optional)) returns(list of hash(qname, hname, qlen, hlen, signif, bit hdesc, qcover, hcover, hstrand) );
sub hitHashesFH {
my $self = shift(@_);
my $fh=shift(@_);
my $format = 'blast';
my $temp = shift(@_);
if($temp) { $format = $temp; }
#print(STDERR "in topHit $filename\n");
  my $searchio = new Bio::SearchIO( -format => $format,
				    -fh   => $fh
                                  );
  my @retlist;
while( my $result = $searchio->next_result() ) {
  unless($result) { return();}
  my $qname=$result->query_name();
  my $qdesc=$result->query_description();
  my $qlen=$result->query_length();
  while (my $hit = $result->next_hit()) {
    unless($hit->next_hsp()) { last; }
    my $hname=$hit->name();
    my $num_hsps = $hit->num_hsps();
    my $hlen=$hit->length();
    my $frac_id = sprintf("%.3f", $hit->frac_identical());
    my $hdesc=$hit->description();
    my $signif=$hit->significance();
    my $laq=$hit->length_aln('query');
    my $qcover = sprintf("%.3f", $laq/$qlen);
    my $lah=$hit->length_aln('hit');
    my $hcover = sprintf("%.3f", $lah/$hlen);
    my $qstart = $hit->start('query');
    my $qend = $hit->end('query');
    my $hstart = $hit->start('hit');
    my $hend = $hit->end('hit');
    my $bitScore = $hit->bits();
    my $strand = $hit->strand('hit');
    my $cA = chr(01);
    $hdesc=~s/$cA/ /g;
    my %rethash = (qname => $qname, hname => $hname, qlen => $qlen, hlen => $hlen,
                   signif => $signif, bit => $bitScore, hdesc => $hdesc,
                   qcover => $qcover, hcover => $hcover, hstrand => $strand,
                   qcov => $qcover, hcov => $hcover,
                   qstart => $qstart, qdesc => $qdesc,
                   qend => $qend, hstart => $hstart, hend => $hend, alnlen => $laq,
                   fracid => $frac_id, numhsps => $num_hsps);
    push(@retlist, { %rethash });
#    return($qname, $hname, $signif, $qcover, $hcover, $frac_id, $hlen);
}
push(@retlist, "\/\/");
}
return(@retlist);
}
# }}}

# {{{ mkblastpdb (hash(inlist, outfile, title, format, faaname)) returns(nothing) 
# writes a blast database with the outfile as basename.
sub mkblastpdb {
  my $self = shift(@_);
  my %args = @_;
  my @files = @{$args{inlist}};
#print(join(" ", @files), "\n");
  unless(ref($args{inlist})) {
    croak qq(A reference to a list of filenames is required.);
  }
  unless($args{outfile}) {
    croak qq(An outfilename is required.);
  }
  my $faaWanted;
  if($args{faaname}) { $faaWanted = $args{faaname}; }

  my($fh, $fn)=tempfile($template, DIR => $tempdir, SUFFIX => '.faa');
  my $protCnt = 0;

  if($args{format}!~m/fasta/i) {
    my $seqout = Bio::SeqIO->new(-fh => $fh, -format => 'fasta');
    my %ids; # Need to check that ids are unique;
    foreach my $infile (@files) {
      my $seqio = Bio::SeqIO->new(-file => $infile);
      while(my $seqobj=$seqio->next_seq()) {
        foreach my $feature ($seqobj->all_SeqFeatures()) {
          if($feature->primary_tag() eq 'CDS') {
            my $annotation;
            my @temp;
            my @temp1;
            my $id;
            my @tags = $feature->get_all_tags();
            foreach my $tag (@tags) {
              if($tag=~m/product|gene|note/i) {
                push(@temp, $feature->get_tag_values($tag));
              }
              if($tag=~m/locus_tag|systematic_id/) {
                push(@temp1, $feature->get_tag_values($tag));
              }
            }
            unless(@temp1) {
              if($feature->has_tag("gene")) {
                push(@temp1, $feature->get_tag_values("gene"));
              }
            }
            $id = $temp1[0];
            if($id) { 
              if(exists($ids{$id})) { $id .= "_1"; }
              $ids{$id} += 1;
              $annotation = join(" ", @temp);
              my $aa_obj = _feat_translate($feature);
              $aa_obj->display_name($id);
              $aa_obj->description($annotation);
              $seqout->write_seq($aa_obj);
              $protCnt += 1;
            }
          }
        }
      }
    }
    close($fh);
  }

  else {
    my $seqout = Bio::SeqIO->new(-fh => $fh, -format => 'fasta');
    foreach my $infile (@files) {
      my $incomingIsTemp = 0;
      if($infile =~ m/\.gz$/) {
        my @temp = split(/\./, $infile);
        my $temp1 = $temp[-2];
        unless($temp1) { $temp1 = '.tmp'; }
        my($infh, $infn)=tempfile($template, DIR => $tempdir, SUFFIX => $temp1);
        unless(gunzip $infile => $infh, AutoClose => 1) {
          close($infh); unlink($infn);
          carp("gunzip failed for $infile $GunzipError\n");
          next;
        }
        $infile = $infn;
        $incomingIsTemp = 1;
      }
      my $seqio = Bio::SeqIO->new(-file => $infile);
      while(my $seqobj=$seqio->next_seq()) {
        $seqout->write_seq($seqobj);
        $protCnt += 1;
      }
      if($incomingIsTemp) { unlink($infile); }
    }
    close($fh);
  }

  if($protCnt) {
    my $xstr = qq(/usr/local/bin/makeblastdb -in $fn -dbtype prot);
    $xstr .= qq( -title "$args{title}" -out $args{outfile} -parse_seqids);
    my $mbdbout = qx($xstr);
    if($faaWanted) {
      copy($fn, $faaWanted);
    }
    unlink($fn);
    return(1);
  }
  else {
    carp(qq(Protein count zero\n));
    unlink($fn);
    return(0);
  }
}
# }}}

# {{{ internal sub _feat_translate
sub _feat_translate {
  my $feature=shift(@_);
  my $codon_start=1;
  if($feature->has_tag('codon_start')) {
      ($codon_start) = $feature->get_tag_values('codon_start');
      }
      my $offset=1;
      if($codon_start > 1) { $offset = $codon_start;}
      my $featobj=$feature->spliced_seq(-nosort => '1');
      my $aaobj=$featobj->translate(-offset => $offset, -complete => 1);
  return($aaobj);
}
# }}}

# {{{ blastpxml (hash(query, db, expect, outfh, threads, maxtargs, outfmt)).
sub blastpxml {
  my $self = shift(@_);
  my %args = @_;
  my $query = $args{query};
  my $ofh = $args{outfh};
  my $db = $args{db};
  my $evalue = $args{expect};
  unless ($evalue) { $evalue = $args{evalue} }
  unless ($evalue) { $evalue = 1e-4; }
  my $outfmt = 5;
  my $n_threads = 4;
  if($args{threads}) { $n_threads = $args{threads}; }
  my $maxtargs = 500;
  if($args{maxtargs}) { $maxtargs = $args{maxtargs}; }

# Thu 25 Jul 2013
# Changed -comp_based_stats to 2 in the blastp calls below
# because that is the default in blastp and this is how
# blastp calls were made when the ortholog tables were made.
# This option can influence which protein comes out as the
# top hit in a blastp run. This in turn can affect whether
# a pair are reciprocal best hits or not.

  my($fh1, $fn1)=tempfile($template, DIR => $tempdir, SUFFIX => '.xml');
  if(ref($query)) {
    my($fh, $fn)=tempfile($template, DIR => $tempdir, SUFFIX => '.faa');
    my $seqout = Bio::SeqIO->new(-fh => $fh, -format => 'fasta');
    $seqout->write_seq($query);
    close($fh);
    qx($blastbindir/blastp -max_target_seqs $maxtargs -num_threads $n_threads -outfmt $outfmt -query $fn -db $db -evalue $evalue -out $fn1 -comp_based_stats 2 -seg no);
    unlink($fn);
  }
  elsif(-e $query and -r $query) {
    qx($blastbindir/blastp -num_threads $n_threads -max_target_seqs $maxtargs -outfmt $outfmt -query $query -db $db -evalue $evalue -out $fn1 -comp_based_stats 2 -seg no);
  }
  close($fh1);
  open(BL, "<$fn1");
  while(<BL>) {
    print($ofh $_);
  }
  close(BL);
  close($ofh);
  unlink($fn1);
  return(1);
}
# }}}

# {{{ blastp (hash(query, db, expect, outfh, threads, naln, ndesc, outfmt)).
sub blastp {
  my $self = shift(@_);
  my %args = @_;
  my $query = $args{query};
  my $ofh = $args{outfh};
  my $db = $args{db};
  my $evalue = $args{expect};
  unless ($evalue) { $evalue = $args{evalue} }
  unless ($evalue) { $evalue = 1e-4; }
  my $outfmt = 0;
  if(exists($args{outfmt})) { $outfmt = $args{outfmt}; }
  my $ndesc = 500;
  if($args{ndesc}) { $ndesc = $args{ndesc}; }
  my $naln = 500;
  if($args{naln}) { $naln = $args{naln}; }
  my $n_threads = 4;
  if($args{threads}) { $n_threads = $args{threads}; }
  # print(STDERR "In blastp with $query\n");

# Thu 25 Jul 2013
# Changed -comp_based_stats to 2 in the blastp calls below
# because that is the default in blastp and this is how
# blastp calls were made when the ortholog tables were made.
# This option can influence which protein comes out as the
# top hit in a blastp run. This in turn can affect whether
# a pair are reciprocal best hits or not.

  my($fh1, $fn1)=tempfile($template, DIR => $tempdir, SUFFIX => '.blast');
  close($fh1);
  if(ref($query)) {
    # print(STDERR "Doing reference $query\n");
    my($fh, $fn)=tempfile($template, DIR => $tempdir, SUFFIX => '.faa');
    my $seqout = Bio::SeqIO->new(-fh => $fh, -format => 'fasta');
    $seqout->write_seq($query);
    close($fh);
    qx($blastbindir/blastp -num_descriptions $ndesc -num_threads $n_threads -num_alignments $naln -outfmt $outfmt -query $fn -db $db -evalue $evalue -out $fn1 -comp_based_stats 2 -seg no);
    unlink($fn);
  }
  elsif(-e $query and -r $query) {
    # print(STDERR "Doing file $query\n");
    qx($blastbindir/blastp -num_descriptions $ndesc -num_alignments $naln  -outfmt $outfmt -query $query -db $db -evalue $evalue -out $fn1 -comp_based_stats 2 -seg no);
  }
  else {
    print(STDERR "Something wrong with the query file $query\n");
  }
  open(BL, "<$fn1");
  while(<BL>) {
    print($ofh $_);
  }
  close(BL);
  close($ofh);
  unlink($fn1);
  return(1);
}
# }}}

# {{{ tblastn (hash(query, outfh, db, expect, outfmt, evalue, task, maxtargets)).
# Returns a filename which you must unlink.
sub tblastn {
 my $self = shift(@_);
 my %args = @_;
 my $query = $args{query};
 my $ofh = $args{outfh};
 my $db = $args{db};
 my $evalue = $args{expect};
 unless ($evalue) { $evalue = $args{evalue} }
 unless ($evalue) { $evalue = 1e-4; }
 my $outfmt = 0;
 if(exists($args{outfmt})) { $outfmt = $args{outfmt}; }
 my $task = $args{task};
 unless($task) { $task = 'tblastn'; }
 my $ntarg;
 if(exists($args{maxtargets})) {
 $ntarg = "-max_target_seqs $args{maxtargets} "; 
 }


 my($fh1, $fn1)=tempfile($template, DIR => $tempdir, SUFFIX => '.blast');
 if(ref($query)) {
   my($fh, $fn)=tempfile($template, DIR => $tempdir, SUFFIX => '.fna');
   my $seqout = Bio::SeqIO->new(-fh => $fh, -format => 'fasta');
   $seqout->write_seq($query);
   close($fh);
   qx/tblastn -outfmt $outfmt $ntarg -query $fn -db $db -evalue $evalue -out $fn1 -seg no/;
   unlink($fn);
 }
 elsif(-e $query and -r $query) {
   qx/tblastn -outfmt $outfmt $ntarg -query $query -db $db -evalue $evalue -out $fn1 -seg no/;
 }
 if($ofh) {
  open(BL, "<$fn1");
  while(<BL>) {
    print($ofh $_);
  }
  close(BL);
  close($ofh);
  unlink($fn1);
  return(1);
 }
 else {
   return($fn1);
 }

}
# }}}

# {{{ blastn (hash(query, db, task, expect, outfmt, filter)) returns two filenames.
sub blastn {
 my $self = shift(@_);
 my %args = @_;
 my $query = $args{query};
 my $db = $args{db};
 my $evalue = $args{expect};
 unless ($evalue) { $evalue = $args{evalue} }
 unless ($evalue) { $evalue = 1e-4; }
 my $outfmt = 0;
 if(exists($args{outfmt})) { $outfmt = $args{outfmt}; }
 my $task = $args{task};
 unless($task) { $task = 'blastn'; }

 my $filter = "no";
 if($args{filter}) {
   $filter = "yes";
 }

 if(ref($query)) {
   my($fh, $fn)=tempfile($template, DIR => $tempdir, SUFFIX => '.fna');
   my $seqout = Bio::SeqIO->new(-fh => $fh, -format => 'fasta');
   $seqout->write_seq($query);
   close($fh);
   my($fh1, $fn1)=tempfile($template, DIR => $tempdir, SUFFIX => '.blast');
   my($fh2, $fn2)=tempfile($template, DIR => $tempdir, SUFFIX => '.blerr');
   qx/blastn -num_alignments 500 -outfmt $outfmt -query $fn -db $db -task $task -evalue $evalue -dust $filter -out $fn1 2> $fn2/;
   unlink($fn);
   return($fn1, $fn2);
 }
 elsif(-e $query and -r $query) {
   my($fh1, $fn1)=tempfile($template, DIR => $tempdir, SUFFIX => '.blast');
   my($fh2, $fn2)=tempfile($template, DIR => $tempdir, SUFFIX => '.blerr');
   qx/blastn -num_alignments 500 -outfmt $outfmt -query $query -db $db -task $task -evalue $evalue -dust $filter -out $fn1 2> $fn2/;
   return($fn1, $fn2);
 }
}
# }}}

# {{{ longhit (hash (file, id)) returns( {qname, hname, qlen, hlen, signif, bit hdesc, qcover, hcover, hstrand } );
sub longhit {
my $self = shift(@_);
my %args = @_;
my $filename = $args{file};
my $idThresh = $args{id};
unless($filename and $idThresh) {
  croak ("longhit giving up: Filename: $filename idThresh: $idThresh");
}
#print(STDERR "$filename\n");
  my $searchio = new Bio::SearchIO( -format => 'blast',
				    -file   => $filename );

  my $result = $searchio->next_result();
  unless($result) { return();}
  my $qname=$result->query_name();
  my $qlen=$result->query_length();
  my @hits = $result->hits();
  unless(@hits) { return; }
  my @hitsByLen = sort _sortHitsByLength @hits;

  foreach my $hit (@hitsByLen) {
    my $hname=$hit->name();
    my $hlen=$hit->length();
    my $frac_id = sprintf("%.3f", $hit->frac_identical());
    if($frac_id < $idThresh) { next; }
    my $hdesc=$hit->description();
    my $signif=$hit->significance();
    my $laq=$hit->length_aln('query');
    my $qcover = sprintf("%.3f", $laq/$qlen);
    my $lah=$hit->length_aln('hit');
    my $hcover = sprintf("%.3f", $lah/$hlen);
    my $qstart = $hit->start('query');
    my $qend = $hit->end('query');
    my $hstart = $hit->start('hit');
    my $hend = $hit->end('hit');
    my $bitScore = $hit->bits();
    my $strand = $hit->strand('hit');
    my %rethash = (qname => $qname, hname => $hname, qlen => $qlen, hlen => $hlen,
                   signif => $signif, bit => $bitScore, hdesc => $hdesc,
                   qcover => $qcover, hcover => $hcover, hstrand => $strand,
                   frac_id => $frac_id, qstart => $qstart, qend => $qend,
                   hstart => $hstart, hend => $hend, alen => $laq);
    return(%rethash);
#    return($qname, $hname, $signif, $qcover, $hcover, $frac_id, $hlen);
  }
  return();
}
# }}}

# {{{ hspHashes (blastOutputFileName, format) returns(list of hashes(qname, hname, qlen, hlen, signif, bit hdesc, qcover, hcover, hstrand) );
sub hspHashes {
  my $self = shift(@_);
  my $filename=shift(@_);
  my $format = 'blast';
  my $temp = shift(@_);
  if($temp) { $format = $temp; }
  my $cA = chr(1);
#print(STDERR "in topHit $filename\n");
  my $searchio = Bio::SearchIO->new( -format => $format,
      -file   => $filename );
  my @retlist;
  while( my $result = $searchio->next_result() ) {
    unless($result) { return();}
    my $qname=$result->query_name();
    my $qdesc=$result->query_description();
    my $qlen=$result->query_length();
    while (my $hit = $result->next_hit()) {
      my $num_hsps = $hit->num_hsps();
      while (my $hsp = $hit->next_hsp()) {
        my $hname=$hit->name();
        my $hlen=$hit->length();
        my $frac_id = sprintf("%.3f", $hsp->frac_identical());
        my $frac_conserved = sprintf("%.3f", $hsp->frac_conserved());
        my $num_id = $hsp->num_identical();
        my $num_conserved = $hsp->num_conserved();
        my $hdesc=$hit->description();
        $hdesc =~ s/$cA/ /g;
        my $signif=$hsp->significance();
        my $laq=$hsp->length('query');
        my $lah=$hsp->length('hit');
        my $qcov = sprintf("%.3f", $laq/$qlen);
        my $hcov = sprintf("%.3f", $lah/$hlen);
        my $qstart = $hsp->start('query');
        my $qend = $hsp->end('query');
        my $hstart = $hsp->start('hit');
        my $hend = $hsp->end('hit');
        my $hframe = $hsp->frame('hit');
        my $bitScore = $hsp->bits();
        my $strand = $hsp->strand('hit');
        my %rethash = (qname => $qname, hname => $hname, qlen => $qlen, hlen => $hlen,
            signif => $signif, bit => $bitScore, qdesc => $qdesc, hdesc => $hdesc,
            hstrand => $strand, qstart => $qstart, hframe => $hframe,
            qend => $qend, hstart => $hstart, hend => $hend, alnlen => $laq, lah => $lah,
            fracid => $frac_id, fracsim => $frac_conserved, qcov => $qcov,
            qstr => $hsp->query_string(), numhsps => $num_hsps,
            numid => $num_id, num_id => $num_id, numconserved => $num_conserved,
            num_conserved => $num_conserved, hstr => $hsp->hit_string(), expect => $signif,
            homolstr => $hsp->homology_string(),
            hcov => $hcov);
        push(@retlist, {%rethash});
      }
      push(@retlist, '//');
    }
    push(@retlist, '////');
  }
  return(@retlist);
}
# }}}

# {{{ n_hspHashes (hash(filename, nhit, nhsp, nresult)) returns(list of hashes(qname, hname, qlen, hlen, signif, bit hdesc, qcover, hcover, hstrand) );
sub n_hspHashes {
  my $self = shift(@_);
  my %args = @_;
  my $filename=$args{filename};
  my $format = 'blast';
  if($args{format}) { $format = $args{format}; }
# print(STDERR "in n_hspHashes $filename\n");
  my $searchio = Bio::SearchIO->new(-format => $format,
      -file   => $filename );
  my @retlist;

  my $nresult = 0;
  my $nhit = 0;
  my $nhsp = 0;

  if($args{nresult}) { $nresult = $args{nresult}; }
  if($args{nhit}) { $nhit = $args{nhit}; }
  if($args{nhsp}) { $nhsp = $args{nhsp}; }

  my $resultC = 0;
  while( my $result = $searchio->next_result() ) {
    unless($result) { return();}
    my $qname=$result->query_name();
    my $qlen=$result->query_length();
    my $hitC = 0;
    while (my $hit = $result->next_hit()) {
      my $hspC = 0;
      while (my $hsp = $hit->next_hsp()) {
        my $hname=$hit->name();
        my $hlen=$hit->length();
        my $frac_id = sprintf("%.3f", $hsp->frac_identical());
        my $frac_conserved = sprintf("%.3f", $hsp->frac_conserved());
        my $hdesc=$hit->description();
        my $signif=$hsp->significance();
        my $laq=$hsp->length('query');
        my $lah=$hsp->length('hit');
        my $qcov = sprintf("%.3f", $laq/$qlen);
        my $hcov = sprintf("%.3f", $lah/$hlen);
        my $qstart = $hsp->start('query');
        my $qend = $hsp->end('query');
        my $hstart = $hsp->start('hit');
        my $hend = $hsp->end('hit');
        my $hframe = $hsp->frame('hit');
        my $bitScore = $hsp->bits();
        my $strand = $hsp->strand('hit');
        my %rethash = (qname => $qname, hname => $hname, qlen => $qlen, hlen => $hlen,
            signif => $signif, bit => $bitScore, hdesc => $hdesc,
            hstrand => $strand, qstart => $qstart, hframe => $hframe,
            qend => $qend, hstart => $hstart, hend => $hend, alnlen => $laq,
            fracid => $frac_id, fracsim => $frac_conserved, qcov => $qcov,
            hcov => $hcov);
        push(@retlist, {%rethash});
        $hspC += 1;
        # print(STDERR join("\t", $resultC, $hitC, $hspC), "\n");
        if($nhsp and $hspC >= $nhsp) { last; }
      }
      $hitC += 1;
      if($nhit and $hitC >= $nhit) { last; }
    }
    push(@retlist, "\/\/");
    $resultC += 1;
    if($nresult and $resultC >= $nresult) { last; }
  }
  return(@retlist);
}
# }}}

# {{{ topHSPs (blastOutputFileName) returns(list of hashes(qname, hname, qlen, hlen, signif, bit hdesc, qcover, hcover, hstrand) );
# Gives the top HSP only of each hit in each blast result in a file
sub topHSPs {
my $self = shift(@_);
my $filename=shift(@_);
my $format = 'blast';
my $temp = shift(@_);
if($temp) { $format = $temp; }
#print(STDERR "in topHit $filename\n");
  my $searchio = new Bio::SearchIO( -format => $format,
				    -file   => $filename );
  my @retlist;
while( my $result = $searchio->next_result() ) {
  unless($result) { return();}
  my $qname=$result->query_name();
  my $qlen=$result->query_length();
  while (my $hit = $result->next_hit()) {
    my $hsp = $hit->next_hsp();
    if($hsp) {
    my $hname=$hit->name();
    my $hlen=$hit->length();
    my $frac_id = sprintf("%.3f", $hsp->frac_identical());
    my $hdesc=$hit->description();
    my $signif=$hsp->significance();
    my $laq=$hsp->length('query');
    my $lah=$hsp->length('hit');
    my $qcov = sprintf("%.3f", $laq/$qlen);
    my $hcov = sprintf("%.3f", $lah/$hlen);
    my $qstart = $hsp->start('query');
    my $qend = $hsp->end('query');
    my $hstart = $hsp->start('hit');
    my $hend = $hsp->end('hit');
    my $hframe = $hsp->frame('hit');
    my $bitScore = $hsp->bits();
    my $strand = $hsp->strand('hit');
    my %rethash = (qname => $qname, hname => $hname, qlen => $qlen, hlen => $hlen,
                   signif => $signif, bit => $bitScore, hdesc => $hdesc,
                   hstrand => $strand, qstart => $qstart, hframe => $hframe,
                   qend => $qend, hstart => $hstart, hend => $hend, alnlen => $laq,
                   fracid => $frac_id, qcov => $qcov, hcov => $hcov);
    push(@retlist, {%rethash});
}
  }
}
return(@retlist);
}
# }}}

# {{{ topHitHSPstrings (blastOutputFileName) returns(a list of hashes);
# Gives the top HSP only of each hit in each blast result in a file
sub topHitHSPstrings {
my $self = shift(@_);
my $filename=shift(@_);
my $format = 'blast';
my $temp = shift(@_);
if($temp) { $format = $temp; }
#print(STDERR "in topHit $filename\n");
  my $searchio = new Bio::SearchIO( -format => $format,
				    -file   => $filename );
  my @retlist;
  my $result = $searchio->next_result();
  unless($result) { return();}
  my $qname=$result->query_name();
  my $qlen=$result->query_length();

  while(my $hit = $result->next_hit()) {
  my $hsp = $hit->next_hsp();
    my $hname=$hit->name();
    my $hlen=$hit->length();
    my $frac_id = sprintf("%.3f", $hsp->frac_identical());
    my $hdesc=$hit->description();
    my $signif=$hsp->significance();
    my $laq=$hsp->length('query');
    my $lah=$hsp->length('hit');
    my $qcov = sprintf("%.3f", $laq/$qlen);
    my $hcov = sprintf("%.3f", $lah/$hlen);
    my $qstart = $hsp->start('query');
    my $qend = $hsp->end('query');
    my $hstart = $hsp->start('hit');
    my $hend = $hsp->end('hit');
    my $hframe = $hsp->frame('hit');
    my $bitScore = $hsp->bits();
    my $strand = $hsp->strand('hit');
    my $qstr = $hsp->query_string();
    my $hitstr = $hsp->hit_string();
    my $homstr = $hsp->homology_string();
    my %rethash = (qname => $qname, hname => $hname, qlen => $qlen, hlen => $hlen,
                   signif => $signif, bit => $bitScore, hdesc => $hdesc,
                   hstrand => $strand, qstart => $qstart, hframe => $hframe,
                   qend => $qend, hstart => $hstart, hend => $hend, alnlen => $laq,
                   fracid => $frac_id, qcov => $qcov, hcov => $hcov, qstr => $qstr,
                   hitstr => $hitstr, homstr => $homstr);
    push(@retlist, {%rethash});
  }
  return(@retlist);
}
# }}}

# {{{ hitHashes (blastOutputFileName format(optional)) returns(list of hash(qname, hname, qlen, hlen, signif, bit hdesc, qcover, hcover, hstrand) );
sub hitHashes {
my $self = shift(@_);
my $filename=shift(@_);
my $format = 'blast';
my $temp = shift(@_);
if($temp) { $format = $temp; }
#print(STDERR "in topHit $filename\n");
  my $searchio = new Bio::SearchIO( -format => $format,
				    -file   => $filename
                                  );
  my @retlist;
while( my $result = $searchio->next_result() ) {
  unless($result) { return();}
  my $qname=$result->query_name();
  my $qdesc=$result->query_description();
  my $qlen=$result->query_length();
  while (my $hit = $result->next_hit()) {
    unless($hit->next_hsp()) { last; }
    my $hname=$hit->name();
    my $num_hsps = $hit->num_hsps();
    my $hlen=$hit->length();
    my $frac_id = sprintf("%.3f", $hit->frac_identical());
    my $hdesc=$hit->description();
    my $signif=$hit->significance();
    my $laq=$hit->length_aln('query');
    my $qcover = sprintf("%.3f", $laq/$qlen);
    my $lah=$hit->length_aln('hit');
    my $hcover = sprintf("%.3f", $lah/$hlen);
    my $qstart = $hit->start('query');
    my $qend = $hit->end('query');
    my $hstart = $hit->start('hit');
    my $hend = $hit->end('hit');
    my $bitScore = $hit->bits();
    my $strand = $hit->strand('hit');
    my $cA = chr(01);
    $hdesc=~s/$cA/ /g;
    my %rethash = (qname => $qname, hname => $hname, qlen => $qlen, hlen => $hlen,
                   signif => $signif, bit => $bitScore, hdesc => $hdesc,
                   qcover => $qcover, hcover => $hcover, hstrand => $strand,
                   qcov => $qcover, hcov => $hcover,
                   qstart => $qstart, qdesc => $qdesc,
                   qend => $qend, hstart => $hstart, hend => $hend, alnlen => $laq,
                   fracid => $frac_id, numhsps => $num_hsps);
    push(@retlist, { %rethash });
#    return($qname, $hname, $signif, $qcover, $hcover, $frac_id, $hlen);
}
push(@retlist, "\/\/");
}
return(@retlist);
}
# }}}

# {{{ topHitMulti (blastOutputFileName, ethresh)
# returns(list of hash(qname, hname, qlen, hlen, signif, bit hdesc, qcover, hcover, hstrand) );
sub topHitMulti {
my $self = shift(@_);
my $filename=shift(@_);
my $ethresh = shift(@_);
my $format = 'blast';
my $temp = shift(@_);
if($temp) { $format = $temp; }
#print(STDERR "in topHit $filename\n");
  my $searchio = new Bio::SearchIO( -format => $format,
				    -file   => $filename );
  my @retlist;
while( my $result = $searchio->next_result() ) {
  unless($result) { return();}
  my $qname=$result->query_name();
  my $qdesc=$result->query_description();
  my $qlen=$result->query_length();
  while(my $hit = $result->next_hit()) {
  if($hit) {
    my $num_hsps = $hit->num_hsps();
    my $hname=$hit->name();
    my $hlen=$hit->length();
    my $frac_id = sprintf("%.3f", $hit->frac_identical());
    my $hdesc=$hit->description();
    my $signif=$hit->significance();
    my $laq=$hit->length_aln('query');
    my $qcover = sprintf("%.3f", $laq/$qlen);
    my $lah=$hit->length_aln('hit');
    my $hcover = sprintf("%.3f", $lah/$hlen);
    my $qstart = $hit->start('query');
    my $qend = $hit->end('query');
    my $hstart = $hit->start('hit');
    my $hend = $hit->end('hit');
    my $bitScore = $hit->bits();
    my $strand = $hit->strand('hit');
    my $cA = chr(01);
    $hdesc=~s/$cA/ /g;
    my %rethash = (qname => $qname, hname => $hname, qlen => $qlen, hlen => $hlen,
                   signif => $signif, bit => $bitScore, qdesc => $qdesc, hdesc => $hdesc,
                   qcover => $qcover, hcover => $hcover, hstrand => $strand, qstart => $qstart,
                   qcov => $qcover, hcov => $hcover, numhsps => $num_hsps,
                   qend => $qend, hstart => $hstart, hend => $hend, alnlen => $laq,
                   fracid => $frac_id);
    if($ethresh) {
      if($signif <= $ethresh) {
        push(@retlist, { %rethash });
      }
    }
    else {
        push(@retlist, { %rethash });
    }
#    return($qname, $hname, $signif, $qcover, $hcover, $frac_id, $hlen);
  }
}
}
return(@retlist);
}
# }}}

# {{{ tophit (blastOutputFileName, format) returns( {qname, hname, qlen, hlen, signif, bit hdesc, qcover, hcover, hstrand } );
sub tophit {
my $self = shift(@_);
my $filename=shift(@_);
my $format = shift(@_);

unless($format) { $format = 'blast'; }

#print(STDERR "$filename\n");
  my $searchio = Bio::SearchIO->new( -format => $format,
				    -file   => $filename );
#print(STDERR "$searchio\n");

  my $result = $searchio->next_result();
#print(STDERR "$result\n");
  unless($result) { return();}
  my $qname=$result->query_name();
  my $qdesc=$result->query_description();
  my $qlen=$result->query_length();
  my $hit = $result->next_hit();
  if($hit) {
    my $hname=$hit->name();
    my $hlen=$hit->length();
    my $num_hsps = $hit->num_hsps();
    my $frac_id = sprintf("%.3f", $hit->frac_identical());
    my $frac_cons = sprintf("%.3f", $hit->frac_conserved());
    my $hdesc=$hit->description();
    my $signif=$hit->significance();
    my $laq=$hit->length_aln('query');
    my $qcover = sprintf("%.3f", $laq/$qlen);
    my $lah=$hit->length_aln('hit');
    my $hcover = sprintf("%.3f", $lah/$hlen);
    my $qstart = $hit->start('query');
    my $qend = $hit->end('query');
    my $hstart = $hit->start('hit');
    my $hend = $hit->end('hit');
    my $bitScore = $hit->bits();
    my $strand = $hit->strand('hit');
    my $cA = chr(01);
    $hdesc=~s/$cA/ /g;
    if($format =~ m/blastxml/i) {
      my @temp = split(/\s+/, $hdesc);
      #$hname = $temp[0];
      @temp = split(/\s+/, $qdesc);
      $qname = $temp[0];
    }
    my %rethash = (qname => $qname, hname => $hname, qlen => $qlen,
        hlen => $hlen, qdesc => $qdesc, hit => $hit,
        signif => $signif, bit => $bitScore, hdesc => $hdesc,
        qcover => $qcover, hcover => $hcover, hstrand => $strand,
        qcov => $qcover, hcov => $hcover, fraccons => $frac_cons,
        fracid => $frac_id, frac_cons => $frac_cons, numhsps => $num_hsps,
        frac_id => $frac_id, qstart => $qstart, qend => $qend,
        hstart => $hstart, hend => $hend);
    return(%rethash);
#    return($qname, $hname, $signif, $qcover, $hcover, $frac_id, $hlen);
  }
  else {
    return();
  }
}
# }}}

# {{{ sub AUTOLOAD

=head2 sub AUTOLOAD()

Does its trick whenever a non-existent sub is called. Note that
if a non-existent variable is asked for then I<undef> is returned
without any warning. Caller should check the value received.

Remember to put an C<our($AUTOLOAD)> in the file scope to keep
C<use strict;> happy.

=cut

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

# {{{ _sortHitsByLength {
sub _sortHitsByLength {
return($b->length() - $b->end('hit') <=> $a->length() - $a->end('hit'));
}
# }}}

# {{{ multiResultTopHits (blastOutputFileName) returns( a list of hashrefs {qname, hname, qlen, hlen, signif, bit hdesc, qcover, hcover, hstrand } );
sub multiResultTopHits {
my $self = shift(@_);
my $filename=shift(@_);
#print(STDERR "$filename\n");
  my $searchio = new Bio::SearchIO( -format => 'blast',
				    -file   => $filename );
  my @retlist;
  while(my $result = $searchio->next_result()) {
  unless($result) { return();}
  my $qname=$result->query_name();
  my $qlen=$result->query_length();
  my $hit = $result->next_hit();
  if($hit) {
    my $hname=$hit->name();
    my $hlen=$hit->length();
    my $frac_id = sprintf("%.3f", $hit->frac_identical());
    my $hdesc=$hit->description();
    my $signif=$hit->significance();
    my $laq=$hit->length_aln('query');
    my $qcover = sprintf("%.3f", $laq/$qlen);
    my $lah=$hit->length_aln('hit');
    my $hcover = sprintf("%.3f", $lah/$hlen);
    my $qstart = $hit->start('query');
    my $qend = $hit->end('query');
    my $hstart = $hit->start('hit');
    my $hend = $hit->end('hit');
    my $bitScore = $hit->bits();
    my $strand = $hit->strand('hit');
    my $cA = chr(01);
    $hdesc=~s/$cA/ /g;
    my %rethash = (qname => $qname, hname => $hname, qlen => $qlen, hlen => $hlen,
                   signif => $signif, bit => $bitScore, hdesc => $hdesc,
                   qcover => $qcover, hcover => $hcover, hstrand => $strand,
                   frac_id => $frac_id, qstart => $qstart, qend => $qend, hstart => $hstart,
                   hend => $hend, lah => $lah, laq => $laq);
    push(@retlist, {%rethash});
#    return($qname, $hname, $signif, $qcover, $hcover, $frac_id, $hlen);
  }
  else {
    return();
  }
  }
  return(@retlist);
}
# }}}

# {{{ sub rpsblast (%(file or seqobj, db, expect)) returns(Filename)
# Callee has to unlink the returned filename.
sub rpsblast {
  my $self = shift(@_);
  my %args = @_;
  my $db;
  if(exists($args{db})) {
    $db = $args{db};
  }
  else {
    $db = qw(/home/sco/ncbiftp/Cdd);
  }
  my $expect;
  if(exists($args{expect})) {
    $expect = $args{expect};
  }
  else {
    $expect = 1e-4;
  }

 
  my $rpsblastdb = $db;
  # print(STDERR join(", ", $db, $rpsblastdb, $expect, %args), "\n");

  if(exists($args{file})) {
  my($fh1, $fn1)=tempfile($template, DIR => $tempdir, SUFFIX => '.blast');
  close($fh1);
  `/usr/local/bin/rpsblast -query $args{file} -db $rpsblastdb -evalue $expect -out $fn1`;
  return($fn1);
  }

  elsif(exists($args{seqobj})) {
  my($fh, $fn)=tempfile($template, DIR => $tempdir, SUFFIX => '.faa');
  my $seqout = Bio::SeqIO->new(-fh => $fh, -format => 'fasta');
  $seqout->write_seq($args{seqobj});
  close($fh);
  my($fh1, $fn1)=tempfile($template, DIR => $tempdir, SUFFIX => '.blast');
  close($fh1);
  `/usr/local/bin/rpsblast -query $fn -db $rpsblastdb -evalue $expect -out $fn1`;
  unlink($fn);
  return($fn1);
  }
}
# }}}

# {{{ asciiLine hash(width, protLen, start, end) returns a string.
sub asciiLine {
my $self = shift(@_);
my %args = @_;
my $oneDash = $args{protLen} / $args{width};
my $asciiStart = int($args{start} / $oneDash);
my $asciiEnd = int($args{end} / $oneDash);
my $leadSpaces = " " x ($asciiStart - 1);
my $dashes = "-" x (($asciiEnd - $asciiStart) + 1);
my $retstr = $leadSpaces . $dashes;
return($retstr);
}
# }}}

# {{{ asciiLine2 hash(width, protLen, start, end) returns a string.
sub asciiLine2 {
my $self = shift(@_);
my %args = @_;

my $spaces = " " x $args{width};

my $oneDash = $args{protLen} / $args{width};
my $asciiStart = int($args{start} / $oneDash);
my $asciiEnd = int($args{end} / $oneDash);

my $dashes = "-" x (($asciiEnd - $asciiStart) + 1);
if(length($dashes) > $args{width}) {
my $temp = substr($dashes, 0, $args{width});
$dashes = $temp;
}

my $dlen = length($dashes);
substr($spaces, $asciiStart, $dlen, $dashes);
if(length($spaces) > $args{width}) {
my $temp = substr($spaces, 0, $args{width});
$spaces = $temp;
}

my $retstr = "|" . $spaces . "|";

return($retstr);
}
# }}}

# {{{ reciblastp (hash(query, refdb, db, biodb, expect)).
# Returns two hashrefs or one hashref and undef.
# query is a faa file with a single sequence or a protein seqobj.
# db is the subject blast database
# refdb is the reference blast database. i.e. of the organism from which
#       the query comes.
# biodb is the Bio::DB::Fasta object from which the hit id can be retrieved.
# expect is the evalue threshold
sub reciblastp {
  my $self = shift(@_);
  my %args = @_;
  my $query = $args{query};
  my $db = $args{db};
  my $biodb = $args{biodb};
  my $refdb = $args{refdb};
  my $evalue = $args{expect};
  unless ($evalue) { $evalue = 1; }
  my $outfmt = 0;

# Note the use of -comp_based_stats to 2 in the blastp calls below.

  my($fh1, $fn1)=tempfile($template, DIR => $tempdir, SUFFIX => '.blast');
  close($fh1);
  if(ref($query)) {
    my($fh, $fn)=tempfile($template, DIR => $tempdir, SUFFIX => '.faa');
    my $seqout = Bio::SeqIO->new(-fh => $fh, -format => 'fasta');
    $seqout->write_seq($query);
    close($fh);
    qx($blastbindir/blastp -outfmt $outfmt -query $fn -db $db -evalue $evalue -out $fn1 -comp_based_stats 2 -seg no);
    unlink($fn);
  }
  elsif(-e $query and -r $query) {
    qx($blastbindir/blastp -outfmt $outfmt -query $query -db $db -evalue $evalue -out $fn1 -comp_based_stats 2 -seg no);
  }
  my %forward = tophit($self, $fn1, "blast");
  unlink($fn1);

  my %reverse;
  if(%forward) {
  my $fhname = $forward{hname};
  # carp("Blast.pm: $fhname");
  my $revquery = $biodb->get_Seq_by_id($fhname);
    my($fh2, $fn2)=tempfile($template, DIR => $tempdir, SUFFIX => '.faa');
    my $seqout2 = Bio::SeqIO->new(-fh => $fh2, -format => 'fasta');
    $seqout2->write_seq($revquery);
    close($fh2);
  my($fh3, $fn3)=tempfile($template, DIR => $tempdir, SUFFIX => '.blast');
    qx($blastbindir/blastp -outfmt $outfmt -query $fn2 -db $refdb -evalue $evalue -out $fn3 -comp_based_stats 2 -seg no);
  %reverse = tophit($self, $fn3, "blast");
  unlink($fn2);
  unlink($fn3);
  }
  else { return(); }

  if(%forward and %reverse) {
  return(\%forward, \%reverse);
  }
  elsif(%forward) {
    return(\%forward, undef);
  }
  else {
    return();
  }
}
# }}}

# {{{ mf_reciblastp (hash(query, refdb, db, biodb, expect)).
# Returns two hashrefs or one hashref and undef.
# query is a faa file with a single sequence or a protein seqobj.
# db is the subject blast database
# refdb is the reference blast database. i.e. of the organism from which
#       the query comes.
# biodb is the Bio::DB::Fasta object from which the hit id can be retrieved.
# expect is the evalue threshold
sub mf_reciblastp {
  my $self = shift(@_);
  my %args = @_;
  my $query = $args{query};
  my $db = $args{db};
  my $biodb = $args{biodb};
  my $refdb = $args{refdb};
  my $evalue = $args{expect};
  unless ($evalue) { $evalue = 1; }
  my $outfmt = 0;

# Note the use of -comp_based_stats to 2 in the blastp calls below.

  my($fh1, $fn1)=tempfile($template, DIR => $tempdir, SUFFIX => '.blast');
  close($fh1);
  if(ref($query)) {
    my($fh, $fn)=tempfile($template, DIR => $tempdir, SUFFIX => '.faa');
    my $seqout = Bio::SeqIO->new(-fh => $fh, -format => 'fasta');
    $seqout->write_seq($query);
    close($fh);
    qx($blastbindir/blastp -outfmt $outfmt -query $fn -db $db -evalue $evalue -out $fn1 -comp_based_stats 2 -seg no);
    unlink($fn);
  }
  elsif(-e $query and -r $query) {
    qx($blastbindir/blastp -outfmt $outfmt -query $query -db $db -evalue $evalue -out $fn1 -comp_based_stats 2 -seg no);
  }
  my @forward = topHitMulti($self, $fn1, $evalue);
  unlink($fn1);
  my @reverse;
  for my $forward (@forward) {
    if(ref($forward)) {
      my %forward = %{$forward};
      my %reverse;
      if(%forward) {
        my $fhname = $forward{hname};
# carp("Blast.pm: $fhname");
        my $revquery = $biodb->get_Seq_by_id($fhname);
        my($fh2, $fn2)=tempfile($template, DIR => $tempdir, SUFFIX => '.faa');
        my $seqout2 = Bio::SeqIO->new(-fh => $fh2, -format => 'fasta');
        $seqout2->write_seq($revquery);
        close($fh2);
        my($fh3, $fn3)=tempfile($template, DIR => $tempdir, SUFFIX => '.blast');
        qx($blastbindir/blastp -outfmt $outfmt -query $fn2 -db $refdb -evalue $evalue -out $fn3 -comp_based_stats 2 -seg no);
        %reverse = tophit($self, $fn3, "blast");
        if(keys(%reverse)) { push(@reverse, \%reverse); }
        unlink($fn2);
        unlink($fn3);
      }
    }
  }
  if(@forward and @reverse) {
    return(\@forward, \@reverse);
  }
  elsif(@forward) {
    return(\@forward, []);
  }
  else {
    return([], []);
  }
}
# }}}


# {{{ reciblastpn (hash(query, refdb, db, biodb, expect)).
# Returns two hashrefs or one hashref and undef.
# query is a faa file with a single sequence or a protein seqobj. This organism
# is the reference.
#
# refdb is the reference blastp database. i.e. of the organism from which
#       the query comes.
#
# db is the subject blastn database
#
# hitobj is a listref of Bio::Seq nucleotide objects from which the subsequence of the hit
# can be retrieved (for reciprocal blasting).
# expect is the evalue threshold

sub reciblastpn {
  my $self = shift(@_);
  my %args = @_;
  my $query = $args{query};
  my $db = $args{db};
  my $hitbiodb = $args{biodb};
  my $refdb = $args{refdb};
  my $evalue = $args{expect};
  unless ($evalue) { $evalue = 1; }
  my $outfmt = 0;

# Note the use of -comp_based_stats to 2 in the blastp calls below.

  my($fh1, $fn1)=tempfile($template, DIR => $tempdir, SUFFIX => '.blast');
  close($fh1);
  if(ref($query)) {
    my($fh, $fn)=tempfile($template, DIR => $tempdir, SUFFIX => '.faa');
    my $seqout = Bio::SeqIO->new(-fh => $fh, -format => 'fasta');
    $seqout->write_seq($query);
    close($fh);
    qx($blastbindir/tblastn -outfmt $outfmt -query $fn -db $db -evalue $evalue -out $fn1 -comp_based_stats 2 -seg no);
    unlink($fn);
  }
  elsif(-e $query and -r $query) {
    qx($blastbindir/tblastn -outfmt $outfmt -query $query -db $db -evalue $evalue -out $fn1 -comp_based_stats 2 -seg no);
  }
    if($main::forblout) {
    copy($fn1, $main::forblout);
    }
  my @tophsps = topHSPs($self, $fn1, "blast");
  my %forward;
  if(ref($tophsps[0])) {
    %forward = %{$tophsps[0]};
  }
  unlink($fn1);

  my %reverse;
  if(%forward) {
    my $fhname = $forward{hname};
    my $sobj = $hitbiodb->get_Seq_by_id($forward{hname});
    my $revquery;
    if($forward{hstrand} == -1) {
      $revquery = $sobj->trunc($forward{hstart}, $forward{hend})->revcom();
    }
    else {
      $revquery = $sobj->trunc($forward{hstart}, $forward{hend});
    }
    my($fh2, $fn2)=tempfile($template, DIR => $tempdir, SUFFIX => '.fna');
    my $seqout2 = Bio::SeqIO->new(-fh => $fh2, -format => 'fasta');
    $seqout2->write_seq($revquery);
    close($fh2);
    if($main::revquery) {
    copy($fn2, $main::revquery);
    }
    my($fh3, $fn3)=tempfile($template, DIR => $tempdir, SUFFIX => '.blast');
    close($fh3);
    qx($blastbindir/blastx -outfmt $outfmt -query $fn2 -db $refdb -evalue $evalue -out $fn3 -comp_based_stats 2 -seg no);
    if($main::revblout) {
    copy($fn3, $main::revblout);
    }
    my @tophsps = topHSPs($self, $fn3, "blast");
    if(ref($tophsps[0])) {
      %reverse = %{$tophsps[0]};
    }
# %reverse = tophit($self, $fn3, "blast");
    unlink($fn2);
    unlink($fn3);
  }
  else { return(); }

  if(%forward and %reverse) {
    return(\%forward, \%reverse);
  }
  elsif(%forward) {
    return(\%forward, undef);
  }
  else {
    return();
  }
}
# }}}

# {{{ sub justASub
sub justASub {
my $self = shift(@_);
my %args = @_;
print(join(" ", %args), "\n");
return();
}
# }}} 

# {{{ oldblastp (hash(query, db, expect, outfh, threads, naln, ndesc)).
sub oldblastp {
  my $self = shift(@_);
  my %args = @_;
  my $blastbindir = qq(/usr/local/blast/bin);
  my $query = $args{query};
  my $ofh = $args{outfh};
  my $db = $args{db};
  my $evalue = $args{expect};
  unless ($evalue) { $evalue = $args{evalue} }
  unless ($evalue) { $evalue = 1e-4; }
  my $outfmt = 0;
  if(exists($args{outfmt})) { $outfmt = $args{outfmt}; }
  my $ndesc = 500;
  if($args{ndesc}) { $ndesc = $args{ndesc}; }
  my $naln = 500;
  if($args{naln}) { $naln = $args{naln}; }
  my $n_threads = 4;
  if($args{threads}) { $n_threads = $args{threads}; }


  my($fh1, $fn1)=tempfile($template, DIR => $tempdir, SUFFIX => '.blast');
  if(ref($query)) {
    my($fh, $fn)=tempfile($template, DIR => $tempdir, SUFFIX => '.faa');
    my $seqout = Bio::SeqIO->new(-fh => $fh, -format => 'fasta');
    $seqout->write_seq($query);
    close($fh);
    qx($blastbindir/blastall -p blastp -m $outfmt -i $fn -d $db -e $evalue -o $fn1 -C T -F F);
    unlink($fn);
  }
  elsif(-e $query and -r $query) {
    qx($blastbindir/blastall -p blastp -m $outfmt -i $query -d $db -e $evalue -o $fn1 -C T -F F);
  }
  close($fh1);
  open(BL, "<$fn1");
  while(<BL>) {
    print($ofh $_);
  }
  close(BL);
  close($ofh);
  unlink($fn1);
  return(1);
}
# }}}


# {{{ subroutines tablist, linelist, tabhash for printing lists and hashes.
# The E versions are for printing to STDERR.

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


return(1);

__END__

