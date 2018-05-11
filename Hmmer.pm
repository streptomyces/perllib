package Sco::Hmmer;
use 5.16.0;
use Carp;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::Seq;
use File::Temp qw(tempfile tempdir);
my $tempdir = qw(/home/sco/volatile);
my $template="HmmerXXXXX";

our($AUTOLOAD);

sub new {
        my($class, $self);
        $class=shift(@_);
        $self={};
        bless($self, $class);
        return($self);
}

### more subs go below ###

# {{{ hspHashes (blastOutputFileName, format) returns(list of hashes(qname, hname, qlen, hlen, signif, bit hdesc, qcover, hcover, hstrand) );
sub hspHashes {
  my $self = shift(@_);
  my $filename=shift(@_);
  my $format = 'hmmer';
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
    while (my $hit = $result->next_hit()) {
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
        # my $hcov = sprintf("%.3f", $lah/$hlen);
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
            qend => $qend, hstart => $hstart, hend => $hend, alnlen => $laq,
            fracid => $frac_id, fracsim => $frac_conserved, qcov => $qcov,
            hcov => "replaced");
        push(@retlist, {%rethash});
      }
    }
    push(@retlist, "\/\/");
  }
  return(@retlist);
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
my $format = 'hmmer';
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
  my $qlen=$result->query_length();
  while (my $hit = $result->next_hit()) {
#    unless($hit->next_hsp()) { last; }
    my $hname=$hit->name();
    my $hlen=$hit->length();
    my $frac_id = sprintf("%.3f", $hit->frac_identical());
    my $hdesc=$hit->description();
    my $signif=$hit->significance();
    my $laq=$hit->length_aln('query');
    my $qcover = sprintf("%.3f", $laq/$qlen);
    my $lah=$hit->length_aln('hit');
#    my $hcover = sprintf("%.3f", $lah/$hlen);
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
                   qcov => $qcover, hcov => "invalid", hstrand => $strand, qstart => $qstart,
                   qend => $qend, hstart => $hstart, hend => $hend, alnlen => $laq,
                   fracid => $frac_id);
    push(@retlist, { %rethash });
#    return($qname, $hname, $signif, $qcover, $hcover, $frac_id, $hlen);
}
push(@retlist, "\/\/");
}
return(@retlist);
}
# }}}

# {{{ topHitMulti (blastOutputFileName, ethresh) returns(list of hash(qname, hname, qlen, hlen, signif, bit hdesc, qcover, hcover, hstrand) );
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
  my $qlen=$result->query_length();
  while(my $hit = $result->next_hit()) {
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
                   qcover => $qcover, hcover => $hcover, hstrand => $strand, qstart => $qstart,
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
  my $searchio = new Bio::SearchIO( -format => $format,
				    -file   => $filename );

  my $result = $searchio->next_result();
  unless($result) { return();}
  my $qname=$result->query_name();
  my $qdesc=$result->query_description();
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
    if($format =~ m/blastxml/i) {
      my @temp = split(/\s+/, $hdesc);
      $hname = $temp[0];
      @temp = split(/\s+/, $qdesc);
      $qname = $temp[0];
    }
    my %rethash = (qname => $qname, hname => $hname, qlen => $qlen, hlen => $hlen,
                   signif => $signif, bit => $bitScore, hdesc => $hdesc,
                   qcover => $qcover, hcover => $hcover, hstrand => $strand,
                   fracid => $frac_id,
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





return(1);

__END__

