package Sco::Fasta;
use 5.14.0;
use Bio::SeqIO;
use Bio::Seq;
our $AUTOLOAD;
use lib qw(/home/sco /home/sco/perllib);
use parent qw(Sco::Global);
use Bio::DB::Fasta;
use File::Spec;
use File::Temp qw(tempfile tempdir);
my $tempdir = qw(/home/sco/volatile);
my $template="fastapmXXXXX";
use IO::Uncompress::Gunzip qw(gunzip $GunzipError) ;

my $soukDataDir = '/home/nouser/souk/data';
my $blastbindir = qq(/usr/local/bin);

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

# {{{ regionTranslate (hash(file, start, end, strand)) returns(proteinSeq)
# for file you can pass sco, sve sgr sav, ssc, scl and it will find the respective chromosomes.
sub regionTranslate {
my $self=shift(@_);
my %args=@_;
my $file;
if(-e $args{file}) {
  $file = $args{file};
}
else {
$file = &_strepNtFastaFile($args{file});
}

my $seqIO = Bio::SeqIO->new(-file => $file);
my $seqObj = $seqIO->next_seq();
my $subObj;
if($args{strand} == 1) {
$subObj = $seqObj->trunc($args{start}, $args{end});
}
elsif($args{strand} == -1) {
$subObj = $seqObj->trunc($args{start}, $args{end})->revcom();
}
my $aaObj = $subObj->translate();
my $aaSeq = $aaObj->seq();
return($aaSeq);
}
# }}}


# {{{ regionTranslateM (hash(file, seqid, start, end, strand)) returns(proteinSeq)
# for file you can pass sco, sve sgr sav, ssc, scl and it will find the respective chromosomes.
sub regionTranslateM {
my $self=shift(@_);
my %args=@_;
my $file;
if(-e $args{file}) {
  $file = $args{file};
}
else {
$file = &_strepNtFastaFile($args{file});
}

my $seqIO = Bio::SeqIO->new(-file => $file);
while( my $seqObj = $seqIO->next_seq() ) {
if($seqObj->display_name() eq $args{seqid}) {
my $subObj;
if($args{strand} == 1) {
$subObj = $seqObj->trunc($args{start}, $args{end});
}
elsif($args{strand} == -1) {
$subObj = $seqObj->trunc($args{start}, $args{end})->revcom();
}
my $aaObj = $subObj->translate();
my $aaSeq = $aaObj->seq();
return($aaSeq);
}
}
return();
}
# }}}

# {{{ sub breaklines (string) returns a hash
sub breaklines {
  my $self = shift(@_);
  my $seq = shift(@_);
  my @lines;
  while(my $line = substr($seq, 0, 60, "")) {
    push(@lines, $line);
  }
  return(@lines);
}
# }}}

# {{{ sub markovBackground (hash(file)) returns a hash
sub markovBackground {
  my $self = shift(@_);
  my %args = @_;
  my $infile = $args{file};
  print(STDERR "$infile\n");
  my $seqio=Bio::SeqIO->new(-file => $infile);

  my %pair_freqs;
  my $pair_cnt;
  my $done_cnt = 0;
  while(my $seqobj=$seqio->next_seq()) {
#print(STDERR ++$done_cnt, "   \r");
    my $seq = $seqobj->seq();
    my $pos=0;
    while(my $pair = lc(substr($seq, $pos, 2))) {
      if(length($pair) == 2) {
        $pair_freqs{$pair} += 1;
        $pair_cnt += 1;
      }
      $pos+=1;
    }
  }

  my %rethash;
  foreach my $key(keys %pair_freqs) {
    if($key=~m/[acgt][acgt]/) {
      $rethash{$key} = $pair_freqs{$key} / $pair_cnt; 
#print($key, " ", $pair_freqs{$key} / $pair_cnt, "\n");
    }
  }

  $seqio=Bio::SeqIO->new(-file => $infile);
  my %freqs;
  my $cnt;
  while(my $seqobj=$seqio->next_seq()) {
    my $seq = $seqobj->seq();
    my $pos=0;
    while(my $nt = lc(substr($seq, $pos, 1))) {
      $freqs{$nt} += 1;
      $cnt += 1;
      $pos+=1;
    }
  }

  foreach my $key(keys %freqs) {
    if($key=~m/[acgt]/) {
      $rethash{$key} = $freqs{$key} / $cnt; 
# print($key, " ", $freqs{$key} / $cnt, "\n");
    }
  }
  return(%rethash);
}
# }}}

# {{{ sub reSearch (hash(file, regexp)) returns a list of ranges.
sub reSearch {
  my $self = shift(@_);
  my %args = @_;
  my $infile = $args{file};
  my $regexp = $args{regexp};
  #print(STDERR "$infile\n");
  my $seqio=Bio::SeqIO->new(-file => $infile);

  my $done_cnt = 0;
  my @finds;
  while(my $seqobj=$seqio->next_seq()) {
    my $seq = $seqobj->seq();
    while($seq=~m/($regexp)/ig) {
      my $foundlen = length($1);
      my $posret = pos($seq);
      my $start = $posret - $foundlen;
      my $end = $posret;
      push(@finds, [$start, $end, $1]);
    }
  }
return(@finds);
}
# }}}

# {{{ subseq (hash(file, start, end, strand)) returns(ntSeq)
# for file you can pass sco, sve sgr sav, ssc, scl and it will find the respective chromosomes.
sub subseq {
my $self=shift(@_);
my %args=@_;
my $file;
if(-e $args{file}) {
  $file = $args{file};
}
else {
$file = &_strepNtFastaFile($args{file});
}

my $seqIO = Bio::SeqIO->new(-file => $file);
my $seqObj = $seqIO->next_seq();
my $subObj;
if($args{strand} == 1) {
$subObj = $seqObj->trunc($args{start}, $args{end});
}
elsif($args{strand} == -1) {
$subObj = $seqObj->trunc($args{start}, $args{end})->revcom();
}
my $ntSeq = $subObj->seq();
return($ntSeq);
}
# }}}

# {{{ somesub (args) returns(retvals)
sub somesub {
my $self=shift(@_);
my %args=@_;
return('something'); # change this
}
# }}}

# {{{ sub getId ( hash(file, id) ). Returns a single Bio::Seq object.
sub getId {
  my $self = shift(@_);
  my %args = @_;
  my $id = exists($args{id}) ? $args{id} : $args{name};
  my $seqio = Bio::SeqIO->new(-file => $args{file});
  while(my $seqobj = $seqio->next_seq()) {
    my $name = $seqobj->display_name();
    if($name eq $id) {
      return($seqobj);
    }
  }
return();
}
# }}}

# {{{sub lengthAndGC (hash(file), format)
sub lengthAndGC {
  my $self = shift(@_);
  my %args = @_;
  my @retlist;
  my $format = 'fasta';
  if($args{format}) { $format = $args{format}; }
  my $seqio = Bio::SeqIO->new(-file => $args{file}, -format => $format);
  while(my $seqobj = $seqio->next_seq()) {
    my $seqlen = $seqobj->length();
    my $gcFrac;
    if($seqlen == 0) {
      $gcFrac = 0;
    }
    else {
    $gcFrac = $self->gc_frac($seqobj->seq());
    }
    my $seqid = $seqobj->display_name();
    push(@retlist, [$seqid, $seqlen, $gcFrac]);
  }
  return(@retlist);
}
# }}}

# {{{ sub identifiers (hash(file)) returns a list of identifiers;
sub identifiers {
  my $self = shift(@_);
  my %args = @_;
  my @retlist;
  my $seqio = Bio::SeqIO->new(-file => $args{file});
  while(my $seqobj = $seqio->next_seq()) {
    my $seqid = $seqobj->display_name();
    push(@retlist, $seqid);
  }
  return(@retlist);
}
# }}}

# {{{ sub id_count (hash(file)) returns a hash (id => count);
sub id_count {
  my $self = shift(@_);
  my %args = @_;
  my %retlist;
  my $seqio = Bio::SeqIO->new(-file => $args{file});
  while(my $seqobj = $seqio->next_seq()) {
    my $seqid = $seqobj->display_name();
    $retlist{$seqid} += 1;
  }
  return(%retlist);
}
# }}}


# {{{ sub idsAndDescs (hash(file)) returns a hash(identifiers, descriptions);
sub idsAndDescs {
  my $self = shift(@_);
  my %args = @_;
  my %rethash;
  my $seqio = Bio::SeqIO->new(-file => $args{file});
  while(my $seqobj = $seqio->next_seq()) {
    my $seqid = $seqobj->display_name();
    my $desc = $seqobj->description();
    $rethash{$seqid} = $desc;
  }
  return(%rethash);
}
# }}}

# {{{ sub _strepNtFastaFile
sub _strepNtFastaFile {
my $name = lc(shift(@_));
my $dir = $soukDataDir . '/' . $name;
my $file = $name . '_chr';
return($dir . '/' . $file . '.fas');
}
# }}}

# {{{ sub fasta2blastpDB %(file, bldbname, title, faafile)
# returns %(bldbname, infile, faafile);
sub fasta2blastpDB {
my $self = shift(@_);
my %args = @_;
#_dumphash(%args);

  my $filename = $args{file};
  my $incomingIsTemp = 0;

  my($faafh, $faafn);
  if($filename =~ m/\.gz$/) {
    if($args{faafile}) {
      $faafn = $args{faafile};
      open($faafh, ">", $faafn);
    }
    else {
      ($faafh, $faafn)=tempfile($template, DIR => $tempdir, SUFFIX => '.faa');
    }

    unless(gunzip $filename => $faafh, AutoClose => 1) {
      close($faafh); unlink($faafn);
      die "gunzip failed: $filename $GunzipError\n";
    }
    $filename = $faafn;
    $incomingIsTemp = 1;
  }

my $mkbldbbin = File::Spec->catfile($blastbindir, "makeblastdb");
qx($mkbldbbin -in $filename -title "$args{title}" -dbtype prot -out $args{bldbname});
# print(STDERR qq($mkbldbbin -in $filename -title "$args{title}" -dbtype prot -out $args{bldbname}), "\n");
if($incomingIsTemp) {
  unless($args{faafile}) {
    unlink($filename);
  }
}
return(bldbname => $args{bldbname}, infile => $args{file}, faafile => $faafn);
}
# }}}

# {{{ sub fasta2blastnDB %(file, bldbname, title, fnafile)
# returns %(bldbname, infile, fnafile);
sub fasta2blastnDB {
my $self = shift(@_);
my %args = @_;
#_dumphash(%args);

  my $filename = $args{file};
  my $incomingIsTemp = 0;

  my($fnafh, $fnafn);
  if($filename =~ m/\.gz$/) {
    if($args{fnafile}) {
      $fnafn = $args{fnafile};
      open($fnafh, ">", $fnafn);
    }
    else {
      ($fnafh, $fnafn)=tempfile($template, DIR => $tempdir, SUFFIX => '.fna');
    }

    unless(gunzip $filename => $fnafh, AutoClose => 1) {
      close($fnafh); unlink($fnafn);
      die "gunzip failed: $filename $GunzipError\n";
    }
    $filename = $fnafn;
    $incomingIsTemp = 1;
  }

my $mkbldbbin = File::Spec->catfile($blastbindir, "makeblastdb");
my $xstr = qq($mkbldbbin -in $filename -title "$args{title}" -dbtype nucl -out $args{bldbname});
qx($xstr);
if($?) {
print(STDERR $xstr, "\n");
croak($xstr);
}
if($incomingIsTemp) {
  unless($args{fnafile}) {
    unlink($filename);
  }
}
return(bldbname => $args{bldbname}, infile => $args{file}, fnafile => $fnafn);
}
# }}}

# {{{ sub selectEntries. hash(infile, ofh, \@ids)
sub selectEntries {
my $self = shift(@_);
my %args = @_;
my $infile = $args{infile};
my $ofh = $args{ofh};
my $listref = $args{ids};

my $seqout = Bio::SeqIO->new(-fh => $ofh, -format => "fasta");

# tie my %fas, 'Bio::DB::Fasta', $infile;
my $biodb = Bio::DB::Fasta->new($infile);

my $outcnt = 0;
for my $id (@{$listref}) {
my $seq = $biodb->seq($id);
my $head = $biodb->header($id);
my $desc = $head;
$desc =~ s/^\w+\s+//;
my $seqobj = Bio::Seq->new(-seq => $seq, -display_id => $id);
$seqobj->description($desc);
if($seqobj) {
$seqout->write_seq($seqobj);
  $outcnt += 1;
}
}
return($outcnt);
}


# {{{ sub id2seqobj ( hash(file, id) ). Returns a single Bio::Seq object.
sub id2seqobj {
  my $self = shift(@_);
  my %args = @_;
  my $id = exists($args{id}) ? $args{id} : $args{name};
  my $seqio = Bio::SeqIO->new(-file => $args{file});
  while(my $seqobj = $seqio->next_seq()) {
    my $name = $seqobj->display_name();
    if($name eq $id) {
      return($seqobj);
    }
  }
return();
}
# }}}

# {{{ sub allProteins %(ifh, ofh, minlen, maxlen maxperstop) returns();
# Nucleotide positions are zero based.
sub allProteins {
  my $self = shift(@_);
  my %args = @_;
  my $ifh = $args{ifh};
  my $ofh = $args{ofh};
  my $minlen = 0;
  if($args{minlen}) { $minlen = $args{minlen}; }
  my $maxlen = 0;
  if($args{maxlen}) { $maxlen = $args{maxlen}; }
  my $maxperstop = 0;
  if($args{maxperstop}) { $maxperstop = $args{maxperstop}; }

  my $seqio=Bio::SeqIO->new(-fh => $ifh);
  my $seqout=Bio::SeqIO->new(-fh => $ofh, -format => 'fasta');

  my $forobj=$seqio->next_seq();
  my $revobj = $forobj->revcom();
  my @strands = ("forward", "reverse");

  for my $seqobj ($forobj, $revobj) {
    my $seqlen = $seqobj->length();
    my $strand = shift(@strands);

    my $seq=uc($seqobj->seq());
    my(@start_posns, @stop_posns);
    @start_posns=();
    @stop_posns=();

    foreach my $start_codon ('ATG', 'GTG', 'TTG') {
      my $pos=-1;
      while (($pos=index($seq, $start_codon, $pos)) > -1) {
        my $frame=$pos % 3;
        push(@start_posns, [$pos, $frame, $start_codon]);
        $pos+=1;
      }
    }

    foreach my $stop_codon ('TAA', 'TAG', 'TGA') {
      my $pos=-1;
      while(($pos=index($seq, $stop_codon, $pos)) > -1) {
        my $frame=$pos%3;
        push(@stop_posns, [$pos, $frame, $stop_codon]);
        $pos+=1;
      }
    }

    my @sort_starts=sort {$a->[0] <=> $b->[0]} (@start_posns);
    my @sort_stops=sort {$a->[0] <=> $b->[0]} (@stop_posns);


    for my $wfr (0, 1, 2) {
      my @pair = (undef, undef);
      for my $stop (@sort_stops) {
        if($stop->[1] != $wfr) { next; }
        shift(@pair);
        push(@pair, $stop);
        if((defined $pair[0]) and (defined $pair[1]) ) {
          my @orf = ($pair[0]->[0] + 3, $pair[1]->[0] - 1);
          my $orseq = substr($seq, $pair[0]->[0] + 3, ($pair[1]->[0] - 1) - ($pair[0]->[0] + 3) + 1);


          my @intStarts;
          for my $start (@sort_starts) {
            if($start->[1] != $wfr) { next; }
            if($start->[0] >= $orf[0] and $start->[0] < $orf[1]) {
              push(@intStarts, $start->[0]);
            }
            if($start->[0] > $orf[1]) { last; }
          }
#print(join("\t", @orf, @intStarts), "\n", $orseq, "\n");

          my @temp = @intStarts; @intStarts = ();

          my $cnt = 0;
          for my $el (@temp) {
            push(@intStarts, $el);
            $cnt += 1;
            if($maxperstop and $cnt >= $maxperstop) { last; }
          }

          for my $inst (@intStarts) {
            my $cslen = ($orf[1] - $inst) + 1;
            my $cs = substr($seq, $inst, $cslen);
#print("$cs\n");
            my $aaobj = _translate($cs);
#print($aaobj->seq(), "\n\n");
            my $protId;
            if($strand eq "forward") {
              $protId = join("_", $strand, $inst + 1, $orf[1] + 1, $wfr+1);
            }
            elsif($strand eq "reverse") {
              $protId = join("_", $strand, $seqlen - $inst, $seqlen - $orf[1], $wfr+1);
            }
            $aaobj->display_id($protId); 
            if($aaobj->length() > $minlen) {
              if($maxlen) {
                if($aaobj->length() <= $maxlen) {
                  $seqout->write_seq($aaobj);
                }
              }
              else {
                $seqout->write_seq($aaobj);
              }
            }
          }

        }

      }
    }
  }

}
# }}}

# {{{ sub allORFs
sub allORFs {
  my $self = shift(@_);
  my %args = @_;
  my $ifh = $args{ifh};
  my $ofh = $args{ofh};
  my $seqio=Bio::SeqIO->new(-fh => $ifh);
  my $seqout=Bio::SeqIO->new(-fh => $ofh, -format => 'fasta');
  while(my $seqobj = $seqio->next_seq()) {
    my $forwSeq = uc($seqobj->seq());
    my $revSeq = uc($seqobj->revcom()->seq());
    my @strands = qw(forward reverse);

    for my $seq ($forwSeq, $revSeq) {
      my $strand = shift(@strands);

      my @stops;
      for my $stop_codon ('TAA', 'TAG', 'TGA') {
        my $pos=-1;
        while(($pos=index($seq, $stop_codon, $pos)) > -1) {
          my $frame=$pos%3;
          push(@stops, {pos => $pos, frame => $frame, codon => $stop_codon});
          $pos+=1;
        }
      }

      my @sorted = sort {$a->{pos} <=> $b->{pos}} @stops;
      @stops = @sorted;

      FRAME: for my $frame (0, 1, 2) {
        my $fastaSerial = 0;
        my @pair = (undef, undef);
        for my $stop (@stops) {
          if($stop->{frame} == $frame) {
            my $discard = shift(@pair);
            push(@pair, $stop);
          if(ref($pair[0]) and ref($pair[1])) {
            my $protStart = $pair[0]->{pos} + 3;
            my $protEnd = $pair[1]->{pos} - 1;
            my $seqlen = $protEnd - $protStart + 1;
            if($seqlen < 15) { next; }
            my $orfSeq = substr($seq, $protStart, $seqlen);
            my $orfProtObj = _translate($orfSeq);
            $fastaSerial += 1;
            # if($fastaSerial >= 20) { next FRAME; }
            my $id = join("_", $strand, $frame, $fastaSerial);
            my $desc;
            if($strand eq "forward") {
            $desc = join(" ", $frame, $protStart, $protEnd);
            }
            elsif($strand eq "reverse") {
            $desc = join(" ", $frame, length($seq) - $protEnd - 1,
            length($seq) - $protStart - 1 );
            }
            $orfProtObj->display_id($id);
            $orfProtObj->description($desc);
            $seqout->write_seq($orfProtObj);
          }
          }
        }
      }
    }
  }
}
# }}}

# {{{ sub sizeFilter. Hash(infile, ofh, minlen);
sub sizeFilter {
  my $self = shift(@_);
  my %args = @_;
  my $infas = $args{infile};
  my $ofh = $args{ofh};
  my $minlen = $args{minlen};
  my $seqio = Bio::SeqIO->new(-file => $infas, -format => 'fasta');
  my $seqout = Bio::SeqIO->new(-fh => $ofh, -format => 'fasta');
  my $ocnt = 0;
  while(my $seqobj = $seqio->next_seq()) {
    if($seqobj->length() >= $minlen) {
      $seqout->write_seq($seqobj);
      $ocnt += 1;
    }
  }
  return($ocnt);
}
# }}}

# {{{ Internal Subroutines

# {{{ sub _translate (ntstring). Returns aaobj.
sub _translate {
my $seq = shift(@_);
my $obj = Bio::Seq->new(-seq => $seq);
my $aaobj = $obj->translate();
return($aaobj);
}
# }}}


# }}}


# {{{ sub allProteinsMax %(ifh, ofh, minlen, maxlen) returns();
# Nucleotide positions are zero based.
sub allProteinsMax {
  my $self = shift(@_);
  my %args = @_;
  my $ifh = $args{ifh};
  my $ofh = $args{ofh};
  my $minlen = 0;
  if($args{minlen}) { $minlen = $args{minlen}; }
  my $maxlen = 0;
  if($args{maxlen}) { $maxlen = $args{maxlen}; }

  my $seqio=Bio::SeqIO->new(-fh => $ifh);
  my $seqout=Bio::SeqIO->new(-fh => $ofh, -format => 'fasta');

  my $seqCnt = 0;
  while(my $forobj=$seqio->next_seq()) {
    $seqCnt += 1;
    my $idprefix = "fasentry" . sprintf("%03d", $seqCnt);
    
  my $revobj = $forobj->revcom();
  my @strands = ("forward", "reverse");

  for my $seqobj ($forobj, $revobj) {
    my $seqlen = $seqobj->length();
    my $strand = shift(@strands);

    my $seq=uc($seqobj->seq());
    my(@start_posns, @stop_posns);
    @start_posns=();
    @stop_posns=();

    foreach my $start_codon ('ATG', 'GTG', 'TTG') {
      my $pos=-1;
      while (($pos=index($seq, $start_codon, $pos)) > -1) {
        my $frame=$pos % 3;
        push(@start_posns, [$pos, $frame, $start_codon]);
        $pos+=1;
      }
    }

    foreach my $stop_codon ('TAA', 'TAG', 'TGA') {
      my $pos=-1;
      while(($pos=index($seq, $stop_codon, $pos)) > -1) {
        my $frame=$pos%3;
        push(@stop_posns, [$pos, $frame, $stop_codon]);
        $pos+=1;
      }
    }

    my @sort_starts=sort {$a->[0] <=> $b->[0]} (@start_posns);
    my @sort_stops=sort {$a->[0] <=> $b->[0]} (@stop_posns);


    for my $wfr (0, 1, 2) {
      my @pair = (undef, undef);
      for my $stop (@sort_stops) {
        if($stop->[1] != $wfr) { next; }
        shift(@pair);
        push(@pair, $stop);
        if((defined $pair[0]) and (defined $pair[1]) ) {
          my @orf = ($pair[0]->[0] + 3, $pair[1]->[0] - 1);
          my $orseq = substr($seq, $pair[0]->[0] + 3, ($pair[1]->[0] - 1) - ($pair[0]->[0] + 3) + 1);


          my @intStarts;
          for my $start (@sort_starts) {
            if($start->[1] != $wfr) { next; }
            if($start->[0] >= $orf[0] and $start->[0] < $orf[1]) {
              push(@intStarts, $start->[0]);
            }
            if($start->[0] > $orf[1]) { last; }
          }

          my $cnt = 0;
          for my $inst (@intStarts) {
            my $cslen = ($orf[1] - $inst) + 1;
            my $cs = substr($seq, $inst, $cslen);
#print("$cs\n");
            my $aaobj = _translate($cs);
#print($aaobj->seq(), "\n\n");
            my $protId;
            if($strand eq "forward") {
              $protId = join("_", $idprefix, $strand, $inst + 1, $orf[1] + 1, $wfr+1);
            }
            elsif($strand eq "reverse") {
              $protId = join("_", $idprefix, $strand, $seqlen - $inst, $seqlen - $orf[1], $wfr+1);
            }
            $aaobj->display_id($protId); 
            if($aaobj->length() > $minlen) {
              if($maxlen) {
                if($aaobj->length() <= $maxlen) {
                  $seqout->write_seq($aaobj);
                  $cnt += 1;
                }
              }
              else {
                $seqout->write_seq($aaobj);
                $cnt += 1;
              }
            }
            if($cnt) { last; }
          }
        }
      }
    }
  }
}
}
# }}}


sub _dumphash {
  my %inh = @_;
  for my $key (keys %inh) {
    print("$key\t$inh{$key}\n");
  }
}



return(1);




