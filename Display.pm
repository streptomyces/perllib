package Sco::Display;
use common::sense;
use Carp;
use Bio::SeqIO;
use Bio::SearchIO;
use Bio::Seq;
use File::Temp qw(tempfile tempdir);
my $tempdir = qw(/home/sco/volatile);
my $template="BlastXXXXX";

our($AUTOLOAD);
my $blastbindir = qq(/usr/local/blastplus/bin);

sub new {
        my($class, $self);
        $class=shift(@_);
        $self={};
        bless($self, $class);
        return($self);
}

### more subs go below ###

# {{{ tab2bed (hash(file, outfile)).
#track	name=Annotation	description="Features"	useScore=0	visibility=1	itemRgb=on	color=100,20,20
sub blastp {
  my $self = shift(@_);
  my %args = @_;
  my $query = $args{query};
  return(1);
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

# Bio::Search::Result::ResultI
# Bio::Search::Hit::HitI
# Bio::Search::Hit::GenericHit
# Bio::Search::HSP::HSPI
# Bio::Search::HSP::GenericHSP

package ScoBlast;   # change this to the file name
use common::sense;
our $AUTOLOAD;
use Carp;
use Bio::SeqIO;
use Bio::Seq;
use Bio::SearchIO;
use File::Temp qw(tempfile tempdir);
my $tempdir = qw(/home/sco/volatile);
my $template="ScoBlastXXXXX";

my $strepBlastDbDir = qw(/home/nouser/souk/blast_databases);
my $blastBinDir = qw(/usr/local/blastplus/bin);


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


# {{{ blastn (hash(organism, query, evalue)) returns(retvals)
sub blastn {
my $self=shift(@_);
my %args=@_;
my($fh, $fn)=tempfile($template, DIR => $tempdir, SUFFIX => '.fas');
my $seqout = Bio::SeqIO->new(-fh => $fh, -format => 'fasta');
$seqout->write_seq($args{query});
my $evalue = 10;
if($args{evalue}) {$evalue = $args{evalue};}


my($fh1, $fn1)=tempfile($template, DIR => $tempdir, SUFFIX => '.blast');
close($fh1);
my $blastbin = $blastBinDir . '/' . 'blastn';
my $blastdb = $strepBlastDbDir . '/' . "$args{organism}/$args{organism}";
my $filter;
my @syslist = ($blastbin, '-query', $fn, '-db', $blastdb, '-dust', 'no', '-out', $fn1,
               '-evalue', $evalue, '-task', 'blastn');
system(@syslist);
unlink($fn);
return($fn1);
}
# }}}


# {{{ topHit (blastOutputFileName) returns( hash(qname, hname, qlen, hlen, signif, bit hdesc, qcover, hcover, hstrand) );
sub topHit {
my $self = shift(@_);
my $filename=shift(@_);
my $format = 'blast';
my $temp = shift(@_);
if($temp) { $format = $temp; }
#print(STDERR "in topHit $filename\n");
  my $searchio = new Bio::SearchIO( -format => $format,
				    -file   => $filename );

  my $result = $searchio->next_result();
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
    my $qcover = $laq/$qlen;
    my $lah=$hit->length_aln('hit');
    my $hcover = $lah/$hlen;
    my $qstart = $hit->start('query');
    my $qend = $hit->end('query');
    my $hstart = $hit->start('hit');
    my $hend = $hit->end('hit');
    my $bitScore = $hit->bits();
    my $strand = $hit->strand('hit');
    my %rethash = (qname => $qname, hname => $hname, qlen => $qlen, hlen => $hlen,
                   signif => $signif, bit => $bitScore, hdesc => $hdesc,
                   qcover => $qcover, hcover => $hcover, hstrand => $strand, qstart => $qstart,
                   qend => $qend, hstart => $hstart, hend => $hend, alnlen => $laq,
                   fracid => $frac_id);
    return(%rethash);
#    return($qname, $hname, $signif, $qcover, $hcover, $frac_id, $hlen);
  }
  else {
    return();
  }
}
# }}}


# {{{ topHSP (blastOutputFileName) returns( hash(qname, hname, qlen, hlen, signif, bit hdesc, qcover, hcover, hstrand) );
sub topHSP {
my $self = shift(@_);
my $filename=shift(@_);
#print(STDERR "$filename\n");
  my $searchio = new Bio::SearchIO( -format => 'blast',
				    -file   => $filename );

  my $result = $searchio->next_result();
  unless($result) { return();}
  my $qname=$result->query_name();
  my $qlen=$result->query_length();
  my $hit = $result->next_hit();
  if($hit) {
  my $hsp = $hit->next_hsp();
  if($hsp) {
    my $hname=$hit->name();
    my $hlen=$hit->length();
    my $frac_id = sprintf("%.3f", $hsp->frac_identical());
    my $hdesc=$hit->description();
    my $signif=$hsp->significance();
    my $laq=$hsp->length('query');
    my $lah=$hsp->length('hit');
    my $qstart = $hsp->start('query');
    my $qend = $hsp->end('query');
    my $hstart = $hsp->start('hit');
    my $hend = $hsp->end('hit');
    my $bitScore = $hsp->bits();
    my $strand = $hsp->strand('hit');
    my %rethash = (qname => $qname, hname => $hname, qlen => $qlen, hlen => $hlen,
                   signif => $signif, bit => $bitScore, hdesc => $hdesc,
                   hstrand => $strand, qstart => $qstart,
                   qend => $qend, hstart => $hstart, hend => $hend, alnlen => $laq,
                   fracid => $frac_id);
    return(%rethash);
#    return($qname, $hname, $signif, $qcover, $hcover, $frac_id, $hlen);
  }
  else {
    return();
  }
  }
  else {
    return();
  }
}
# }}}




return(1);
