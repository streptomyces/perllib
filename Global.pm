package Sco::Global;
use strict;
use File::Basename;
our($AUTOLOAD);
use Carp;
use Cwd qw(abs_path);

my %months = (
jan => 1,
feb => 2,
mar => 3,
apr => 4,
may => 5,
jun => 6,
jul => 7,
aug => 8,
sep => 9,
oct => 10,
nov => 11,
dec => 12
);

my %sorting_pars = ();


sub new {
        my($class, $self);
        $class=shift(@_);
        $self={};
        bless($self, $class);
        return($self);
}

### more subs go below ###

# {{{ sub charToNum. Hash(char, base). Returns a number.
sub charToNum {
  my $self = shift(@_);
  my %args = @_;
  my $char = lc($args{char});
  my $base;
  if($args{base}) { $base = $args{base}; }
  else { $base = 0; }
  
  my @lcalpha = ("a".."z");
  my %lhash;
  for my $ndx ($base..$base+25) {
    $lhash{$lcalpha[$ndx-$base]} = $ndx;
  }
  if(exists($lhash{$char})) {
    return($lhash{$char});
  }
  else {
    return(undef);
  }
}
# }}}

# {{{ sub randomizeList. Takes a list (not a listref), returns a list.
sub randomizeList {
  my $self = shift(@_);
  my @list = @_;
  my @retlist = ();
  while(@list) {
    my $randex = int(rand($#list + 1));
    push(@retlist, splice(@list, $randex, 1));
  }
  return(@retlist);
}
# }}}

# {{{ sub catFile
sub catFile {
my $self = shift(@_);
my %args = @_;
my @argl = @_;
my $file;
if(exists($args{file})) {
  $file = $args{file};
}
else {
$file = shift(@argl);
}

open(FILE, "<$file") or croak("Failed to open $file\n");
while(<FILE>) {
print($_);
}
return(0);
}
# }}}

# {{{ listCompare (listref, listref) returns(common_ref unique1_ref unique2_ref)
sub listCompare {
  my $self = shift(@_);
  my $lr1 = shift(@_);
  my $lr2 = shift(@_);
  my @common;
  my @unique1;
  my @unique2;
  foreach my $el (@{$lr1}) {
    if(grep {$_ eq $el} @{$lr2}) {
      push(@common, $el);
    }
    else {
      push(@unique1, $el);
    }
  }
  foreach my $el (@{$lr2}) {
    unless(grep {$_ eq $el} @{$lr1}) {
      push(@unique2, $el);
    }
  }
  return(\@common, \@unique1, \@unique2);
}
# }}}

# {{{ listCompareNumeric (listref, listref) returns(common_ref unique1_ref unique2_ref)
sub listCompareNumeric {
  my $self = shift(@_);
  my $lr1 = shift(@_);
  my $lr2 = shift(@_);
  my @common;
  my @unique1;
  my @unique2;
  foreach my $el (@{$lr1}) {
    if(grep {$_ == $el} @{$lr2}) {
      push(@common, $el);
    }
    else {
      push(@unique1, $el);
    }
  }
  foreach my $el (@{$lr2}) {
    unless(grep {$_ == $el} @{$lr1}) {
      push(@unique2, $el);
    }
  }
  return(\@common, \@unique1, \@unique2);
}
# }}}

# {{{ sub file2list (filename, colnum) returns a list
sub file2list {
my $self = shift(@_);
my $file = shift(@_);
my $colnum = shift(@_); # zero based
my $temp = shift(@_);
my $separator;
if($temp) {
  $separator = $temp;
}
else {
  $separator = 't';
}
#print(STDERR "separator is $separator\n");
my @retlist;
open(IN, "<$file");
while(<IN>) {
my $line=$_;
chomp($line);
if($line=~m/^\s*\#/ or $line=~m/^\s*$/) {next;}
my @llist;
if($separator=~m/^t/i) {
  @llist=split(/\t/, $line);
}
elsif($separator=~m/^s/i) {
  @llist=split(/\s+/, $line);
}
push(@retlist, $llist[$colnum]);
}
close(IN);
return(@retlist);
}

# }}}

# {{{ sub prog (indicator, line_end = "", timestamp_bool = 0) prints indicator to STDERR with the specified line_end.
sub prog {
  my $self = shift(@_);
  my $indicator = shift(@_);
  my $delim = shift(@_);
  my $tsbool = shift(@_);
  unless ($indicator) { $indicator = '. '; }
  my $endl;
  if($delim=~m/^n/i) { $endl = "\n";}
  elsif($delim=~m/^r/i) { $endl = "\r";}
  else { $endl = "";}
  my $ts = "";
  if($tsbool) {
    $ts = localtime();
    print(STDERR "$ts\t", $indicator, $endl);
  }
  else {
    print(STDERR $indicator, $endl);
  }
}
# }}}

# {{{ sub prerr (message) prints message to STDERR.
sub prerr {
my $self = shift(@_);
my $message = shift(@_);
chomp($message);
print(STDERR $message, "\n");
}
# }}}

# {{{ sub htmlend
sub htmlend {
my $timestamp=`date`;
my $abspath = abs_path();
my $retval = <<"HTMLEND";
<hr />
<p class="timestamp">
Generated : $timestamp by $0 in $abspath
</p>
</body>
</html>
HTMLEND
return($retval);
}
# }}}

# {{{ htmlbegin (hash(title, pagehead, stylefile)) returns(htmlString)
sub htmlbegin {
  my $self = shift(@_);
  my %args = @_;
  unless(exists($args{stylefile})) {
    $args{stylefile} = '/home/sco/perllib/Sco/style.css';
  }
my $style;
open(STYLE, "<", $args{stylefile}) or warn("Failed to open $args{stylefile}");
while(<STYLE>) {
$style .= $_;
}
close(STYLE);

my $retval = <<"HTMLBEGIN";
<!DOCTYPE "html">
<html lang="en">
<head>
<meta charset="utf-8" />
<title>$args{title}</title>
$style
</head>
<body>
<script language="javascript" type="text/javascript">
// function mopup (shod, tabId, divId, inht)
function mopup (shod, tabId, divId, inht) {
  if(shod == "show") {
    var divh = document.getElementById(divId);
    var inh = document.getElementById(tabId);
    inh.innerHTML = inht;
    //inh.style.height = 400;
    divh.style.visibility = "visible";
    //alert("Div height " + divht);
    //debugprint("height of div: " + divht);
  }
  else if(shod == "hide") {
    var divh = document.getElementById(divId);
    divh.style.visibility = "hidden";
  }
}
</script>
<div id = "mopup">
<table id = "pu">
<!-- just some hidden text.-->
</table>
</div>
<h4>$args{pagehead}</h4>
HTMLBEGIN
return($retval);
}
# }}}

# {{{ sub random_nt_seq (length = 10000, gcfration = 0.5) returns (ntseq)
sub random_nt_seq {                       ### change name here
my $self=shift(@_);
my $reqlen = shift(@_);
my $gcfrac = shift(@_);
unless($reqlen) {$reqlen = 10000;}
unless($gcfrac) {$gcfrac = 0.5;}

my $nG = int(($reqlen*$gcfrac)/2);
my $nC = $nG;
my $nA = int(($reqlen - ($nC + $nG))/2);
my $nT = $nA;
my @basebank;
my $ntstr = '';

for(my $cnt = 0; $cnt < $nG; $cnt+=1) {
  push(@basebank, 'G');
}
for(my $cnt = 0; $cnt < $nC; $cnt+=1) {
  push(@basebank, 'C');
}
for(my $cnt = 0; $cnt < $nA; $cnt+=1) {
  push(@basebank, 'A');
}
for(my $cnt = 0; $cnt < $nT; $cnt+=1) {
  push(@basebank, 'T');
}

while(@basebank) {
  my $random = int(rand(scalar(@basebank)));
  my $rant = splice(@basebank, $random, 1);
  $ntstr .= $rant;
}
return($ntstr);
}

# }}}

# {{{ gc_frac (ntseq_string) returns (GC_count, Total_count, GC_fraction)
sub gc_frac {
my $self = shift(@_);
my $string = shift(@_);
my $gc_count;
my $tot_count;
my $pos=0;
while(my $nt = lc(substr($string, $pos, 1))) {
if($nt eq 'g' or $nt eq 'c' or $nt eq 's') {
$gc_count+=1;
}
if($nt=~m/[acgtsw]/) {
  $tot_count += 1;
}
$pos+=1;
}
my $gcf = $gc_count / $tot_count;
#print(STDERR "$string\t$gc_count\t$tot_count\t$retval\n");
return($gc_count, $tot_count, $gcf);
}

# }}} ################

# {{{ sub nt_freqs (ntseq) returns a hash like (a => 20, c => 50)

=head2 Sub nt_freqs

Arguments: A string of nucleotides

Side effect: None

Returns: A hash

=cut

sub nt_freqs {
my $self = shift(@_);
my $string = shift(@_);
my $gc_count;
my $tot_count;
my $pos=0;
my %rethash;
while(my $nt = lc(substr($string, $pos, 1))) {
  $rethash{$nt} += 1;
$pos+=1;
}
return(%rethash);
}

# }}} 

# {{{ sub parse_ptt_file {
# Location   Strand   Length   PID   Gene   Synonym   Code   COG   Product
# Start  End   Strand   Length   PID   Gene   Synonym   Code   COG ccref   Product
#   0     1       2        3      4      5       6        7     8     9     10
sub parse_ptt_file {
  my $self = shift(@_);
  my $pttfile = shift(@_);
#  print(STDERR "$pttfile\n");
  my @ptts;
  open(PTT, "<$pttfile");
  while(<PTT>) {
    my $line = $_;
    chomp($line);
    unless($line=~m/^\d+\.\.\d+/) {next;}
    my @llist = split(/\t/, $line);
    my $location = $llist[0];
    my ($t1, $t2) = split(/\.\./, $location);
    my $start = $t1 < $t2 ? $t1 : $t2;
    my $end = $t1 < $t2 ? $t2 : $t1;
    my $strand = $llist[1] eq '-' ? -1 : 1;
    my ($cog, $ccstr);
    my $ccref;
    if($llist[7]=~m/COG/) {
    $cog = $llist[7];
    $cog=~s/[A-Z]+$//;
    $ccstr = $llist[7];
    $ccstr=~s/^COG\d+//;
    $ccref= [split("", $ccstr)]; 
    }
    push(@ptts, [$start, $end, $strand, @llist[2,3,4,5,6],$cog, $ccref, $llist[8]]); 
  }
  close(PTT);
  return(\@ptts);
}
# }}}

# {{{ sub file2hash (filename, key_column, value_column) returns a hashref 

=head2 Sub file2hash

Arguments: Filename, column to use as key, column to use as values

Side effect: None

Returns: A hashref

=cut

sub file2hash {
  my $self = shift(@_);
  my $infile = shift(@_);
  my $keycol = shift(@_);
  my @valcols = @_;
  my %rethash;
  open(INFILE, "<$infile");
  while(<INFILE>) {
    my $line=$_;
    chomp($line);
    if($line=~m/^\s*\#/ or $line=~m/^\s*$/) {next;}
    my @llist=split(/\t/, $line);
    my $val;
    foreach my $vk (@valcols) {
      $val .= $llist[$vk];
    }
    $rethash{$llist[$keycol]} = $val; 
  }
  close(INFILE);
  return(\%rethash);
}

# }}} ################

# {{{ file2err (filename) prints all lines to STDERR
sub file2err {
  my $self = shift(@_);
  my $infile = shift(@_);
  open(INFILE, "<$infile");
  while(<INFILE>) {
    print(STDERR $_);
  }
  close(INFILE);
}

# }}} ################

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

# {{{ tophit (blastOutputFileName) returns( {qname, hname, qlen, hlen, signif, bit hdesc, qcover, hcover, hstrand } );
sub tophit {
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
                   qcover => $qcover, hcover => $hcover, hstrand => $strand);
    return(%rethash);
#    return($qname, $hname, $signif, $qcover, $hcover, $frac_id, $hlen);
  }
  else {
    return();
  }
}
# }}}

# {{{ feat_translate (Bio::Feature) returns(Bio::Seq)
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

# {{{ sub dateStamp
sub dateStamp {
my $self = shift(@_);
my $lts = localtime();
my @ts = split(/\s+/, $lts);
#my 0214
my $retstr=join(" ", @ts[0,2,1,4]);
return($retstr);
}
# }}}

# {{{ sub monthToNum 
sub monthToNum {
my $self = shift(@_);
my $month = lc(shift(@_));
#print(STDERR "\n\n m2n  $month\t$months{$month}\t$months{dec}\t stuff\n\n");
return($months{$month});
}

# }}}

# {{{ html_sequence (hash (sequence, width, start, end, class)) returns a html string.
sub html_sequence {
my $self = shift(@_);
my %args = @_;
my $seq = $args{sequence};
my $pos = 0;
my $width = 60; if(exists($args{width})) { $width = $args{width}; }
my $class = qq(coding); if(exists($args{class})) { $class = $args{class}; }
my $histart;  if(exists($args{start})) { $histart = $args{start}; }
my $hiend;  if(exists($args{end})) { $hiend = $args{end}; }
if(defined($histart)) { $histart -= 1; }

my @llist;
while(my $subseq = substr($seq, $pos, $width, "")) {
  push(@llist, $subseq);
}
my $wint = int($histart/$width);
my $wpos = $histart % $width;
substr($llist[$wint], $wpos, 0, qq(<span class="$class">));
$wint = int($hiend/$width);
$wpos = $hiend % $width;
substr($llist[$wint], $wpos, 0, qq(</span>));

my $retstr = join("<br/>\n", @llist); 
return($retstr);
}
# }}}

# {{{ sub cgi_header {
sub cgi_header {
my $self = shift(@_);
return(qq(Content-type: text/html\n\n));
}
# }}}

# {{{ sub hammingDistance hash(string1, string2) returns an integer.
sub hammingDistance {
my $self = shift(@_);
my %args = @_; 
my $string1 = $args{string1};
my $string2 = $args{string2};

if(length($string1) != length($string2)) {
carp("Unequal string lengths ", length($string1), "  ", length($string2), "\n");
return(0);
}
my $lastdx = length($string1) - 1;
my $distance = 0;
foreach my $idx (0..$lastdx) {
if(substr($string1, $idx, 1) ne substr($string2, $idx, 1)) {
$distance += 1;
}
}
return($distance);
}

# }}}

# {{{ sub revcom hash(string) returns as string
sub revcom {
my $self = shift(@_);
my %args = @_;
my $idx = -1;
my $retString;
while(my $nt = uc(substr($args{string}, $idx, 1))) {
$idx -= 1;
if($nt eq 'A') { $retString .= 'T'; }
elsif($nt eq 'C') { $retString .= 'G'; }
elsif($nt eq 'G') { $retString .= 'C'; }
elsif($nt eq 'T') { $retString .= 'A'; }
}
return($retString);
}
# }}}

# {{{ sub fileSort (hash(file, sep, numeric, order, column)) returns a listref.
# sep defaults to tab. 
sub fileSort {
my $self = shift(@_);
my %args = @_;
#_printHash(%args);
open(INFILE, "<$args{file}");

my @retlist;
my @unsorted;

while(<INFILE>) {
my $line=$_;
chomp($line);
if($line=~m/^\s*\#/ or $line=~m/^\s*$/) {next;}

my @llist;
if($args{sep}=~m/^w|^s/i) {
@llist = split(/\s+/, $line);
}
else {
@llist = split(/\t/, $line);
}
push(@unsorted, [@llist]);
}
close(INFILE);

if($args{numeric}) {
$sorting_pars{numeric} = 1;
}
else {
  $sorting_pars{numeric} = 0;
}
$sorting_pars{order} = $args{order};
$sorting_pars{column} = $args{column};
my @sorted = sort _sorter @unsorted;

return(\@sorted);
}
# }}}

# {{{ sub _sorter
sub _sorter {
my %sp = %sorting_pars;
my ($alpha, $beta);
if($sp{order}=~m/^d/) {
$alpha = $b; $beta = $a;
} else { $alpha = $a; $beta = $b; }

if($sp{numeric}) {
return ($alpha->[$sp{column}] <=> $beta->[$sp{column}]);
}
else {
return ($alpha->[$sp{column}] cmp $beta->[$sp{column}]);
}
}
# }}}

# {{{ Internal sub _printHash
sub _printHash {
my %hash = @_;
foreach my $key (keys %hash) {
print(STDERR "$key\t$hash{$key}\n");
}
}
# }}}

# {{{ sub dsfn hash(pre, post, ext). Returns a string.
sub dsfn {
  my $self = shift(@_);
  my %args = @_;
  if($args{ext}) { $args{ext} =~ s/\.//g; }

  unless($args{pre}) { $args{pre} = "pre"; }

my @wdays = qw(Sun Mon Tue Wed Thu Fri Sat);
my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $pday = $wdays[$wday];
my $pdate = sprintf("%02d", $mday);
my $pmonth = sprintf("%02d", $mon);
my $pyear = $year + 1900;
my $phour = sprintf("%02d", $hour);
my $pmin = sprintf("%02d", $min);
my $psec = sprintf("%02d", $sec);

my $fn = $args{pre} . '_' . $pyear . $pmonth . $pdate . '_' . $phour .
$pmin . $psec;
if($args{post}) {
$fn .= '_' . $args{post}; 
}
if($args{ext}) {
$fn .= '.' . $args{ext}; 
}

return($fn);
}
# }}}

# {{{ sub errlogfn ( hash(unique) ) returns a hash(logfile errfile);
sub errlogfn {
my $self = shift(@_);
my %args = @_;
my($scriptname, $dir, $extension) = fileparse($0, qr/\.[^.]*/);
my @wdays = qw(Sun Mon Tue Wed Thu Fri Sat);
my @months = qw(Jan Feb Mar Apr May Jun Jul Aug Sep Oct Nov Dec);
my ($sec,$min,$hour,$mday,$mon,$year,$wday,$yday,$isdst) = localtime(time);
my $pday = $wdays[$wday];
my $pdate = sprintf("%02d", $mday);
my $pmonth = sprintf("%02d", $mon+1);
my $pyear = $year + 1900;
my $phour = sprintf("%02d", $hour);
my $pmin = sprintf("%02d", $min);
my $psec = sprintf("%02d", $sec);
my ($errorfile, $logfile);
if($args{unique}) {
$errorfile = $scriptname . '_' . $pyear . $pmonth . $pdate . '_' . $phour .
$pmin . $psec . ".err";
$logfile = $scriptname . '_' . $pyear . $pmonth . $pdate . '_' . $phour .
$pmin . $psec . ".log";
}
else {
$errorfile = $scriptname . ".err";
$logfile = $scriptname . ".log";
}
return(logfile => $logfile, errfile => $errorfile);
}
# }}}


return(1);

