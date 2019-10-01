package Sco::ChIP_Seq;   # change this to the file name
use 5.14.0;
our $AUTOLOAD;
use lib qw(/home/sco/perllib);

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

sub DESTROY {
my $self = shift(@_);
#$handle->disconnect();
return(1);
}


### more subs here ###

sub lir {
  my $self = shift(@_);
  my %args = @_;
  my $dbh = $args{dbh};
  my $left = $args{left};
  my $right = $args{right};
  my $pritag = $args{pritag};
  my $qstr1 = qq/select max(end_pos) from features where end_pos < $left/;
  $qstr1 .= qq/ and pritag = '$pritag'/;
  my ($lpos) = $dbh->selectrow_array($qstr1);
  my $qstr2 = qq/select * from features where end_pos == $lpos/;
  my $stmt = $dbh->prepare($qstr2);
  $stmt->execute();
  my $left_hr = $stmt->fetchrow_hashref();
  
  my $qstr3 = qq/select min(start_pos) from features where start_pos > $right/;
  $qstr3 .= qq/ and pritag = '$pritag'/;
  my ($rpos) = $dbh->selectrow_array($qstr3);
  my $qstr4 = qq/select * from features where start_pos == $rpos/;
  my $stmtt = $dbh->prepare($qstr4);
  $stmtt->execute();
  my $right_hr = $stmtt->fetchrow_hashref();

  my $qstr5 = qq/select start_pos, end_pos from features/;
  $qstr5 .= qq/ where (start_pos <= $right and end_pos >= $right) /;
  $qstr5 .= qq/ or (start_pos <= $left and end_pos >= $left) /;
  $qstr5 .= qq/ and pritag = '$pritag'/;
  # carp($qstr5);
  # print(STDERR $qstr5);
  ($rpos, $lpos) = $dbh->selectrow_array($qstr5);
  my $in_hr;
  if($rpos and $lpos) {
  my $qstr6 = qq/select * from features where start_pos == $rpos and end_pos == $lpos/;
  my $stmttt = $dbh->prepare($qstr6);
  $stmttt->execute();
  $in_hr = $stmttt->fetchrow_hashref();
  }
  return($left_hr, $in_hr, $right_hr);
}


# {{{ featHash (featFile) returns(%annos)
sub featHash {
my $self=shift(@_);
my %args=@_;
my %annos;
open(my $feat, "<", $args{featFile}) or croak("Could not open $args{featFile}");
my $lineCnt = 0;
while(<$feat>) {
my $line=$_;
chomp($line);
my @llist=split(/\t/, $line);
my $lt = $llist[1];
my $anno = $llist[5];
my $geneLen = ($llist[3] - $llist[2]) + 1;
$annos{$lt} = {anno => $anno, len => $geneLen};
$lineCnt += 1;
}
close($feat);
print(STDERR "$lineCnt features read\n");
return(%annos);
}
# }}}

# {{{ sub bed2featList (bedfile)
sub bed2featList {
my $self = shift(@_);
my %args = @_;
my @retlist;
open(BED, "<", $args{bedfile});
my $trackLine = readline(BED);
while(<BED>) {
my $line = $_;
chomp($line);
my @llist = split(/\t/, $line);
my $strand;
if($llist[5] eq '+') { $strand = 1 }
elsif($llist[5] eq '-') { $strand = -1 }

push(@retlist, [ $llist[3], $llist[1], $llist[2], $strand ]);
}
close(BED);
return(@retlist);
}


# {{{ somesub (args) returns(retvals)
sub somesub {
my $self=shift(@_);
my %args=@_;
return('something'); # change this
}
# }}}






return(1);

