package Sco::Souk;
use 5.14.0;
use Bio::SeqIO;
use Bio::Seq;
use DBI;
our $AUTOLOAD;
use lib qw(/home/sco /home/sco/perllib);
use Scoglobal;
our @ISA = qw(Scoglobal);

my $soukDataDir = '/home/nouser/souk/data';

my $dbname='souk';
my $dbhost='n108377.nbi.ac.uk';
my $handle=DBI->connect("DBI:Pg:dbname=$dbname;host=$dbhost", 'sco', 'tsinH4x');

my %preAccession = (
SCO  => "AL645882",
SAP1  => "AP005645",
SAV   => "BA000030",
SCP1  => "AL589148",
SCP2  => "AL645771",
SCLAV => "CM000913",
SGR   => "NC_010572",
SCAB  => "Ssc",
SCL4  => "CM000914",
SVEN  => "FR845719",
SVEN15  => "SVEN2015"
);


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

# {{{ newCDS (%(file)) returns(0 or 1)
sub newCDS {
my $self=shift(@_);
my %args=@_;
open(INFILE, "<$args{file}");
my %rec;
while(<INFILE>) {
my $line=$_;
chomp($line);
if($line=~m/^\s*\#/ or $line=~m/^\s*$/) {next;}
if($line=~m/^\/\//) {
my ($maxSerial) = $handle->selectrow_array("select max(serial) from features where accession = '$rec{accession}' and serial < 100000");
my $nextSerial = $maxSerial + 1;

print<<"INSERT";
insert into features
(id, pri_tag, accession, serial, start_pos, end_pos, strand, product, annotation)
values('$rec{id}', '$rec{pri_tag}', '$rec{accession}', $nextSerial, $rec{start_pos}, $rec{end_pos}, $rec{strand}, '$rec{product}', '$rec{annotation}');

INSERT

print<<"PROT";
insert into proteins (accession, serial, gene, aa_seq, nt_seq)
values ('$rec{accession}', $nextSerial, '$rec{gene}', '$rec{aa_seq}', '$rec{nt_seq}');

PROT
print("\n");
%rec=();
}
else {
my @llist=split(/\t/, $line);
$rec{$llist[0]} = $llist[1];
}
}

close(INFILE);
return(1);
}
# }}}

# {{{ featuresAt (%(pre, posn) ) returns(@( [serial, id]))
sub featuresAt {
my $self=shift(@_);
my %args=@_;
my $pre = $args{pre};
my $posn = $args{posn};
my @retlist;
my $qstr = "select id, serial, start_pos, end_pos from features where accession = '$preAccession{$pre}' and 
start_pos <= $posn and end_pos >= $posn order by start_pos";
my $stmt = $handle->prepare($qstr); $stmt->execute();
while(my $ar = $stmt->fetchrow_arrayref()) {
push(@retlist, [@{$ar}]);
}
return(@retlist);
}
# }}}


# {{{ proteins (hash (ids => listref) ) returns list of Bio::Seq objects)
sub proteins {
my $self=shift(@_);
my %args=@_;
my @retlist;
my @ids = @{$args{ids}};
foreach my $id (@ids) {
my $qstr = "select product, aa_seq from vprot where id = '$id'";
my ($product, $aa_seq) = $handle->selectrow_array($qstr);
my $seqobj = Bio::Seq->new(-seq => $aa_seq);
$seqobj->display_name($id); $seqobj->description($product);
push(@retlist, $seqobj);
}
return(@retlist);
}
# }}}


# {{{ featURL (id) returns(a string like http:// ...)
sub featURL {
my $self=shift(@_);
my %args=@_;
my $id = $args{id};
my $qstr = "select accession, id from features where id = '$id'";
my($accession, $fid) = $handle->selectrow_array($qstr);
my $url = "http://strepdb.streptomyces.org.uk/cgi-bin/dc3.pl?accession=$accession&name=$id";
return($url); # change this
}
# }}}

# {{{ somesub (args) returns(retvals)
sub somesub {
my $self=shift(@_);
my %args=@_;
return('something'); # change this
}
# }}}

# {{{ accser2annotation {accession, serial}
sub accser2annotation {
my $self=shift(@_);
my %args=@_;
my $qstr = qq/select annotation, product, id, gene from vprot where
accession = '$args{accession}' and serial = $args{serial}/;
my %rethash;
($rethash{annotation}, $rethash{product},
$rethash{id}, $rethash{gene}) = $handle->selectrow_array($qstr);
return(%rethash);
}
# }}}

# {{{ id2annotation {id}
sub id2annotation {
my $self=shift(@_);
my %args=@_;
my $qstr = qq/select annotation, product, id, gene from vprot where
id = '$args{id}'/;
my %rethash;
($rethash{annotation}, $rethash{product},
$rethash{id}, $rethash{gene}) = $handle->selectrow_array($qstr);
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

# {{{ DESTROY
sub DESTROY {
$handle->disconnect();
}
# }}}


return(1);

