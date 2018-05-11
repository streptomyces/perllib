# Avoid using this module. Try to use Bio::DB::Taxonomy and
# Bio::Taxon instead.
# I am leaving this here just in case we need it for something
# which becomes too difficult or too slow using the Bioperl
# modules.
package Sco::Taxonomy;
use 5.14.0;
use lib qw(/home/sco /home/sco/perllib);
use Sco::Common qw(tablist linelist tablistE linelistE tabhash tabhashE tabvals
                   tablistV tablistVE linelistV linelistVE tablistH linelistH);
our $AUTOLOAD;


my $namesFile = qq(/home/sco/seq/nt/newftp/taxonomy/names.dmp);
my $nodesFile = qq(/home/sco/seq/nt/newftp/taxonomy/nodes.dmp);
my %parent;
my %names;
my %taxid;

# {{{ sub new
sub new {
  my($class, $self);
  $class=shift(@_);
  my %argv = @_;
  $self={};
  foreach my $key (%argv) {
    $self->{$key} = $argv{$key}
  }
  bless($self, $class);
  mkParentHash();
  mkNamesHash();
  return($self);
}
# }}}

### more subs here ###

# {{{ sub getParent
sub getParent {
  my $self = shift(@_);
  my @ids = @_;
  my @retlist;

  for my $id (@ids) {
    if(exists($parent{$id})) {
      push(@retlist, $parent{$id});
    }
    else {
      push(@retlist, undef);
    }
  }
  return(@retlist);
}
# }}}

# {{{ sub names2taxids
sub names2taxids {
  my $self = shift(@_);
  my @names = @_;
  #my @names = map {$_ =~ s/_/ /g} @temp;

  my @retlist;

  for my $name (@names) {
    $name =~ s/_/ /g;
    if(exists($taxid{$name})) {
      push(@retlist, $taxid{$name});
    }
    else {
      push(@retlist, undef);
    }
  }
  return(@retlist);
}
# }}}

# {{{ sub taxids2names
sub taxids2names {
  my $self = shift(@_);
  my @tids = @_;
  #my @names = map {$_ =~ s/_/ /g} @temp;

  my @retlist;

  for my $tid (@tids) {
    if(exists($names{$tid})) {
      push(@retlist, $names{$tid});
    }
    else {
      push(@retlist, undef);
    }
  }
  return(@retlist);
}
# }}}



# {{{ sub taxtrace {
sub taxtrace {
  my $self = shift(@_);
  my $taxid = shift(@_);
  my @trace = ($taxid);;
  my $tempid = $taxid;
  while(1) {
    my $pid = $parent{$tempid};
    if(grep { $_ eq $pid } @trace) { last; }
    else {
      push(@trace, $pid);
      $tempid = $pid;
    }
  }
  return(@trace);
}
# }}}



# {{{ sub mkParentHash. Called from new();
# populates %parent which is declared as our.
sub mkParentHash {
  open(my $nofh, "<$nodesFile") or croak("Could not open $nodesFile");
  while(my $line = readline($nofh)) {
    chomp($line);
    if($line=~m/^\s*\#/ or $line=~m/^\s*$/) {next;}
    my @ll=split(/\t\|\t/, $line);
    my $taxid = $ll[0];
    my $pid = $ll[1];
    my $divid = $ll[4];
  $parent{$taxid} = $pid;
  }
  close($nofh);
}
# }}}

# {{{ sub mkNamesHash. Called from new();
# Populates %names and %taxid both of which are declared with our.
sub mkNamesHash {
  open(my $nmfh, "<$namesFile") or croak("Could not open $namesFile");
  while(my $line = readline($nmfh)) {
    chomp($line);
    if($line=~m/^\s*\#/ or $line=~m/^\s*$/) {next;}
    my @ll=split(/\t\|\t/, $line);
    my $taxid = $ll[0];
    my $name = $ll[1];
    my $uniqName = $ll[2];
    my $classOfName = $ll[3];
    if($classOfName =~ m/scientific\s+name/i) {
      $taxid{$name} = $taxid;
      $names{$taxid} = $name;
    }
  }
  close($nmfh);
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




return(1);
