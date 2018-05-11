package Sco::AWS;   # change this to the file name
use 5.14.0;
use utf8;
use Carp;
use JSON::XS;
use File::Basename;
our $AUTOLOAD;

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

# {{{ sub instances
sub instances {

my $cmd = qq(aws --output json ec2 describe-instances);
my $cmdres = qx($cmd);
my $decoded = decode_json($cmdres);

my @retlist;

for my $reserv (@{$decoded->{Reservations}}) {
for my $instance (@{$reserv->{Instances}}) {
my @instlist;
push(@instlist, [ qq/InstanceId/, $instance->{InstanceId} ]);
push(@instlist, [ qq/ImageId/, $instance->{ImageId} ]);
push(@instlist, [ qq/PublicIpAddress/, $instance->{PublicIpAddress} ]);
push(@instlist, [ qq/VpcId/, $instance->{VpcId} ]);
push(@instlist, [ qq/SubnetId/, $instance->{SubnetId} ]);
push(@instlist, [ qq/SecurityGroup/, $instance->{SecurityGroups}->[0]->{GroupId} ]);
push(@instlist, [ qq/AvailabilityZone/, $instance->{Placement}->{AvailabilityZone} ]);
push(@instlist, [ qq/InstanceType/, $instance->{InstanceType} ]);
push(@instlist, [ qq/PrivateIpAddress/, $instance->{PrivateIpAddress} ]);

my $tr = $instance->{BlockDeviceMappings};
my @bdm = @{$tr};
my (@temp, @temp1);
for my $dev (@bdm) {
push (@temp, $dev->{DeviceName});
push (@temp1, $dev->{Ebs}->{VolumeId});
}
push(@instlist, [ qq/BlockDevices/, join(", ", @temp) ]);
push(@instlist, [ qq/VolumeIds/, join(", ", @temp1) ]);


push(@retlist, [@instlist]);
}

}

return(@retlist);

}
# }}}



# {{{ somesub (args) returns(retvals)
sub somesub {
my $self=shift(@_);
my %args=@_;
return('something'); # change this
}
# }}}


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
# }}}


# {{{ sub DESTROY {
sub DESTROY {
my $self = shift(@_);
#$handle->disconnect();
return(1);
}
# }}}

return(1);
