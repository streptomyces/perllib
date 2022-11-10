=head1 Name

Sco::Protein

=head1 Description

This module contains some functions for proteins sequences.

=head1 Example

=head1 See also

The Perl documentation

=head1 History of changes

=over 2

=item Sometime in 2007: First written.

=item 1 May 2008: Added sub ex_coef.

=back

=head1 Author

Govind Chandra E<lt>govind.chandra@bbsrc.ac.ukE<gt>

=cut

package Sco::Protein;   # change this to the file name
use strict;
use Bio::Seq;

sub new {
  my($class, $self, %argv);
  ($class, %argv)=@_;
  unless($argv{aa_seq}) {
    print(STDERR "AA sequence (aa_seq) is required.\nReturning zero");
    return(0);
  }
  $self={
    'aa_seq' => $argv{aa_seq},
    'nt_seq' => $argv{nt_seq}
  };
  bless($self, $class);
  return($self);
}


### more subs go below ###

# {{{ sub ex_coef ###

sub ex_coef {
my $self=shift(@_);
my $prot=$self->{aa_seq};
my $seq;
if(ref($prot) eq "Bio::Seq") {
  $seq=$prot->seq();
}
elsif(ref($prot) eq "") {
  $seq=$prot;
}
else {
  my $type = ref($prot);
  die <<"STOPPING";
Sco::Protein sub ex_coef
Argument can either be a Bio::Seq object or a scalar.
You have provided $type.
STOPPING
}

#EC = (#Trp W)(5,500) + (#Tyr Y)(1,490) + (#cystine C)(l25).
my %counts=('W' => 0, 'Y' => 0, 'C' => 0);

my $pos=0;
while(my $aa = uc(substr($seq, $pos, 1))) {
  if($aa eq 'W'
      or $aa eq 'Y'
      or $aa eq 'C'
    ) {
    $counts{$aa}+=1;
  }
  $pos+=1;
}

#EC = (#Trp W)(5,500) + (#Tyr Y)(1,490) + (#cystine C)(l25).

my $ec = $counts{W}*5500 + $counts{Y}*1490 + $counts{C}*125;
return($ec, \%counts);
}

# }}}

# {{{ ### sub isoelectric ###

sub isoelectric {
  my $self=shift(@_);
  my($string, $aa, %aacount, $pos, $key, $pH, %pk);
  my $prot=$self->{aa_seq};
  if(ref($prot) eq "Bio::Seq") {
    $string=$prot->seq();
  }
  elsif(ref($prot) eq "") {
    $string=$prot;
  }
  %aacount=(
      'C'=>0,
      'D'=>0,
      'E'=>0,
      'Y'=>0,
      'H'=>0,
      'K'=>0,
      'R'=>0
      );

  %pk=(
      'nter'=>8,
      'C'=>8.5,
      'D'=>4.4,
      'E'=>4.4,
      'Y'=>10,
      'H'=>6.5,
      'K'=>10,
      'R'=>12,
      'cter'=>3.1
      );

  $pos=0;
  while($aa=uc(substr($string,$pos,1))){
    if($aa eq 'C'){
      $aacount{'C'}+=1;
    }
    elsif($aa eq 'D'){
      $aacount{'D'}+=1;
    }
    elsif($aa eq 'E'){
      $aacount{'E'}+=1;
    }
    elsif($aa eq 'Y'){
      $aacount{'Y'}+=1;
    }
    elsif($aa eq 'H'){
      $aacount{'H'}+=1;
    }
    elsif($aa eq 'K'){
      $aacount{'K'}+=1;
    }
    elsif($aa eq 'R'){
      $aacount{'R'}+=1;
    }
    $pos+=1;
  } # while ends here

  foreach $key (sort(keys(%aacount))){
#    print($key,'   ---->   ',$aacount{$key},"\n");
  }

  $pH=7.0;
  my $move=$pH/2;
  while(1){
    my $cr_c=10**($pH-$pk{'C'});
    my $cr_d=10**($pH-$pk{'D'});
    my $cr_e=10**($pH-$pk{'E'});
    my $cr_y=10**($pH-$pk{'Y'});
    my $cr_cter=10**($pH-$pk{'cter'});
    my $cr_h=10**($pk{'H'}-$pH);
    my $cr_k=10**($pk{'K'}-$pH);
    my $cr_r=10**($pk{'R'}-$pH);
    my $cr_nter=10**($pk{'nter'}-$pH);


    my $pch_c=($cr_c/($cr_c+1))*$aacount{'C'};
    my $pch_d=($cr_d/($cr_d+1))*$aacount{'D'};
    my $pch_e=($cr_e/($cr_e+1))*$aacount{'E'};
    my $pch_y=($cr_y/($cr_y+1))*$aacount{'Y'};
    my $pch_cter=($cr_cter/($cr_cter+1));
    my $pch_h=($cr_h/($cr_h+1))*$aacount{'H'};
    my $pch_k=($cr_k/($cr_k+1))*$aacount{'K'};
    my $pch_r=($cr_r/($cr_r+1))*$aacount{'R'};
    my $pch_nter=($cr_nter/($cr_nter+1));

    my $charge=$pch_nter+$pch_h+$pch_k+$pch_r-$pch_cter-$pch_c-$pch_d-$pch_e-$pch_y;

#    print("At pH $pH charge is $charge\n");

    if(abs($charge) < 0.001){
#      print("At pH $pH charge is $charge\n");
      return $pH;
    }
    elsif($charge < 0){
      $pH-=$move;
    }
    elsif($charge > 0){
      $pH+=$move;
    }
    $move=$move/2;
  } # while 1 ends here
} # sub isoelectric ends here

# }}} ###


# {{{ ### sub molwt ###

sub molwt {
  my $self=shift(@_);
  my(%mwhash, $mw, $string, $pos, $aa);
  my $prot=$self->{aa_seq};
  if(ref($prot) eq "Bio::Seq") {
    $string=$prot->seq();
  }
  elsif(ref($prot) eq "") {
    $string=$prot;
  }
  %mwhash=(
      'A'=>'71.079',
      'R'=>'156.188',
      'N'=>'114.104',
      'D'=>'115.089',
      'C'=>'103.144',
      'Q'=>'128.131',
      'E'=>'129.116',
      'G'=>'57.052',
      'H'=>'137.142',
      'I'=>'113.160',
      'L'=>'113.160',
      'K'=>'128.174',
      'M'=>'131.198',
      'F'=>'147.177',
      'P'=>'97.117',
      'S'=>'87.078',
      'T'=>'101.105',
      'W'=>'186.213',
      'Y'=>'163.170',
      'V'=>'99.113'
        );

  $pos=0;
  $mw=18;
  while($aa=uc(substr($string,$pos,1))){
    $mw+=$mwhash{$aa};
    $pos+=1;
  }
#print("molwt is $mw\n");
  return int($mw);
} # sub molwt ends here

# }}} ###

# {{{ sub inframe_tta

=head2 Sub inframe_tta

Returns a list of positions (integers) where an inframe TTA codon is found
in the nucleotide sequence of C<$self>.

=head3 Arguments

None other than $self which is a Micro::Protein object.

=head3 Returns

Returns a list of positions (integers) where an inframe TTA codon is found
in the nucleotide sequence of C<$self>.

=cut

sub inframe_tta {
  my $self=shift(@_);

  my($seq);
  my $cds=$self->{nt_seq};
  if(ref($cds) eq "Bio::Seq") {
    $seq=$cds->seq();
  }
  elsif(ref($cds) eq "") {
    $seq=$cds;
  }
  unless ($seq) {
    print(STDERR "Nucleotide sequence is required\n");
    return();
  }
  $seq=~s/\s+//g;
  my @posns=();
  while($seq=~m/tta/gi) {
    my $pos=pos($seq)-3;
    if($pos%3 == 0) {
      push(@posns, $pos+1);
    }
  }
  if(@posns) {
    return(@posns);
  }
  else{return();}
}
# }}}

return(1);


