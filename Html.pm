package Sco::Html;   # change this to the file name
use strict;
use File::Basename;
our($AUTOLOAD);
use Carp;
use Cwd qw(abs_path);

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


### more subs here ###



# {{{ sub listToRowNcol hash(lr, ncol).
sub listToRowNcol {
  my $self=shift(@_);
  my %arg = @_;
  my $lr = $arg{lr};
  my @ll=@{$lr};
  my $ncol = $arg{ncol};


  my $retval;
    $retval .= qq(<tr>\n);
    my $elCnt = 0;
    for my $el (@ll) {
      if($elCnt == $#ll) {
        my $span = $ncol - $elCnt;
      $retval .= qq(<td colspan=$span>$el</td>);
      }
      else {
      $retval .= qq(<td>$el</td>);
      }
      $elCnt += 1;
    }

    $retval .= qq(\n</tr>\n);
  return($retval);
}
# }}}

# {{{ ### sub listToRow. list.
sub listToRow {
  my $self=shift(@_);
  my @ll=@_;
  my $retval;
    $retval .= qq(<tr>\n);
    for my $el (@ll) {
      $retval .= qq(<td>$el</td>);
    }
    $retval .= qq(\n</tr>\n);
  return($retval); # change this
}
# }}}

# {{{ ### sub listToTableHead. list.
sub listToTableHead {
  my $self=shift(@_);
  my @ll=@_;
  my $retval;
    $retval .= qq(<tr>\n);
    for my $el (@ll) {
      $retval .= qq(<th>$el</th>);
    }
    $retval .= qq(\n</tr>\n);
  return($retval); # change this
}
# }}}

# {{{ ### sub tdfToHtml ###
sub tdfToHtml {
  my $self=shift(@_);
  my @args=@_;
  my $infile = $args[0];
  my $retval;
  open(INFILE, "<$infile") or croak("Could not open $infile");
  my $lineCnt = 0;
  while(my $line = readline(INFILE)) {
    chomp($line);
    if($line=~m/^\s*\#/ or $line=~m/^\s*$/) {next;}
    my @ll=split(/\t/, $line);
    $retval .= qq(<tr>\n);
    for my $el (@ll) {
      $retval .= qq(<td>$el</td>);
    }
    $retval .= qq(\n</tr>\n);
  }
  close(INFILE);
  return($retval); # change this
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

__END__

