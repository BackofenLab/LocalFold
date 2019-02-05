#!/usr/bin/perl
use strict 'vars';
use warnings;
use Getopt::Long;
use Pod::Usage;
use List::Util qw/ min max /;
use Cwd qw/ getcwd abs_path /;
use File::Temp qw/ tempdir /;
use File::Basename;

=head1 NAME

localfold.pl

=head1 SYNOPSIS

localfold.pl -seqfile=FASTAFILE

LocalFold is a modification of the RNAplfold local folding algorithm.
LocalFold only considers base-pairs centred within a window such that
(i,j) are not located within the first and last positions of the window.

Options:

    -seqfile        fasta file
    -W              set window size to W (default: 200)
    -L              compute only basepairs with maximal span L (default: 150)
    -u              compute probabilities for regions of this length (default: 1)
                    are unpaired (default: 1)
    -skipbordernts  do not use the n first and last nucleotides of a window (default: 10)
    -cutoff         report only probabilities > cutoff (default: 0)
    -T              rescale energy parameters to a temperature (default 37C)
    -noLP           run RNAplfold with option -noLP (no lonely base pairs)
    -P              run RNAplfold with option -P (read energy parameters from paramfile)
    -nodot          do not compute dotplots
    -noacc          do not compute accessibilities
    -debug          enable debug output
    -help           brief help message
    -man            full documentation

=head1 DESCRIPTION

=cut

# variables for command line options
my $help;
my $man;
my $nodot;
my $noacc;
my $seqfile;
my $window;
my $maxspan;
my $u;
my $skipbordernts;
my $cutoff;
my $debug;
my $noLP;		# RNAplfold option -noLP
my $P;          # RNAplfold option -P (paramFile)
my $T;

# global vars
my $currDir = getcwd;
my $currUsr = getlogin;
my (undef, $path_to_skript) = fileparse($0);
$path_to_skript = abs_path($path_to_skript);
my $len;
my @pairprobs_wl;
my @PuW;

# standard settings
my $stdwindow = 200;
my $stdmaxspan = 150;
my $stdcutoff = 0;
my $stdskipbordernts = 10;
my $stdu = 1;
my $stdT = 37;

# check if RNAplfold is in PATH environment variable
`RNAplfold --help 2>&1 /dev/null`;
if ($? != 256 && $? != 0) {
	print STDERR "error: please ensure that RNAplfold is installed and in the current PATH\n";
	exit;
}

# create temporary directory that is deleted after exit
my $tmp_template = 'local-fold-XXXXXX';
my $tmp_prefix = '/var/tmp/';
my $tmpdir = tempdir($tmp_template, DIR => $tmp_prefix, CLEANUP => 1);
# clean up temporary files when interrupted
$SIG{'INT'} = 'end_handler';
$SIG{'TERM'} = 'end_handler';
sub end_handler {
	print "signal ", $_[0], " caught, cleaning up temporary files\n";
	File::Temp::cleanup();
	die();
}

###############################################################################
# set pairing probability from RNAplfold
# in: winpos, i, j, pairing probability
###############################################################################
sub setPijWL{
	my ($winpos,$i,$j,$pij)=@_;
	$pairprobs_wl[$winpos][$i][$j]=$pij;
}

###############################################################################
# retrieve pairing probability from RNAplfold
# in: winpos, i, j
# out: pairing probability
###############################################################################
sub getPijWL{
	my ($winpos,$i,$j)=@_;
	if ($pairprobs_wl[$winpos][$i][$j]) { 
		return $pairprobs_wl[$winpos][$i][$j];
	} else { 
		return 0; 
	} 
}

###############################################################################
# set unpaired probability from RNAplfold
# in: offset, position, unpaired probability
###############################################################################
sub setPuW{
	my ($offset, $pos, @Probs)=@_;
	$PuW[$offset][$pos]=\@Probs;
}

###############################################################################
# Parse a RNAplfold unpaired basepairs _lunp file
# in: filename, offset
###############################################################################
sub parse_lunp  {
    my ($filename, $offset) = @_;
    local *IN;
    
    open(IN,$filename) || die "Cannot read $filename for parsing as _lunp file.\n";
   
	<IN>; <IN>; # skip comments    
    while (my $line = <IN>) {
		chomp $line;
		my @probs = split "\t", $line;
		my $pos = (shift @probs) - 1; 	# convert to 0-base
        setPuW($offset, $pos, @probs);
    }
    
    close IN;
}

###############################################################################
# retrieve unpaired probability from RNAplfold
# in: offset, position, unpaired probability
# out: pairing probability
###############################################################################
sub getPuW{
	my ($offset, $pos)=@_;
	
	# load data if needed
	if (not defined $PuW[$offset]) {
		my $flunp = "$tmpdir/p\_$offset\_lunp";
		($debug) and print STDERR "loading PuW at offset $offset\n";
		parse_lunp $flunp, $offset;
	}
	
	if (defined $PuW[$offset][$pos]) {
		return @{$PuW[$offset][$pos]};
	} else {
	    die "error: invalid entry for getPuW: offset=$offset, pos=$pos";
	} 
}

###############################################################################
# clear unpaired probabilities at offset to save memory
# in: offset
###############################################################################
sub clearPuW{
	my ($offset) = @_;
	if (defined $PuW[$offset]) {
		($debug) and print STDERR "clearing PuW at offset $offset\n";
		$PuW[$offset]=();
	}
}

###############################################################################
# parse first id and sequence from fasta file
# in: filename
# out: reference to hash with id => sequence
###############################################################################
sub parse_fasta {
	my ($seqfile) = @_;
	
	my %seq;
	open FASTA, $seqfile or die "$!";

	my $id;
	my $seq;
	my $line;
	IDSTART: while ($line = <FASTA>) {
		chomp $line;
		next if (length($line) == 0); 		# skip empty lines
		if (substr($line, 0, 1) eq '>') {
			$id = substr($line,1);
			while ($line = <FASTA>) {
				chomp $line;
				next if (length($line) == 0);  		# skip empty lines
				# is this a comment?
				if (substr($line, 0, 1) eq '>') {
					($id) = split(" ", $id);
					if (not defined $seq{$id}) {
						$seq{$id} = $seq;
					} else {
						print STDERR "error: id $id is not unique\n";
					}
					$seq = "";
					redo IDSTART;
				} elsif ($line =~ /[^NAGCUTnagcut]/) {
					print STDERR "error, neither comment nor sequence:\n";
					print STDERR $line, "\n";
					die;
				} else {
					$seq .= $line;
				}
			}
		} else {
			print STDERR "error: expecting id\n";
			print STDERR "$line\n";
			die;
		}
	}
	($id) = ($id =~ /\s*(\S+)/);
	if (not defined $seq{$id}) {
		$seq{$id} = $seq;
	} else {
		print STDERR "error: id $id is not unique\n";
	}
	
	return \%seq;
}

###############################################################################
# Parse a "RNAplfold"-generated dotplot postscript file
# in: filename, offset
###############################################################################
sub parse_dp_ps_RNAplfold {
	my ($filename, $offset) = @_;
	local *IN;

	open(IN,$filename) || die "Cannot read $filename for parsing as dp-ps file.\n";
	my $seq="";

	while (my $line=<IN>) {
		if ($line =~ /^\/sequence \{ \(/) {
			while (defined($line = <IN>) && ($line !~  /\} def/  ))  {
				chomp $line;
				$line =~ s/\\$//;
				$seq .= $line;
			}
		}
		if ($line =~ /(\d+)\s+(\d+)\s+(\S+)\s+ubox/) {
			setPijWL($offset, $1-1, $2-1, $3*$3);	# save & convert to 0-base
		}
	}
	close IN;

	$seq ne "" || die "Empty sequence in dp.ps file $filename\n";
}

##############################################################################
# simulate RNAplfold using RNAplfold
# in: sequence, window
##############################################################################
sub callPlfold {
	my ($fasta_seq, $window, $maxlenbp, $u, $outfile_name) = @_;
	
	if($noLP){
		$noLP = "-noLP";
	} else {
		$noLP = "";
	}

	my $lastpos = length($fasta_seq)-$window;

	# create temporary fasta file containing all subsequences of size window
	open TMPFASTA, ">$tmpdir/tmppl.fa" or die "error: couldn't open $tmpdir/tmppl.fa for writing";
	for (my $offset=0; $offset<=$lastpos; $offset++) {
		print TMPFASTA ">p\_$offset\n";
		print TMPFASTA substr($fasta_seq, $offset, $window), "\n";
	}
	close TMPFASTA;

	# compute pair probabilities and ensemble energies
	print STDERR "calling RNAplfold:";
	chdir $tmpdir;
	if ($P){
		system("RNAplfold $noLP -P $P -c 0 -d2 -u $u -W $window -L $maxlenbp -T $T ".
			"< $tmpdir/tmppl.fa > $tmpdir/ensenergy_wl.out");
	} else {
		system("RNAplfold $noLP -c 0 -d2 -u $u -W $window -L $maxlenbp -T $T ".
			"< $tmpdir/tmppl.fa > $tmpdir/ensenergy_wl.out");
	}
	chdir $currDir;
	print STDERR " finished!\n";
	
	# read pairing and unpairedness probabilities
	($debug) and print "scanning probabilities: ";
	for (my $offset=0; $offset<=$lastpos; $offset++) {
		my $fname = "$tmpdir/p\_$offset\_dp.ps";
		($debug) and print $fname, ", ";
		($nodot) or parse_dp_ps_RNAplfold $fname, $offset;
	}
	($debug) and print "\n";
	print "\n";
}

###############################################################################
# add RNAplfold dotplot postscript around coordinates
# in: coordinate file, filename final ps file, $sequence_id, $sequence, $window
###############################################################################
sub wrapPostScript {
	my ($infile, $outfile, $fasta_id, $fasta_seq, $tagline) = @_;
	
	open F, '>', "$tmpdir/ps_stuff";
	print F '270 665 moveto /Helvetica findfont 14 scalefont setfont ($seqfile) show',"\n";
	print F "\n";
	print F '/sequence { (\\',"\n";
	print F "$fasta_seq\\","\n";
	print F ') } def',"\n";
	print F "/winSize ",$maxspan," def","\n";
	close F;
	
	my $exec = "cat $path_to_skript/data/pstemplate.one $tmpdir/ps_stuff $path_to_skript/data/pstemplate.two $infile $path_to_skript/data/pstemplate.three | sed 's/RNA Dot Plot/$tagline/' > $outfile";
	print STDERR $exec,"\n" if ($debug);
	system($exec);
}

###############################################################################
# compute averages using L parameter for maximal span of bps
# in: fasta id, fasta sequence
###############################################################################
sub compute_dots {
    my ($outfile_prefix, $fasta_seq, $plwindow, $maxlenbp, $fasta_id) = @_;
    
    # average probabilities
    open CENTERED, ">$tmpdir/centered_wl.coordinates";

	# set boundaries for centered approach
	my $first_admissible = $skipbordernts;
	my $last_admissible = $window - $skipbordernts -1;
    
    # compute for each basepair (i,j) that fits into window
    for (my $i = 0; $i < $len-3; $i++) {
        for (my $j = $i+2; $j < min($i+$maxlenbp, $len); $j++) {
        	# centered probability
        	my $PCij = 0;
        	# number of centered windows
        	my $WCij = 0;

		    # compute first and last window a la plfold
		    my $beg=max(0,$j-$plwindow+1);
            my $end=min($i,$len-$plwindow);

            # sum up probabilities and weights
            for (my $offset = $beg; $offset<=$end; $offset++) {
            	# centered
            	my $relstart = $i-$offset;
            	my $relend = $j-$offset;
            	if (
					# use only basepairs that do not begin or end at the border
					($relstart >= $first_admissible and $relend <= $last_admissible) or
					# handle sequence borders separately to avoid NAs
					($i < $skipbordernts) or ($j > $len - $skipbordernts -1)
					) {
						# centered probability
						$PCij += getPijWL($offset, $i-$offset, $j-$offset);
						# number of centered windows
						$WCij++;
				}
            }
            
            # compute final value
        	$PCij = $PCij/$WCij;
            
            # write dotplot            
            if ($PCij >= $cutoff) {
            	my $float = sprintf("%.4f", sqrt($PCij));
			    print CENTERED $i+1, " ", $j+1, " ", $float, " ubox", "\n";
           	}
        }
    }
    close CENTERED;

    wrapPostScript(
		"$tmpdir/centered_wl.coordinates",
		"${outfile_prefix}\_W$window\_L$maxspan\_skip$skipbordernts.ps",
		$fasta_id, $fasta_seq,
		"$fasta_id. LocalFold W$window L$maxspan skip$skipbordernts");
}

###############################################################################
# average unpaired probabilities using Boltzmann-weights
# in: fasta id, window size, filename suffix for secondary
#     structure element accessibilities
###############################################################################
sub compute_accs {
    my ($fasta_id, $plwindow) = @_;

    # open output files
	my $centered_acc_file = "${fasta_id}\_W$window\_L$maxspan\_skip$skipbordernts.acc";
	open CACC, ">$centered_acc_file" or die "error: couldn't open $centered_acc_file";
	
	# set boundaries for centered approach
	my $first_admissible = $skipbordernts;
	my $last_admissible = $window - $skipbordernts -1;

    # for all positions in sequence
    for (my $i=0; $i < $len; $i++) {
		
		# remove unneeded PuW data (setting of upper bound of offset _very_ generous)
		for (my $offset = 0; $offset < ($i - $plwindow); $offset++) {
			clearPuW($offset);
		}
		
		my @cprobs = ();	# centered probabilities

		# for all ranges of unpaired bps
		foreach my $nup (1..$u) {
			# compute first position that is incorporated into prob
			my $firstpos = $i-$nup+1;
			# skip if this position is invalid
			if ($firstpos < 0) {
				push @cprobs, "NA";
				next;
			}
			# compute positions of first and last admissible windows
			my $beg = max(0,$i-$plwindow+1);
			my $end = min($firstpos,$len-$plwindow);
			# compute probabilities
			my $cprob=0;	# centered probability
			my $cwins=0;	# number of centered windows
			for (my $offset = $beg; $offset<=$end; $offset++) {
				
				my @pus = getPuW($offset, $i-$offset);
				
				my $relstart = $i-$offset;
				if (
					# use only unpaired stretches that do not begin or end at the border
					(($relstart - $u + 1 >= $first_admissible) and
					($relstart <= $last_admissible)) or
					# handle sequence borders separately to avoid NAs
					($i -$u +1 < $skipbordernts) or ($i > $len - $skipbordernts -1)
					) {
						# add up centered probabilities
						$cprob += $pus[$nup-1];
						$cwins++;
				}
				
			}
			# add final probabilities
			push @cprobs, $cprob/$cwins;
		}
		# save to file
		print CACC join ("\t", $i+1, @cprobs),"\n";

	}

	# clean up
	close CACC;
}

###############################################################################
# parse command line options
###############################################################################

# save command line options for later use
my $options = join(" ", @ARGV);
print STDERR "complete line of options supplied was: '$options'\n" if ($debug);
my $result = GetOptions (	"help"			=> \$help,
							"man"			=> \$man,
							"seqfile=s" 	=> \$seqfile,
							"W|window=i" 	=> \$window,
							"L|maxspan=i" 	=> \$maxspan,
							"u=i"       	=> \$u,
							"skipbordernts=i" => \$skipbordernts,
							"cutoff=f" 		=> \$cutoff,
							"nodot"			=> \$nodot,
							"noacc"			=> \$noacc,   
							"debug"			=> \$debug,
							"noLP"			=> \$noLP,
							"P=s"             	=> \$P,
							"T=f"			=> \$T);
							
							
pod2usage(-exitstatus => 1, -verbose => 1) if $help;
pod2usage(-exitstatus => 0, -verbose => 2) if $man;
($result) or pod2usage(2);
($seqfile) or pod2usage("error: seqfile is mandatory");
(-f $seqfile) or pod2usage("error: no such file '$seqfile'");
($window) or ($window=$stdwindow);
($window > 1) or pod2usage("error: window size too low");
($maxspan) or ($maxspan=$stdmaxspan);
($skipbordernts) or ($skipbordernts = $stdskipbordernts);
($skipbordernts >= 0) or pod2usage("error: skipbordernts has to be >= 0");
($u) or ($u = $stdu);
($cutoff) or ($cutoff = $stdcutoff);
($cutoff >= 0 and $cutoff < 1) or pod2usage("cutoff has to be in [0,1)");
($T) or $T = $stdT;

# parse input filename
my ( $fasta_prefix, $path, $suffix ) = fileparse( $seqfile, "\.[^.]*" );

# parse input fasta
my $hashref = parse_fasta($seqfile);
my %seqs = %{$hashref};
# if fasta contains one sequence use its filename for the output files
# instead of the sequence id
my $use_fastafname_for_output = (scalar keys %seqs == 1);

# main loop
while (my ($fasta_id,$fasta_seq) = each %seqs) {
	my $outfile_name = $use_fastafname_for_output ? $fasta_prefix : $fasta_id;
	$len = length($fasta_seq);

	# set window and span for short sequences
	my $seqwin = $window;
	my $seqspan = $maxspan;
	if ($seqwin > $len) {
		print STDERR "warning: sequence shorter than window, decreasing window size\n";
		$seqwin = $len;
		if ($seqspan > $len){
			print STDERR "warning: span too large, decreasing span\n";
			$seqspan = $len - 2 * $skipbordernts;
		}
	}
	
	print STDERR "sequence id: $fasta_id\n";
	print STDERR "sequence length: $len\n";
	print STDERR "window: $seqwin\n";
	print STDERR "maxspan: $seqspan\n";
	print STDERR "u: $u\n";
	print STDERR "skipbordernts: $skipbordernts\n";
	print STDERR "temperature: $T\n";
	
	# compute probabilities
	mkdir $tmpdir;
	callPlfold($fasta_seq, $seqwin, $seqspan, $u, $outfile_name);

	# compute dotplots
	if (not $nodot) {
		compute_dots($outfile_name, $fasta_seq, $seqwin, $seqspan, $fasta_id);
	}

	# compute accessibilities
	if (not $noacc) {
		compute_accs($outfile_name, $seqwin);
	}
		
	# new sequence: clean out everything
	rmdir $tmpdir;
	@PuW = ();
	@pairprobs_wl = ();
}
