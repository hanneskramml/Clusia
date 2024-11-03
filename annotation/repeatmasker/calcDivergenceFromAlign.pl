#!/lisc/app/perl/5.38.2/bin/perl
##---------------------------------------------------------------------------##
##  File:
##      @(#) calcDivergenceFromAlign.pl
##  Author:
##      Robert M. Hubley   rhubley@systemsbiology.org
##  Description:
##      A utility script to calculate a new divergence measure on the
##      RM alignment files.
##
#******************************************************************************
#* Copyright (C) Institute for Systems Biology 2003-2009 Developed by
#* Arian Smit and Robert Hubley.
#*
#* This work is licensed under the Open Source License v2.1.  To view a copy
#* of this license, visit http://www.opensource.org/licenses/osl-2.1.php or
#* see the license.txt file contained in this distribution.
#*
#******************************************************************************
#
# ChangeLog
#
#     $Log$
#
###############################################################################
#
# To Do:
#

=head1 NAME

calcDivergenceFromAlign.pl - Recalculate/Summarize the divergences in an align file.

=head1 SYNOPSIS

  calcDivergenceFromAlign.pl [-version] [-s <summary_file>] [-noCpGMod]
                             [-a <new_align_file>] 
                             *.align[.gz]

  Typical use case: generate *.divsum file for createRepeatLandscape.pl
      
       ./calcDivergenceFromAlign.pl -s mygenome.divsum  mygenome.fa.align.gz

  Note: The "-a" parameter is only necessary if you want to save the
        per-alignment divergence values in addition to the family summaries.

=head1 DESCRIPTION

  A utility script to calculate a new divergence measure on the
  RM alignment files.  Currently RepeatMasker only calculates the
  standard Kimura 2-Parameter divergence metric for the *.align file. 

  The new divergence metric is modification of the Kimura 2-Parameter
  model where "CG" dinucleotide sites in the consensus sequence are
  treated specially:

    - Two transition mutations are counted as a single 
      transition,
    - One transition is counted as 1/10 of a standard
      transition, and 
    - Transversions are counted normally (as they would 
      outside of a CpG site).  

  This modification to the Kimura 2 parameter model accounts for the 
  extremely high rate of mutations in at a CpG locus. Note this is 
  applicable to organisms for which CpG methylation is applicable.


The options are:

=over 4

=item -version

Displays the version of the program

=item -s <summary_file>

Generate a file containing the summarized divergences for each family.  This 
is needed by the createRepeatLandscape.pl tool.

=item -a <new_align_file>

Optionally generate a new *.align file with each alignment labeled with 
the recalculated divergence.

=item -noCpGMod

Do not modify the transition counts at CpG sites. 

=back

=head1 SEE ALSO

=head1 COPYRIGHT

Copyright 2013-2023 Robert Hubley, Institute for Systems Biology

=head1 AUTHOR

Robert Hubley <rhubley@systemsbiology.org>

=cut

#
# Module Dependence
#
use strict;
use Getopt::Long;
use Data::Dumper;
use FileHandle;

## TODO: Remove this
use lib "/home/rhubley/projects/RepeatMasker";
use FindBin;
use lib $FindBin::RealBin;
use lib "$FindBin::Bin/..";
use RepeatMaskerConfig;
use SearchResult;
use CrossmatchSearchEngine;

#
# Version
#
my $Version = $RepeatMaskerConfig::VERSION;

#
# Magic numbers/constants here
#  ie. my $PI = 3.14159;
#
my $DEBUG = 0;

#
# Option processing
#  e.g.
#   -t: Single letter binary option
#   -t=s: String parameters
#   -t=i: Number paramters
#
my @getopt_args = (
                    '-version',    # print out the version and exit
                    '-noCpGMod',
                    '-a=s',
                    '-s=s'
);

my %options = ();
Getopt::Long::config( "noignorecase", "bundling_override" );
unless ( GetOptions( \%options, @getopt_args ) ) {
  usage();
}

sub usage {
  print "$0 - $Version\n\n";
  exec "pod2text $0";
  exit;
}

if ( $options{'version'} ) {
  print "$Version\n";
  exit;
}

if ( !$options{'s'} && !$options{'a'} ) {
  print
"\n===\n=== Error: One or more of the options '-a' or '-s' must be supplied!\n===\n";
  usage();
}

if ( ! -s $ARGV[ 0 ] ) {
  print "\n===\n=== Error: Missing alignment file parameter!\n===\n";
  usage();
}

my $alignFile      = $ARGV[ 0 ];
my $maxDiv         = 70;
my $cntAlign       = 0;
my %repeatMuts     = ();
my %classDivWCLen  = ();
my $prevQueryName  = "";
my $prevQueryBegin = "";
my $prevQueryEnd   = "";
my $prevHitName    = "";
my $prevDiv        = "";
my $prevClass      = "";

my $searchResultsFH = new FileHandle;

if ( $alignFile =~ /.+\.gz/ ) {
  open $searchResultsFH, "gunzip -c $alignFile|"
      or die
      "calcDivergenceFromAlign: Could not open gunzip for reading $alignFile: $!\n";
}
else {
  open $searchResultsFH, "<$alignFile"
      or die "calcDivergenceFromAlign: Could not open $alignFile for reading: $!\n";
}

if ( $options{'s'} ) {
  open SOUT, ">$options{'s'}"
      or die "Error: Could not open $options{'s'} for writing!\n";
#  if ( $options{'noCpGMod'} ) {
#    print SOUT "Jukes/Cantor and Kimura subsitution levels\n";
#    print SOUT "==========================================\n";
#  }
#  else {
#    print SOUT
#        "Jukes/Cantor and Kimura subsitution levels adjusted for CpG sites\n";
#    print SOUT
#        "=================================================================\n";
#  }
#  print SOUT "File: " . $alignFile . "\n";
}

#
# Process the alignment file
#
my $outAlign = 0;
if ( $options{'a'} ) {
  open COUT, ">$options{'a'}"
      or die "Could not open $options{'a'} for writing!\n";
  $outAlign = 1;
}

CrossmatchSearchEngine::parseOutput( searchOutput => $searchResultsFH,
                                     callback     => \&processAlignment );

if ( $options{'a'} ) {
  close COUT;
}

if ( $options{'s'} ) {

#  print SOUT "Weighted average Kimura divergence for each repeat family\n";
#  print SOUT "chr\trepeat\tabslen\twellcharlen\tkimura\n";
#  print SOUT "-----\t------\t------\t-----------\t-------\n";

   print SOUT "chr\trepeat\tbin\tbp\n";

   foreach my $seqName ( sort keys %repeatMuts ) {
     foreach my $id ( sort keys %{ $repeatMuts{$seqName} } ) {
#       my $kimura = 100;
#       if ( $repeatMuts{$seqName}->{$id}->{'wellCharLen'} > 0 ) {
#         $kimura = sprintf( "%4.2f",
#                            $repeatMuts{$seqName}->{$id}->{'sumdiv'} /
#                                $repeatMuts{$seqName}->{$id}->{'wellCharLen'} );
#
#         $kimura = $maxDiv if ( $kimura > $maxDiv );
#       }
#
#       if ( $seqName =~ /Simple|Low_complexity|ARTEFACT/ ) {
#         print SOUT "$seqName\t$id\t"
#             . $repeatMuts{$seqName}->{$id}->{'absLen'} . "\t"
#             . $repeatMuts{$seqName}->{$id}->{'wellCharLen'}
#             . "\t----\n";
#       }
#       else {
#         print SOUT "$seqName\t$id\t"
#             . $repeatMuts{$seqName}->{$id}->{'absLen'} . "\t"
#             . $repeatMuts{$seqName}->{$id}->{'wellCharLen'}
#             . "\t$kimura\n";
#       }

        my $j = 0;
        while ( $j <= $maxDiv ) {
          my $label = "$id $j";
          $classDivWCLen{$label} = 0 unless $classDivWCLen{$label};
          print SOUT "$seqName\t$id\t$j\t$classDivWCLen{$label}\n";
          ++$j;
        }

     }
   }

#  print SOUT "\n\n";
#
#  print SOUT "Coverage for each repeat class and divergence (Kimura)\n";
#  print SOUT "Div ";
#  foreach my $class ( sort keys %repeatMuts ) {
#    print SOUT "$class ";
#  }
#  print SOUT "\n";
#
#  my $j = 0;
#  while ( $j <= $maxDiv ) {
#    print SOUT "$j ";
#    foreach my $seqName ( sort keys %repeatMuts } ) {
#      my $label = "$seqName $j";
#      $classDivWCLen{$label} = 0 unless $classDivWCLen{$label};
#      print SOUT "$classDivWCLen{$label} ";
#    }
#    print SOUT "\n";
#    ++$j;
#  }
  close SOUT;
}

exit;

######################## S U B R O U T I N E S ############################

##-------------------------------------------------------------------------##
## Use: my processAlignment( $parameter => value );
##
##      $parameter       : A parameter to the method
##
##  Returns
##
##-------------------------------------------------------------------------##
sub processAlignment {
  my $result = shift;

  return if ( !$result );

  my $hitname;
  my $class;
  my $subjName = $result->getSubjName();
  if ( $subjName =~ /(\S+)\#(\S+)/ ) {
    $hitname = $1;
    $class   = $2;
  }
  else {
    $hitname = $subjName;
    $class   = $result->getSubjType();
  }

  # JR 20200116: Some combinations of RepeatMasker/Dfam can produce
  # an empty class name, particularly in the UCON4 family.
  if ($class eq "") {
    $class = "Unspecified";
  }

  my $seqName    = $result->getQueryName();
  my $queryStart = $result->getQueryStart();
  my $queryEnd   = $result->getQueryEnd();

  # Simple repeats, low complexity and artefacts should not be counted
  #if ( $class =~ /Simple|Low_complexity|ARTEFACT/ )
  #{
  #  if ( $outAlign )
  #  {
  #    print COUT ""
  #      . $result->toStringFormatted( SearchResult::AlignWithQuerySeq ) . "\n";
  #  }
  #  return;
  #}

  print STDERR "." if ( $cntAlign++ % 1000 == 0 );

  #return if ( !( $seqName =~ /^(CMU01|CMU02|CMU03|CMU04|CMU06|CMU07|CMU10|CMU11|CMU12|CMU14|CMU16|CMU19|CMU20|CMU25)$/ ) );  # H1
  #return if ( !( $seqName =~ /^(CMU09|CMU05|CMU22|CMU15|CMU08|CMU28|CMU24|CMU21|CMU17|CMU13|CMU18|CMU27|CMU29|CMU23|CMU26)$/ ) );  # H2
  #return if ( !( $seqName eq "CMU24" ) );

  my ( $div, $transi, $transv, $wellCharBases, $numCpGs );

  my $alen = $queryEnd - $queryStart + 1;
  $wellCharBases = $alen - int( $alen * ( $result->getPctInsert() / 100 ) );

  if ( $class =~ /Simple|Low_complexity|ARTEFACT/ ) {
    $div     = 100;
    $hitname = "combined";
  }
  else {

    # Obtain divergence from modern *.align files directly
    $div = $result->getPctKimuraDiverge();

    if ( $div eq "" ) {

      # Calculate divergence on the fly
      ( $div, $transi, $transv, $wellCharBases, $numCpGs ) =
          $result->calcKimuraDivergence( divCpGMod => 1 );
      $result->setPctKimuraDiverge( sprintf( "%4.2f", $div ) );
    }
  }

  if ( $prevQueryName eq $seqName ) {
    if ( $prevQueryEnd > $queryStart ) {

      # Overlap
      my $overlapAbsLen = $prevQueryEnd - $queryStart + 1;
      if ( $prevQueryEnd >= $queryEnd ) {
        if ( $outAlign ) {
          print COUT ""
              . $result->toStringFormatted( SearchResult::AlignWithQuerySeq )
              . "\n";
        }
        return;
      }
      if ( $div > $prevDiv ) {

        # Previous gets overlap bases - subtract overlap from this hit
        $wellCharBases -= $overlapAbsLen;
        $alen          -= $overlapAbsLen;
      }
      else {

        # Current gets overlap bases - subtract overlap from previous

        # TODO: (JR 20191003) This is not quite correct, because it assumes
        # the previous hit got the overlap. That isn't necessarily always
        # the case, for example where there are 3 or more overlapping hits.
        # This usually introduces only relatively minor errors, however.

        my $key = "$prevHitName $prevDiv";
        $classDivWCLen{$key}                        -= $overlapAbsLen;
        $repeatMuts{$prevQueryName}->{$prevHitName}->{'sumdiv'} -=
            $prevDiv * $overlapAbsLen;
        $repeatMuts{$prevQueryName}->{$prevHitName}->{'wellCharLen'} -= $overlapAbsLen;
        $repeatMuts{$prevQueryName}->{$prevHitName}->{'absLen'}      -= $overlapAbsLen;
      }
    }
  }
  $prevQueryName  = $seqName;
  $prevDiv        = $div;
  $prevQueryBegin = $queryStart;
  $prevQueryEnd   = $queryEnd;
  $prevHitName    = $hitname;
  $prevClass      = $class;

  $repeatMuts{$seqName}->{$hitname}->{'sumdiv'}      += $div * $wellCharBases;
  $repeatMuts{$seqName}->{$hitname}->{'wellCharLen'} += $wellCharBases;
  $repeatMuts{$seqName}->{$hitname}->{'absLen'}      += $alen;
  $div = int( $div );
  my $key = "$hitname $div";
  $classDivWCLen{$key} += $wellCharBases;

  if ( $outAlign ) {
    print COUT ""
        . $result->toStringFormatted( SearchResult::AlignWithQuerySeq ) . "\n";
  }

}

1;
