package JCMT::DAS;

=head1 NAME

JCMT::DAS - Utility routines for the JCMT DAS correlator

=head1 SYNOPSIS

  use JCMT::DAS qw/ das_merge /;

  ($out, $freq, $vel) = das_merge( \@spectrum, \@f_cen, \@f_inc,
				   trunc => 0,
				   merge => 1,
				   clip => 0,
				 );

=head1 DESCRIPTION

Routines useful for manipulating DAS data.

=cut

use 5.006;
use strict;
use warnings;
use Carp;

use vars qw/ $VERSION @EXPORT_OK /;
$VERSION = '0.1';

use base qw(Exporter DynaLoader);

@EXPORT_OK = qw( das_merge );

bootstrap JCMT::DAS;

=head1 FUNCTIONS

=over 4

=item B<das_merge>

Despike and merge (no DC adjustment) DAS spectra.


  ($out, $freq, $vel) = das_merge( \@spectrum, \@f_cen, \@f_inc,
				   trunc => 0,
				   merge => 1,
				   clip => 0,
				 );

Arguments are:

 @spectrum : Reference to array containing the spectrum data
 @f_cen    : Reference to array contining the centre frequencies
             for each band in GHz
 @f_inc    : Frequency step between each subband in MHz

Optional hash arguments are:

 trunc     : number of channels/overlap to truncate
             Defaults to 50%
             0 < n < 1 keep n fraction of overlap
             n >= 1: drop n channels from ends

 merge     : If true, merge subbands. Defaults to no merge

 clip      : Amplitude cutoff. Defaults to 1000

Output arrays are:

  The merged spectrum
  The frequency axis
  The velocity axis

For normal DAS GSD data the spectrum can be read from component
C13DAT, the central frequencies can be read from C12CF and the
frequency increment can be read from C12FR.

=cut

sub das_merge {
  my $sp_in = shift;
  my $f_cen = shift;
  my $f_inc = shift;

  # Defaults
  my %def = ( trunc => 0, merge => 0, clip => 0 );

  # Read optional args
  my %opt  = ( %def, @_ );

  # Call the internal dasmerge routine
  my (@freq, @sp_out, @vel_out, $numvals);

  dasmerge( scalar(@$sp_in), scalar(@$f_cen), $f_cen, $f_inc,
	    $sp_in, $opt{trunc}, $opt{merge}, $opt{clip},
	    \@freq, \@sp_out, \@vel_out, $numvals
	  );

  # sanity check
  croak "Error in dasmerge: output spectrum has different size (".
    scalar(@sp_out).") to the reported size from the C routine ($numvals)"
      if  $numvals != scalar(@sp_out);

  # And return the results
  return (\@sp_out, \@freq, \@vel_out);

}

=back


=head1 AUTHORS

Tim Jenness E<lt>t.jenness@jach.hawaii.eduE<gt> for the perl
interface.

Remo Tilanus for the DAS merge C routine.

Array packing code supplied by Karl Glazebrook.

=cut


1;
