package Statistics::SDT;

use strict;
use warnings;
use Carp qw(croak);
use vars qw($VERSION @EXPORT_OK);
use Exporter qw( import );
$VERSION = 0.02;
@EXPORT_OK = qw(d_sensitivity a_sensitivity area_d_sensitivity likelihood_bias log_likelihood_bias decision_bias griers_bias dprime aprime area_dprime);

use Math::Cephes qw(:dists);

sub new {
	my ($class, $args) = @_;
	my $self = {};
	bless $self, $class;
    
    if (scalar keys %{$args}) {
        
        foreach (keys %{$args}) {
            $self->{$_} = $args->{$_};
        }
        
        ($self->{'hr'}, $self->{'far'}) = _init_rates($args);
    }

	return $self;
}

# --------------------
# Sensitivity measures:
# --------------------

sub d_sensitivity {
    my ($h, $f, $m) = _init_rates(@_);
    # Use Smith's algorithms if there are more than 2 alternatives:
    if ($m && $m > 2) {
        return _dprime_fc($h, $m);
    }
    else {
        # Assume d' = 0 if both rates = 0 or both = 1:
        if ( (!$h && !$f) || ($h == 1 && $f == 1) ) {
            return 0;
        }
        else {
            return ndtri($h) - ndtri($f);
        }
    }
}

# Smith's (1982) algorithms for forced-choice alternatives (m) greater than 2:
sub _dprime_fc {
    my ($h, $m) = @_;
   
    if ($m < 12) {
        my $km = .86 - .085 * log($m - 1);
        my $lm = ( ($m - 1) * $h) / (1 - $h);
        return $km * log($lm);
    }
    else {
        my $A = ( -4 + sqrt(16 + 25 * log($m - 1) ) )/3;
        my $B = sqrt( ( log($m - 1) + 2) / (log($m - 1) + 1) );
        return $A + ($B * ndtri($h));
    }
}

sub a_sensitivity {
    my ($h, $f) = _init_rates(@_);
    
    if ($h >= $f) {
        return (.5 + ( ($h - $f) * (1 + $h - $f) ) / ( 4 * $h * (1 - $f) ) );
    }
    else {
        return (.5 + ( ($f - $h) * (1 + $f - $h) ) / ( 4 * $f * (1 - $h) ) );
    }
    
}

sub area_d_sensitivity {
    my ($h, $f) = _init_rates(@_);
    
    # Assume A(d') = 0.5 if both rates = 0 or both = 1:
    if ( (!$h && !$f) || ($h == 1 && $f == 1) ) {
        return 0.5;
    }
    else {
        return ndtri(d_sensitivity($h, $f, 2) / sqrt(2));
    }
}

# --------------------
# Bias measures:
# --------------------

sub likelihood_bias { # beta
    my ($h, $f) = _init_rates(@_);
    return exp( ( ( (ndtri($f)**2) - (ndtri($h)**2) ) / 2 ) );
}

sub log_likelihood_bias { # ln(beta)
    my ($h, $f) = _init_rates(@_);
    return ( ( (ndtri($f)**2) - (ndtri($h)**2) ) / 2 );
}

sub decision_bias { # c
    my ($h, $f) = _init_rates(@_);
    return -1 *( ( ndtri($h) + ndtri($f) ) / 2 ) ;
}

sub griers_bias { # B''
    my ($h, $f) = _init_rates(@_);
    if ($h >= $f) {
        my $a = ( $h * (1 - $h) );
        my $b = ( $f * (1 - $f) );
#        if ($a + $b) {
            return ( $a - $b ) /  ( $a + $b );
#        }
#        else {
#            return 0;
#        }
    }
    else {
        my $a = ( $f * (1 - $f) );
        my $b = ( $h * (1 - $h) );
 #       if ($a + $b) {
            return ( $a - $b ) / ( $a + $b );
#        }
#        else {
#            return 0;
#        }
    }
}

# --------------------
# Initialise hit and false-alarm rates:
sub _init_rates {
# --------------------
    my $args = shift;
    foreach ([qw/hr hits signal_trials/], [qw/far false_alarms noise_trials/]) {
        if (! exists $args->{$_->[0]}) {
            
            # Need (i) no. of hits and signal trials, and (ii) no. of false alarms and noise trials:
            croak __PACKAGE__, "Number of $_->[1] and $_->[2] needed to calculate $_->[0]"
            if ! exists $args->{$_->[1]} || ! exists $args->{$_->[2]};
            
            # Apply the "loglinear" correction, regardless of values:
            if ($args->{'correct'} and $args->{'correct'} > 1) {
                $args->{$_->[0]} = ($args->{$_->[1]} + .5) / ($args->{$_->[2]} + 1);
            }
            # or get the rate first, applying corrections if needed:
            else {
                $args->{$_->[0]} = $args->{$_->[1]} / $args->{$_->[2]};
                
                if ($args->{'correct'}) {
                    if (! $args->{$_->[0]}) {
                        $args->{$_->[0]} = .5 / $args->{$_->[2]};
                    }
                    elsif ($args->{$_->[0]} == 1) {
                        $args->{$_->[0]} = ($args->{$_->[2]} - .5) / $args->{$_->[2]};
                    }
                }
            }    
        }
    } # end foreach loop
    return ($args->{'hr'}, $args->{'far'}, $args->{'alternatives'});
}

# --------------------
# Function aliases
# --------------------
*dprime = \&d_sensitivity;
*aprime = \&a_sensitivity;
*area_dprime = \&area_d_sensitivity;

1;

__END__

=head1 NAME

Statistics::SDT - Signal detection theory measures of sensitivity and response-bias

=head1 VERSION

This is documentation for Version 0.02 of Statistics-SDT (2006-11-20).

=head1 SYNOPSIS

    use Statistics::SDT;

    $sdt = Statistics::SDT->new(
        { 
            hits => 50,
            signal_trials => 50,
            false_alarms => 17,
            noise_trials => 25,
            correct => 2,
        }
    );

    $d = $sdt->d_sensitivity();
    $c = $sdt->decision_bias();

=head1 DESCRIPTION

Signal Detection Theory algorithms (e.g., of d', A', decision bias), as prescribed by Stanislav & Todorov (1999). Both object- and function-oriented interfaces are provided.

=head1 KEY VALUES

For both object- and function-oriented styles, the following named parameters must be given as a hash-reference: either to the L<new|new> constructor method, or (with the function-oriented style) into each function. Basically, either all of the first four parameters are required (in order to calculate the hit-rate and false-alarm-rate), or the required rates are themselves supplied.

=over 4

=item hits

The number of hits.

=item false_alarms

The number of false alarms.

=item signal_trials

The number of signal trials. The hit-rate is derived by dividing the number of hits by the number of signal trials.

=item noise_trials

The number of noise trials. The false-alarm-rate is derived by dividing the number of false-alarms by the number of noise trials.

=item alternatives

The number of response alternatives. Default = 2 (for the classic signal-detection situation of discriminating between signal+noise and noise-only). If the number of alternatives is greater than 2, the measure of sensitivity, when calling L<d_sensitivity|d_sensitivity>, is based on the Smith (1982) algorithms.

=item correct

A parameter that indicates whether or not to perform a correction on the number of hits and false-alarms as a corrective when the hit-rate or false-alarm-rate equals 0 or 1 (due, e.g., to strong inducements against false-alarms, or easy discrimination between signals and noise). This is relevant to all functions that make use of the I<inverse phi> function (all except L<a_sensitivity|a_sensitivity> and L<griers_bias|griers_bias>). 

If set to greater than 1, the loglinear transformation is applied, i.e., 0.5 is added to both the number of hits and false-alarms, and 1 is added to the number of signal and noise trials. These adjustments are made irrespective of the extremity of the rates themselves.

If set to 1, extreme rates (of 0 and 1, only) are replaced with the number of signal/noise trials, moderated by a value of 0.5 (specifically, where I<n> = number of signal or noise trials: 0 is replaced with 0.5 / I<n>; 1 is replaced with (I<n> - 0.5) / I<n>.

Stanislav and Todorov (1999) advise that the latter correction is the most common method of handling extreme rates, but that it might bias sensitivity measures and not be as satisfactory as the loglinear transformation applied to all hits and false-alarms.

If set to zero (the default), no correction is performed to the calculation of the rates. This should only be used when you are using (1) the parametric measures and are sure the rates are not at the extremes of 0 and 1; or (2) the nonparametric algorithms (L<a_sensitivity|a_sensitivity> and L<griers_bias|griers_bias>). An alternative to these corrections is, indeed, to use the nonparametric measures.

=item hr

This is the hit-rate. Instead of passing the number of hits and signal trials, give the hit-rate directly - but, if doing so, ensure the rate does not equal zero or 1 in order to avoid errors thrown by the inverse-phi function (which will be given as "ndtri domain error").

=item far

This is the false-alarm-rate. Instead of passing the number of false alarms and noise trials, give the false-alarm-rate directly - but, if doing so, ensure the rate does not equal zero or 1 in order to avoid errors thrown by the inverse-phi function (which will be given as "ndtri domain error").

=back

=head1 METHODS

The methods can be accessed by an object- or function-oriented style: 

B<Object-oriented style>: firstly construct a class object by supplying L<KEY VALUES|KEY VALUES> to the L<new|new> method, and then call each method by supplying the class object as the first argument.

 $sdt = Statistics::SDT->new({hr => => 10/12, far => 1/12});
 $bias = $sdt->likelihood_bias();

B<Function-oriented style>: firstly explicitly import the functions to be C<use>d, and then simply call the function, with L<KEY VALUES|KEY VALUES>.

 use Statistics::SDT qw(likelihood_bias);
 $bias = likelihood_bias({hr => => 10/12, far => 1/12});

=head2 Class constructor

=head3 new

This method may be used to construct a class object that holds the values of the parameters, as above, and to call the following methods, without having to resubmit all the values. 

As well as holding the values of the parameters submitted to it, the class-object returned by C<new> will hold two arguments, B<hr>, the hit-rate, and B<far>, the false-alarm-rate.

This method can be skipped, and the following methods used directly (all optionally exported), as long as the parameters, as above, are passed to each one.

You can supply the hit-rate and false-alarm-rate themselves, but ensure that they do not equal zero or 1 in order to avoid errors thrown by the inverse-phi function. The calculation of the hit-rate and false-alarm-rate by the module corrects for this limitation - see the notes on the C<correct> parameter, above. 

=head2 Sensitivity measures

=head3 d_sensitivity

Returns the index of sensitivity, or discrimination, I<d'> (d prime). (C<dprime> is an alias for this function.)

=over 12

=item

I<d'> = phi^-1(hr) - phi^-1(far)

=back

I<d'> is found by subtracting the I<z>-score that corresponds to the false-alarm rate (B<far>) from the I<z>-score that corresponds to the hit rate (B<hr>). Sensitivity is measured in standard deviation units, larger positive values indicating greater sensitivity.

If both the hit-rate and false-alarm-rate are either 0 or 1, then C<d_sensitivity> returns 0.

If there are more than two alternatives (as specified by the parameter named I<alternatives> in the hash-reference passed to the L<new|new> constructor or this method, then Smith's (1982) algorithms are used.

A value of 0 indicates no sensitivity to the presence of the signal, i.e., it cannot be discriminated from noise. 

Values less than 0 indicate a lack of sensitivity that might result from response-confusion.

=head3 a_sensitivity

Returns the nonparametric index of sensitivity, I<A'>. (C<aprime> is an alias for this function.)

Ranges from 0 to 1. Values greater than 0.5 indicate positive discrimination (1 = perfect performance); values less than 0.5 indicate a failure of discrimination (perhaps due to response-confusion); and a value of 0.5 indicates no sensitivity to the presence of the signal, i.e., it cannot be discriminated from noise.

=head3 area_d_sensitivity

The area under the receiver-operating-characteristic (ROC) curve, simply equalling the proportion of correct responses that would have been yielded had the task been a two-alternative forced-choice task rather than a yes/no task. (C<area_dprime> is an alias for this function.)

If both the hit-rate and false-alarm-rate are either 0 or 1, then C<area_d_sensitivity> returns 0.5.

=head2 Bias measures

=head3 likelihood_bias

Returns the I<beta> measure of response bias, based on the ratio of the likelihood the decision variable obtains a certain value on signal trials, to the likelihood that it obtains the value on noise trials.

Values less than 1 indicate a bias toward the I<yes> response, values greater than 1 indicate a bias toward the I<no> response, and the value of 1 indicates no bias toward I<yes> or I<no>.

=head3 log_likelihood_bias

Returns the natural logarithm of the likelihood bias, I<beta>. 

Ranges from -1 to +1, with values less than 0 indicating a bias toward the I<yes> response, values greater than 0 indicating a bias toward the I<no> response, and a value of 0 indicating no response bias.

=head3 decision_bias

Implements the I<c> parametric measure of response bias. Ranges from -1 to +1, with deviations from zero, measured in standard deviation units, indicating the position of the decision criterion with respect to the neutral point where the signal and noise distributions cross over, there is no response bias, and I<c> = 0. 

Values less than 0 indicate a bias toward the I<yes> response (the decision variable exceeds the criterion); values greater than 0 indicate a bias toward the I<no> response (the decision variable is less than the criterion); and a value of 0 indicates no response bias.

=head3 griers_bias

Implements Griers I<B''> nonparametric measure of response bias. 

Ranges from -1 to +1, with values less than 0 indicating a bias toward the I<yes> response, values greater than 0 indicating a bias toward the I<no> response, and a value of 0 indicating no response bias.

=head1 EXAMPLES

=head2 Object-oriented style    

 use strict;

 use Statistics::SDT;

 # 1. Create SDT object, loading it with the key values:
 my $sdt = Statistics::SDT->new(
  { 
     hits           => 50,
     signal_trials  => 50,
     false_alarms   => 17,
     noise_trials   => 25,
     correct        => 2,
   }
 );

 # you could also give the hit and false-alarm rates directly:
 # my $sdt = Statistics::SDT->new({hr => 50/50, far => 17/25});
 # but this will throw an error when calling a parametric measure, 
 # as 1 is being given as the hit-rate.

 # 2. Use the SDT object to access methods (and the rates), e.g.: 
 print 'd = ' . $sdt->d_sensitivity() . "\n"; # d = 1.85864907492633
 print 'c = ' . $sdt->decision_bias() . "\n"; # c = -1.39702333657767
 print 'hr = ' . $sdt->{'hr'} . "\n"; # hr = 0.99

=head2 Function-oriented style

 # 1. Explicitly import the required methods:
 use Statistics::SDT qw(d_sensitivity decision_bias);

 # 2. Access the methods, supplying the key values to each, e.g.:
 my $d = d_sensitivity(
   { 
       hits => 10,
       signal_trials => 12,
       false_alarms => 0,
       noise_trials => 12,
       correct => 1,
    }
 );

 # you could also give the hit and false-alarm rates directly:
 # my $d = d_sensitivity({hr => 10/12, far => 0/12});
 # but this will throw an error, as zero is being given as the false-alarm-rate,
 # and a parametric measure of sensitivity is being called.

 my $bias = decision_bias(
    { 
       hits => 10,
       signal_trials => 12,
       false_alarms => 0,
       noise_trials => 12,
       correct => 1,
    }
 );

 print "d = $d\nbias = $bias\n";
 # d = 2.69908596222395
 # bias = 0.382121415010272

=head1 REFERENCES

Smith, J. E. K. (1982). Simple algorithms for M-alternative forced-choice calculations. I<Perception and Psychophysics>, I<31>, 95-96.

Stanislaw, H., & Todorov, N. (1999). Calculation of signal detection theory measures. I<Behavior Research Methods, Instruments, and Computers>, I<31>, 137-149.

=head1 SEE ALSO

L<Math::Cephes|lib::Math::Cephes> : The present module imports the L<ndtri|lib::Math::Cephes> or I<inverse phi> function from this package, which is used to calculate I<z>-scores from probabilities.

L<Statistics::ROC|lib::Statistics::ROC> : Receiver operating characteristic curves.

=head1 AUTHOR/LICENSE

This program is free software; you can redistribute it and/or modify it under the same terms as Perl itself.

Copyright (C) 2006 Roderick Garton  

=head1 DISCLAIMER

To the maximum extent permitted by applicable law, the author of this module disclaims all warranties, either express or implied, including but not limited to implied warranties of merchantability and fitness for a particular purpose, with regard to the software and the accompanying documentation.

=cut