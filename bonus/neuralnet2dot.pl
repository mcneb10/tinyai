#!/usr/bin/env perl

use strict;
use warnings;
use Carp;
use Getopt::Long;
use File::Slurper 'read_text';
use File::Basename;
use Pod::Usage;

use lib dirname(__FILE__).'/lib';

use TinyAI::NeuralNet;

my $prune_useless = 0;

GetOptions('p|prune-useless!', \$prune_useless) || pod2usage(-exitval => 2, -verbose => 1);

# Must provide input and output file name
pod2usage(-exitval => 2, -verbose => 1) unless(@ARGV == 2);

my ($neural_net_filename, $graphviz_file_name) = @ARGV;

my $neural_net = NeuralNet->new(data => read_text($neural_net_filename));

open my $ofh, '>', $graphviz_file_name;
print $ofh $neural_net->dump_to_graphviz($prune_useless);
close $ofh;

__END__

=head1 NAME

neuralnet2dot.pl - Convert neural nets from tinyai to graphviz

=head1 SYNOPSIS

neuralnet2dot.pl [options] [input neural network] [output graphviz file]

=head1 OPTIONS

=over 4

=item B<--prune-useless, -p>

Hide nodes that don't affect the outputs.
The program never removes Input, Output, or Bias nodes.

=back

=head1 LICENSE

This program is licensed under the same terms as Perl.

=cut