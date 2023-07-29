package NeuralNet;

use strict;
use warnings;
use Carp;

# Unfortunately, i'll have to wait a while for Corinna to be fully implemented

#use v5.38;
#use feature 'class';

#no warnings 'experimental::class';

# In the mean time i'll have to use this

use v5.26;

use Object::Pad;

# This is probably really inefficient but this is an added bonus so you can't complain
# This is script has been rewritten from script in OOP so it can be prettier

=head1 NAME

B<TinyAI::NeuralNet> - Perl representation of the tinyann.hpp neural network from TinyAI

=head2 SYNOPSIS

This package contains classes that represent part of the tinyann neural network.

=head2 CLASSES

=head3 C<Connection>

Represents a single connection between two C<Neuron>s

=cut

class Connection {
=head4 FIELDS

=over 4

=item C<from> => int, required in constructor, reader

The index of the node the connection is going from

=item C<to> => int, required in constructor, reader

The index of the node the connection is going to

=item C<weight> => float, required in constructor, reader

The weight of the connection

=back

=cut
    field $from   :param :reader;
    field $to     :param :reader;
    field $weight :param :reader;
}

=head3 C<Neuron>

Represents a single node/neuron in the neural network

=cut

class Neuron {
=head4 FIELDS

=over 4

=item C<type> => int, required in constructor, reader

The type of the node. Possible values:

=over 4

=item C<0>

A regular node

=item C<1>

An input

=item C<2>

An output

=item C<3>

A bias node. Never seen one in the wild and i'm not quite sure what they do.

=back

=item C<connections_from> => C<Connection> array, writer, reader

The connections going from the node

=item C<connections_to> => C<Connection> array, writer, reader

The connections going to the node

=item C<prune_if_pruning_enabled> => boolean, writer, reader

Should the node be hidden if pruning useless nodes is enabled in C<dump_to_graphviz> in C<NeuralNet>.
This is true if the node isn't a Input, Output, or Bias node and the node doesn't connect anywhere
or isn't connected to.

=item C<activation> => string, writer, reader

The activation function name. Only set if the file was made by a program compiled in CHANGEABLE_ACTIVATION_AND_AGGREGATION
mode and the node isn't an input/bias? node

=item C<aggregation> => string, writer, reader

The aggregation function name. Only set if the file was made by a program compiled in CHANGEABLE_ACTIVATION_AND_AGGREGATION
mode and the node isn't an input/bias? node

=item C<already_visited> => bool, reader

Indicates if the node has been visited by C<visit> and had C<prune_if_pruning_enabled> set yet

=back

=cut

    field $type                     :param  :reader;
    field @connections_from         :writer :reader = ();
    field @connections_to           :writer :reader = ();
    field $prune_if_pruning_enabled :writer :reader = 1;
    field $activation               :writer :reader = "";
    field $aggregation              :writer :reader = "";
    field $already_visited          :reader = 0;

=head4 METHODS

=over 4

=item C<visit>(C<nodes>: array of C<Neuron>s)

Looks at the connections of this C<Neuron> and determines if it's useless and sets C<prune_if_pruning_enabled>
accoringly. Takes an array of C<Neuron>s from the neural network as it's parameter.
Returns the value of C<prune_if_pruning_enabled>.

=back

=cut

    method visit(@nodes) {
        $already_visited = 1;
        if ( $type > 0 ) {
            $prune_if_pruning_enabled = 0;
        }
        else {
            if ( !@connections_from || !@connections_to ) {
                $prune_if_pruning_enabled = 1;
            }
            else {
                for my $node_from (@connections_from) {
                    if ( $nodes[ $node_from->to ]->already_visited ) {
                        if ( $nodes[ $node_from->to ]
                            ->prune_if_pruning_enabled == 0 )
                        {
                            $prune_if_pruning_enabled = 0;
                            return 0;
                        }
                    }
                    else {
                        if ( $nodes[ $node_from->to ]->visit(@nodes) == 0 ) {
                            $prune_if_pruning_enabled = 0;
                            return 0;
                        }
                    }
                }
            }
        }
        return $prune_if_pruning_enabled;
    }
}

=head3 C<NeuralNet>

A class the represents the Neural Network from tinyann.hpp

=cut

class NeuralNet {
=head4 FIELDS

=over 4

=item C<recurrent> => bool, reader

Specifies if the Neural Network is recurrent.

=item C<changeable_funcs> => bool, reader

Specifies if the Neural Network has changeable activation and aggregation functions

=item C<nodes> => array of C<Neuron>s, reader

The neurons/nodes in the Neural Network

=item C<connections> => array of C<Connection>s, reader

The connections between nodes in the Neural Network

=item C<data> => string, required in constructor

The text data from the neural network file, to be "deserialized" into the fields in the C<NeuralNetwork> class

=back

=cut
    field $recurrent :reader;
    field $changeable_funcs :reader = 0;

    field @nodes :reader = ();
    field $_node_count; # Not used
    field @connections :reader = ();

    field $data :param;

    ADJUST {
        # Get signature
        if ( $data =~ s/^(\S+)\s+// ) {
            my $sig = $1;
            croak "Invalid neural net signature" unless $sig =~ /recurrent/i;
            $recurrent = $sig !~ /^non/i;
            if ( $data =~ s/^([a-zA-Z_]+)\s+// ) {
                if ( $1 eq "CHANGEABLE_ACTIVATION_AND_AGGREGATION" ) {
                    $changeable_funcs = 1;
                }
                else {
                    croak "Unknown data after neural net signature";
                }
            }

            # Get neuron count
            $data =~ s/^(\d+)\s+//;
            $_node_count = $1;

            # Get each neuron
            my $neuron_c = 0;
            while ( $data =~ s/^(\d+)\s+(\d+)\s*// ) {
                my ( $type, $input_count ) = ( $1, $2 );

                my $neuron = Neuron->new( type => $type );

                my @connections_to_node = ();

                if ($changeable_funcs) {
                    $data =~ s/^(\S+)\s+(\S+)\s*//;
                    $neuron->set_activation($2);
                    $neuron->set_aggregation($1);
                }

                for ( 1 .. $input_count ) {
                    $data =~ s/^(\d+)\s+([0-9.-]+)\s*//;
                    my $con = Connection->new(
                        from   => $1,
                        to     => $neuron_c,
                        weight => $2
                    );
                    push @connections,         $con;
                    push @connections_to_node, scalar(@connections) - 1;
                }

                $neuron->set_connections_to(@connections_to_node);

                push @nodes, $neuron;

                $neuron_c++;
            }

            for my $con (@connections) {
                my @cf = $nodes[ $con->from ]->connections_from;
                push @cf, $con;
                $nodes[ $con->from ]->set_connections_from(@cf);
            }

            for my $node (@nodes) {
                if ( !$node->already_visited ) {
                    $node->visit(@nodes);
                }
            }

        }
        else {
            croak "Failed to parse file";
        }
    }

=head4 METHODS

=over 4

=item C<dump_to_graphviz>(C<prune_unused_nodes>: boolean, optional)

Dump the neural network do a B<DOT> digraph.
If C<prune_unused_nodes> is specified and set to 1,
useless nodes won't be included in the graph.
Returns the generated B<DOT> code.

=back

=cut

    method dump_to_graphviz( $prune_unused_nodes = 0 ) {
        my $result = "digraph {\n";

        my ( $input_count, $output_count, $bias_count ) = ( 0, 0, 0 );
        for ( my $i = 0 ; $i < scalar(@nodes) ; $i++ ) {
            my $node = $nodes[$i];
            $result .= $i;
            if ( $node->type == 1 ) {

                # Input
                $result .=
"[style=\"filled\",shape=\"box\",label=\"Input $input_count\"]";
                $input_count++;
            }
            elsif ( $node->type == 2 ) {

                # Output
                $result .=
"[style=\"filled\",shape=\"box\",label=\"Output $output_count\"";
                if ($changeable_funcs) {
                    $result .=
                        ",xlabel=\"Activation: "
                      . $node->activation
                      . "\\nAggregation: "
                      . $node->aggregation . "\"";
                }
                $result .= "]";
                $output_count++;
            }
            elsif ( $node->type == 3 ) {

                # Bias
                # Never seen one of these nodes in the wild before
                $result .=
"[style=\"filled\",shape=\"doublecircle\",label=\"Bias $bias_count\"]";
                $bias_count++;
            }
            else {
                $result .= "[";
                if (   $node->prune_if_pruning_enabled
                    && $prune_unused_nodes )
                {
                    $result .= "style=\"invis\"";
                    $result .= "," if $changeable_funcs;
                }
                if ($changeable_funcs) {
                    $result .=
                        "xlabel=\"Activation: "
                      . $node->activation
                      . "\\nAggregation: "
                      . $node->aggregation . "\"";
                }
                $result .= "]";
            }
            $result .= ";\n";
        }

        for my $con (@connections) {
            next
              if $prune_unused_nodes
              && ( $nodes[ $con->to ]->prune_if_pruning_enabled
                || $nodes[ $con->from ]->prune_if_pruning_enabled );

            $result .=
                $con->from . "->"
              . $con->to
              . "[label=\"Weight "
              . $con->weight
              . "\",penwidth=\""
              . ( 0.1 + abs( $con->weight ) )
              . "\"];\n";
        }

        $result .= "}\n";

        return $result;
    }
}

1;

__END__

=head1 LICENSE

This program is licensed under the same terms as Perl.

=cut
