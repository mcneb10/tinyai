use strict;
use warnings;
use Carp;
use File::Slurper 'read_text';
use List::Util 'first';

# This is probably really inefficient but this is an added bonus so you can't complain

croak "Usage: $0 [files]" if ( @ARGV < 1 );

# TODO!!
my $delete__useless;

my (
    $result,     $data,           %ann,
    $recurrent,  $neuron_count,   $type,
    $input_size, $connected_node, $weight
);
my ( @nodes, @input_nodes, @output_nodes, @bias_nodes, @nodes_with_used_output,
    @nodes_that_take_input );

for (@ARGV) {
    $data = read_text($_);

    # Parse the file
    $result .= "digraph{\n";

    # Trim whitespace
    $data =~ s/^\s+//;

    # Get if recurrent
    $data =~ s/
    ^(\S+) # Start at the beginning and get if recurrent
    \s+ # Delete white space
    (\S+) # Get neuron count or if this has its own activation and aggregation functions
    \s+ # Ignore white space
    //x;
    
    ( $recurrent, $neuron_count ) = ( $1, $2 );

    # TODO: figure out what recurrent means lc($tmp) eq 'recurrent';
    for ( 0 .. $neuron_count - 1 ) {
        $data =~ s/
        ^(\d+) # Get number
        \s+ # Remove space
        (\d+) # Get number
        \s+ # Remove space
        //x;
        ( $type, $input_size ) = ( $1, $2 );
        my $i = $_;
        if ( $type == 1 ) {
            push @input_nodes, $_;
        }
        elsif ( $type == 2 ) {
            push @output_nodes, $_;
        }
        elsif ( $type == 3 ) {
            push @bias_nodes, $_;
        }
        my @connected_nodes;
        if ( $input_size > 0 ) {
            for ( 0 .. $input_size - 1 ) {
                $data =~ s/
                ^(\d+) # Get number
                \s+ # Remove space
                ([0-9\-\.]+) # Get floating point number
                \s* # Remove space if any is left
                //x;
                ( $connected_node, $weight ) = ( $1, $2 );
                push @connected_nodes,
                  { connected_node => $connected_node, weight => $weight };
                push @nodes_with_used_output, $i;
                # Will include duplicates, and it doesn't matter
                push @nodes_that_take_input, $connected_node;
            }
        }
        push @nodes, { type => $type, connected_nodes => \@connected_nodes };
    }

    for ( my $i = 0 ; $i < @nodes ; $i++ ) {
        my $n = $nodes[$i];
        $result .= "$i";
        if ( $n->{type} == 1 ) {
            $result .=
              "[shape=box,label=\"Input $i\",style=filled,fillcolor=lightgrey]";
        }
        elsif ( $n->{type} == 2 ) {
            my $idx = grep { $output_nodes[$_] == $i } 0 .. $#output_nodes;
            $result .=
"[shape=box,label=\"Output $idx\",style=filled,fillcolor=lightgrey]";
        }
        elsif ( $n->{type} == 3 ) {
            my $idx = grep { $bias_nodes[$_] == $i } 0 .. $#bias_nodes;
            $result .= "[shape=box,label=\"Bias $idx\"]";
        }
        else {
            unless ( grep( /^$i$/, @nodes_with_used_output ) != 0
                || $n->{type} == 1
                || $n->{type} == 2 )
            {
                # Make invisible
                $result .= "[style=filled,fillcolor=yellow,label=\"ZERO\"]";
            }
            unless ( grep( /^$i$/, @nodes_that_take_input ) != 0
                || $n->{type} == 1
                || $n->{type} == 2 )
            {
                # Make invisible
                $result .= "[style=filled,fillcolor=red,label=\"USELESS\"]";
            }
        }
        $result .= ";\n";
        for ( @{ $n->{connected_nodes} } ) {
            my $width = 0.1 + abs( $_->{weight} );
            my $c     = $_->{connected_node};
            $result .=
              "$c->$i [label=\"weight " . $_->{weight} . "\",penwidth=$width";
            unless ( grep( /^$c$/, @nodes_with_used_output ) != 0
                || grep( /^$c$/, @input_nodes ) != 0 )
            {
                $result .= ",style=dotted";
            }
            unless ( grep( /^$c$/, @nodes_that_take_input ) != 0
                || grep( /^$c$/, @output_nodes ) != 0 )
            {
                $result .= ",style=dotted";
            }
            $result .= "];\n";
        }

    }
    $result .= "}";
    print "$result\n\n\n";
}
