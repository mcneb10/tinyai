# tinyai_independent

This is a fork of tinyai (compact implementation of neat) that can be "independent", meaning it has no outgoing calls to things like the file system.

I wrote this so I can use it on my esp8266

# Modes

If you define the *INDEPENDENT* macro, the library will be compiled without any filesystem calls.
All functions are implemented in *INDEPENDENT* mode

The macro *CHANGEABLE_ACTIVATION_AND_AGGREGATION* allows programs to change the activation and aggregation functions.

*Note*: when initializing ann::neuralnet in *CHANGEABLE_ACTIVATION_AND_AGGREGATION* make sure to
initialize activation_funcs or aggregation_funcs, even if you are using the default list in the
neat::pool object. Otherwise things will crash and burn!

# Bonus!

As an added bonus, i'll add in Perl a script to convert a gene to DOT absolutely free! (coming soon)