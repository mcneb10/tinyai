# tinyai_independent

This is a fork of tinyai (compact implementation of neat) that is "independent" meaning it has no outgoing calls to things like the file system
I wrote this so I can use it on my esp8266

# Modes

If you define the *INDEPENDENT* macro, the library will be compiled without any filesystem calls.
All functions are implemented in *INDEPENDENT* mode

# Bonus!

As an added bonus, i'll add in a script to convert a gene to DOT absolutely free! (coming soon)