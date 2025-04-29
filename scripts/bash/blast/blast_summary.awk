#!/bin/bash

awk '
BEGIN {
    total_count = 0;
    hits_count = 0;
    }
    
# Function to print the buffered block if it has hits
function flush_buffer() {
    if (buffer !~ /No hits found/) {
    	hits_count += 1;
        #print buffer;
    }
    total_count += 1
    buffer = "";  # Clear the buffer
}

/^Query=/ {
    flush_buffer();         # Handle the previous entry
    buffer = $0 "\n";       # Start new buffer with current line
    next;
}

{
    buffer = buffer $0 "\n";  # Add all other lines to buffer
}

END {
    flush_buffer();  # Flush the last entry
    print total_count;
    print hits_count;
}
' $1
