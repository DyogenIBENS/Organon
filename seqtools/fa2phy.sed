#!/bin/sed -nf

:readseq
${H;b out}              # if last line, append to hold, and goto 'out'
1{h;n;b readseq}        # if first, overwrite hold, and start again at 'readseq'
/^>/!{H; n; b readseq}  # if not a sequence label, append to hold, read next line, start again at 'readseq'. Else, it continues to 'out'

:out
x;         # exchange hold content with pattern content
s/^>//;    # substitute the starting '>'
s/\n/  /;  # substitute first newline with 2 spaces
s/\n//g;   # delete any subsequent newline
#s/\n/  /g;
p;         # print pattern buffer

