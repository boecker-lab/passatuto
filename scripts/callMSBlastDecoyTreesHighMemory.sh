#!/bin/bash

export MALLOC_ARENA_MAX=4
java -Xmx35G -classpath /home/i3scke/MSBlast/MSBlast/lib/*:/home/i3scke/MSBlast/MSBlast/bin/ de.unijena.bioinf.msblast.MSBlastDecoyTrees $@

