#include <stdio.h>
#include <string.h>

#include "interface.h"


int main()
{
	printf("running example!");

	enum { kMaxArgs = 64 };
	int argc = 0;
	char *argv[kMaxArgs];

	char commandLine[200] = "f5c index -d ../../test/ecoli_2kb_region/fast5_files/ ../../test/ecoli_2kb_region/reads.fasta";
	// char commandLine[200] = "f5c call-methylation -b ../../test/ecoli_2kb_region/reads.sorted.bam -g ../../test/ecoli_2kb_region/draft.fa -r ../../test/ecoli_2kb_region/reads.fasta --secondary=yes --min-mapq=0 -B 2M > ../../test/ecoli_2kb_region/result.txt";
	// char commandLine[500] = "f5c eventalign -b ../../test/ecoli_2kb_region/reads.sorted.bam -g ../../test/ecoli_2kb_region/draft.fa -r ../../test/ecoli_2kb_region/reads.fasta --secondary=yes --min-mapq=0 -B 2M > ../../test/ecoli_2kb_region/f5c_event_align.txt";
	// char commandLine[200] = "f5c meth-freq -d ../../test/ecoli_2kb_region/fast5_files/ ../../test/ecoli_2kb_region/reads.fasta";
	char *p2 = strtok(commandLine, " ");
	
	while (p2 && argc < kMaxArgs-1)
	  {
	    argv[argc++] = p2;
	    p2 = strtok(0, " ");
	  }
	argv[argc] = 0;
	
	init(argc,argv);

	return 0;
}

