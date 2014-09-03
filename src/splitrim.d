/* 

        Last working version: v1j3e5

        BUGS:
        (1) readCounter in trimEntry() is not working (see first 2 fragments in: /home/traceyf/chienchi/working/largeTestset_noEmptyLines.fastq)

        // TO DO: 
        (1) Read in up to 16 entries (64 lines) or more per thread, and THEN process, rather than one entry (4 lines) at a time.
        (2) Stagger the writes for each thread so they don't naturally overlap
        (3) Use lustre's "lctl" command to spread a FASTQ file across X inodes

        Author:         Tracey Freitas
        Last Update:    2013-01-30
        Last Version:   v1f [single-end version]
        Language:       D

        2013-02-28:
        Fixed ASCII transcoding error where an extra element in the array was output along with the quality line

        2013-03-07:
        Added fixed-length read trimming (after quality trimming) with the "--recycle" option. 

        2013-03-19:
        Fixed an off-by-one error in the rawSeqNum count displayed in the stats file.
        
                
        Example (50-bp read):                   
            Original read                       =================================N================      50-bp
            Splitrimmed read (N's)              ==================33=============|=======16=======      33-bp, 16-bp
            Splitrimmed read (Quals)            ==============28============|==4=|=======16=======      28-bp,  4-bp, 16-bp
            
            ========
            Option 1: minL=10 (keeps 2 frags)   ==============28============      =======16=======      28-bp, 16-bp
            ========
              (keeps 2 frags, 44-bp)
              
            ========
            Option 2: fixL=10                   ====10====                        ====10====            10-bp, 10-bp, 10-bp
            ========                                      ====10====                        XXXXXX
              (keeps 3 frags, 30-bp)                                XXXXXXXX
            
            ========                       
            Option 3: fixL=10 --recycle         ====10====                        =======16=======      10-bp, 18-bp, 16-bp
            ========                                      ========18========
              (keeps 3 frags, 44-bp)


        2013-04-03:
        Since SAM files clip header names at the first whitespace, the trim START:STOP positions do not show up in the SAM file.
        The trim START:STOP positions were moved from the end of the header, as in
                                                                ****
            @HISEQ:148:C0F5FACXX:1:1307:2735:88518 1:N:0:CGATGT 1:30
        
        and appended to the header name just before the first whitespace so that they will show up in SAM files, as in
                                                   ****
            @HISEQ:148:C0F5FACXX:1:1307:2735:88518:1:30 1:N:0:CGATGT

        This will facilitate tracking of (1) the hit count in terms of the actual full-length reads, and (2) the average fraction
        of a read that does map. The latter would be calculated as follows:
        
        ///////////////////////////////////////////////////////////////////////////////////////////////////
        //                                                                                               //
        //  Average Frac. of Read Mapping = { SUM_i=1|N [(# bases mapped)i/(# bases available)i] } / N   //
        //                                                                                               //
        //-----------------------------------------------------------------------------------------------//
        //                                                                                               //
        //  where                    N = no. of full-length reads                                        //
        //                           i = the i-th full-length read                                       //
        //           (# bases mapped)i = no. of bases in i that are in an alignment                      //
        //        (# bases available)i = no. of bases in i that passed the trimming step                 //
        //                                                                                               //
        ///////////////////////////////////////////////////////////////////////////////////////////////////
        
        2013-04-04:
        Reverted the positions of the trim START:STOP positions to the very end of the header line. 
        
        All FASTQ reads are converted to the OLD Illumina format, where "/1" or "/2" is appended to the header line
        just before the first whitespace. Since only single-end reads are considered here, only "/1" is used.
        -> Uniqueness of each fragment derived from a single, full-length read, is enforced by prepending a "_" and
           a counter just before the "/1" designator.

            @HISEQ:148:C0F5FACXX:1:1307:2735:88518:1:30 1:N:0:CGATGT
        to  @HISEQ:148:C0F5FACXX:1:1307:2735:88518 1:N:0:CGATGT 1:30
        to  @HISEQ:148:C0F5FACXX:1:1307:2735:88518/1 1:N:0:CGATGT 1:30
        
        But in a SAM file, the first whitespace-delimited header string is indistinguishable from
            @HISEQ:148:C0F5FACXX:1:1307:2735:88518/1 1:N:0:CGATGT 31:60

        So we prepend a fragment counter before the "/1" designator:
            __________SAM-recognized_Header___________
            @HISEQ:148:C0F5FACXX:1:1307:2735:88518_0/1 1:N:0:CGATGT 1:30
            @HISEQ:148:C0F5FACXX:1:1307:2735:88518_1/1 1:N:0:CGATGT 31:60

        

        TO DO:
        =====
        -> add option to trim N's off the ends of a read only, without any splitting

*/


import std.stdio;
import std.string;
import std.getopt;
import std.path: buildPath, baseName, stripExtension;
import std.algorithm;
import std.math;
import std.c.stdlib;
import std.exception;
import std.datetime;
import std.file;
import std.range;
import std.regex;
import std.parallelism;
import std.conv;

struct seqLenStats {
    double  trimStdev;
    double  trimAvg;
    ushort  trimMax;
    ushort  trimMin;
    double  trimMedian;
};
struct inputOptions {
    ///////////////////////////////////
    // General
    ///////////////////////////////////
    string             progName;
    string             progVersion;
    string             inFastq;
    string             inFastqBasename;

    ///////////////////////////////////
    // Trimming & Filtering
    ///////////////////////////////////
    ushort             userMinL;                // min. length filter
    ushort             userFixL;                // fixed length filter
    bool               recycle;                 // true or false
    uint               threads;       
    ushort             baseLength;              // the non-zero contents of either userMinL or userFixL
    ubyte              userMinQ;                // min. quality fileter
    ushort             userMaxN;                // min. count filter
    float              lowComplexity;           // 0.0 - 1.0
    float              minAveQ;

    ///////////////////////////////////
    // Transcoding
    ///////////////////////////////////
    ubyte              asciiEncoding;           // 33 or 64
    ubyte              outEncoding;
    byte               outOffset;

    ///////////////////////////////////
    // Multithreading
    ///////////////////////////////////
    bool               newIlluminaFastqFormat;
    ulong              inFastqFileSize;         // measured in bytes
    ulong              fastqEntrySize;          // # of bytes in a valid 4-lined FASTQ entry
    ulong              fastqEntries;            // inFastqFileSize/fastqentrySize
    ulong              fastqEntriesPerThread;   // fastqEntries/threads
    float              fastqExactEntries;       // inFastqFileSize/fastqEntrySize
    uint               effThreadCount;          // No. of effective threads after partitioning FASTQ file
    ulong[]            partitionOffsets;        // Contains the byte offsets of each partition
    ulong              flush;                   // When read count hits this value, trimmed entries are flushed to disk

    ///////////////////////////////////
    // Output
    ///////////////////////////////////
    string             outPath;
    string             prefix;
    string             outFastq;
    string             outLenFilename;
    string             outGCFilename;
    string             statsFilename;
    string             outStatsFilename;
    string             outHistoFilename;
    bool               sortLenAsc;
    bool               sortLenDesc;
    bool               verbose;

    ///////////////////////////////////
    // Statistics
    ///////////////////////////////////
    ulong              rawSeqNum;               // No. of input reads
    ulong              totalRawSeqLen;          // Total input read length (bp)
    ulong              trimSeqNum;              // No. of trimmed reads
    ulong              totalTrimSeqLen;         // Total trimmed read length (bp)
    ulong[ushort]      trimSeqLenHash;          // Freq[Length]
    ulong              trashedBases;            // No. of bases trashed from using --fixL without --recycle
};

__gshared   inputOptions    inOpts;

immutable   byte    qBest  = 41;
immutable   string  PREFIX = "QC";

//
// PERFORMANCE:
//
//      minL=1      maxTries=100_000        TIME: 19475 ms      STDOUT
//      minL=1      maxTries=100_000        TIME:   827 ms      FILE
//      minL=12     maxTries=100_000        TIME: 14766 ms      STDOUT
//      minL=12     maxTries=100_000        TIME:   296 ms      FILE
//      minL=50     maxTries=100_000        TIME:   156 ms      FILE
//      minL=50     maxTries=1_000_000      TIME:  1524 ms      FILE
//      minL=50     maxTries=10_000_000     TIME: 15504 ms      FILE
//
//      minL=3      maxTries=1_000_000      TIME: 10049 ms      FILE        empty while loop
//      minL=3      maxTries=1_000_000      TIME: 10295 ms      FILE        i in while loop
//      minL=3      maxTries=10_000_000     TIME: 102898 ms     FILE        empty while loop
//      minL=3      maxTries=10_000_000     TIME: 105442 ms     FILE        i in while loop
    //ushort minQ    = 20;
    //ushort minL    = 3;
    //ulong maxTries = 1; //10_000_000;

/**********************************************************************/
void verifyFastqFormat(std.stdio.File IN, ref bool newIlluminaFormat) {

    ubyte[char] stdNucleotides = [ 'A':1, 'C':1, 'G':1, 'T':1, 'N':1, 'a':1, 'c':1, 'g':1, 't':1, 'n':1 ];
    enum ctr = ctRegex!(`^(\@.+\/\d+)\s(.+)`);          // Looking for /1 or /2 in OLD Illumina formats

    int lineCount = 0;
    while(!IN.eof()) {
        string line = IN.readln().chomp;                // header
        if(line.length == 0) continue;                  // Skip any leading empty lines
        ++lineCount;
        if(lineCount == 1) {                            // Check that the header begins with "@"
            if(line[0] == '@') {    // valid header line
                auto m = match(line, ctr);              // Old Illumina formats terminate their
                string header = m.captures[1];          // headers with either a /1 or /2 for
                                                        // paired-end reads, while single-end 
                                                        // reads terminate with /1
                if(header == "") {  // NEW Illumina format
                    newIlluminaFormat = true;
                }
            }
            else {                  // invalid
                writeln("*FATAL*: Unrecognized FASTQ format!"); 
                exit(1);
            }
            
        }
        if(lineCount == 2) {                            // 2nd valid line better be DNA!
            foreach(base; line) {
                //writeln("Base: ", base);
                if(base !in stdNucleotides) {
                    writeln("*FATAL*: Unrecognized FASTQ format!"); 
                    exit(1);
                }
            }
        }
        if((lineCount == 3) && (line[0] != '+')) {      // 3rd valid line should start with "+"
            writeln("*FATAL*: Unrecognized FASTQ format!"); 
            exit(1);
        }
        if(lineCount == 4) {                            // 4th line is the quality line
            verifyAsciiEncoding(inOpts, line);
        }
    }

    IN.rewind();

}
/**********************************************************************/
bool calcStats(ref seqLenStats myStats, ref ulong[ushort] trimSeqLenHash) {

    if(trimSeqLenHash.length == 0) {
        myStats.trimStdev = 0.0;
        myStats.trimAvg   = 0.0;
        myStats.trimMax   = 0;
        myStats.trimMin   = 0;
        myStats.trimMedian= 0.0;
        return false;
    }

    //foreach(readLen; sort!((a,b) {return a > b; }) (trimSeqLenHash.keys)) {
    //    writeln(readLen," => ", trimSeqLenHash[readLen]);
    //}
    
    //auto max = trimSeqLenArray[$-1];
    //auto min = trimSeqLenArray[0];

    //double median;
    //if(trimSeqLenArray.length & 1) {	// ODD
    //    median = trimSeqLenArray[cast(ulong)(trimSeqLenArray.length/2)];
    //}
    //else {
    //    median = (trimSeqLenArray[trimSeqLenArray.length/2] 
    //           +  trimSeqLenArray[trimSeqLenArray.length/2 - 1]) / 2;
    //}

    // Step1: Find mean
    ulong total1   = 0;
    ulong numElems = 0;
    foreach(ushort len; trimSeqLenHash.keys) { 
        total1 += len*trimSeqLenHash[len];
        numElems += trimSeqLenHash[len];
    }
    
    // EDIT
    //double mean1 = total1/trimSeqLenArray.length;
    double mean1 = total1/numElems;

    // Step2: Find sum of squares of difference between mean & datapoint
    double total2 = 0;
    foreach(len; trimSeqLenHash.keys) {
        total2 += ((mean1-len)^^2) * trimSeqLenHash[len];
    }
    //double mean2 = total2/trimSeqLenArray.length;
    double mean2 = total2/numElems;

    // Step3: Calc. stdev
    double stdev = sqrt(mean2);

    myStats.trimStdev = stdev; //writefln("stdev: %.2f", myStats.trimStdev);
    myStats.trimAvg   = mean1;
    //myStats.trimMax   = max;
    //myStats.trimMin   = min;
    //myStats.trimMedian= median;

    return true;

}
/**********************************************************************/
void processOptions1(ref inputOptions inOpts) {

    // Check that prefix is valid
    if(inOpts.prefix.length == 0) {
        inOpts.prefix = PREFIX;
    }


    // Input file exists?
    if(!inOpts.inFastq.exists()) {
        writeln("*FATAL*: Input FASTQ file \"",inOpts.inFastq,"\" does not exist!");
        exit(0);
    }
    
    // Input file is a file?
    if(!inOpts.inFastq.isFile) {
        writeln("*FATAL*: Input file \"",inOpts.inFastq,"\" is not a file");
        exit(0);
    }

    // Output path exists?
    if(!inOpts.outPath.exists()) {

        if(inOpts.verbose) {
            write("Output directory \"",inOpts.outPath,"\" does not exist. Attempting to create...");
            stdout.flush();
        }

        // try to create if !exists
        try {
            mkdirRecurse(inOpts.outPath);
        }
        catch (Exception e) {
            //throw new FileException(format("Cannot create directory %s.", inOpts.outPath));
            //writeln(e);
            if(inOpts.verbose) writeln("FAILED");
            writeln("*FATAL*: Unable to create output directory: ", e.msg);
            exit(0);
            //stdin.readln();
        }
        if(inOpts.verbose) writeln("successfully created.");

    }
    
    if(!inOpts.outPath.isDir) {
        writeln("*FATAL*: Output path \"",inOpts.outPath,"\" is not a directory.");
        exit(0);
    }

    if(!inOpts.userMinL && !inOpts.userFixL) {
        writeln("One of either \"--minL=<INT>\" or \"--fixL=<INT>\" must be specified.");
        exit(0);
    }

    if((inOpts.userMinL > 0) && (inOpts.userFixL > 0)) {
        writeln("Please specify only one of \"--minL\" or \"--fixL\". Abort.");
        exit(0);
    }
    
    if(inOpts.threads < 1) {
        writeln("*FATAL*: Thread count [",inOpts.threads,"] must be >= 1. Abort.");
        exit(0);
    }
   
}
/**********************************************************************/
void processOptions2(inputOptions inOpts) {

    enforce((inOpts.asciiEncoding == 33) || (inOpts.asciiEncoding == 64), "ASCII Encoding must be 33 or 64!");
    enforce(inOpts.outPath.length > 0,                                    "Invalid output path");
    enforce((inOpts.userMinL > 0) || (inOpts.userFixL > 0),               "Valid read length must be > 0");
    enforce((inOpts.userMinQ >= 0) && (inOpts.userMinQ <= qBest),         "Minimum quality must be between 0-41");
//    enforce(inOpts.userMaxN);
    enforce((inOpts.sortLenAsc & inOpts.sortLenDesc) == false,            "Cannot specify both --sortLenAsc and --sortLenDesc");
    enforce(((inOpts.minAveQ > 0) && (inOpts.minAveQ <= qBest)),          "Minimum average read quality is out of range");
    if(
          (inOpts.userMinL  > 0) 
       && (inOpts.userFixL == 0) 
       && (inOpts.recycle)
      ) 
    writeln("Note: Option \"--recycle\" is meaningless when \"--minL\" is specified. Use \"--fixL\" instead. Ignoring...\n");

}
/**********************************************************************/
void initFiles(ref inputOptions inOpts, ref std.stdio.File IN, 
               ref std.stdio.File OUT,  ref std.stdio.File STATS, 
               ref std.stdio.File HISTO) {

    debug writeln("outStatsFilename=",inOpts.outStatsFilename);
    debug writeln("outHistoFilename=",inOpts.outHistoFilename);

    IN     = std.stdio.File(inOpts.inFastq, "r"); 
    OUT    = std.stdio.File(inOpts.outFastq, "w");
    STATS  = std.stdio.File(inOpts.outStatsFilename, "w");
    HISTO  = std.stdio.File(inOpts.outHistoFilename, "w");
    //GC  = std.stdio.File(inOpts.outGCFilename, "w");
    
}
/**********************************************************************/
void verifyAsciiEncoding(ref inputOptions inOpts, string q) {
    
    // swap encoding if wrong encoding
    bool wrongEncoding = false;
    QUAL: foreach(qual; q) {
    
        //write(qual,"=",cast(int)qual," - 64 =",cast(int)qual-64,"/ ");
    
        // True encoding is 33 if this is negative
        if((cast(int)qual - 64) < 0) {
            // It's 33, but given encoding is 64; swap
            if (inOpts.asciiEncoding == 64) { 
                debug writeln("Wrong encoding (64). Switching to 33!");
                wrongEncoding = true; 
                inOpts.asciiEncoding = 33;
                break QUAL; 
            }
            else {
                debug writeln("Encoding is 33...");
                break QUAL;     // Given encoding is the true encoding, since a
                                // Q64 encoding - 64 will be >= 0;
            }
        }
        // True encoding is 64 if this is > 10
        else if((cast(int)qual - 64) > 10) {
            if(inOpts.asciiEncoding == 33) {
                debug writeln("Wrong encoding (33). Switching to 64!");
                inOpts.asciiEncoding = 64;
                break QUAL;
            }
            else {
                debug writeln("Encoding is 64...");
                break QUAL;
            }
        }
        
    }
    
    if((inOpts.asciiEncoding == 33) && (inOpts.outEncoding == 64)) 
        inOpts.outOffset = (64-33);
    else if((inOpts.asciiEncoding == 64) && (inOpts.outEncoding == 33))
        inOpts.outOffset = -(64-33);

//    write("Press <ENTER>..."); stdin.readln();

}
/**********************************************************************/
void writeHistoFile(std.stdio.File HISTO, in inputOptions inOpts) {


    if(inOpts.sortLenAsc) {
        // sort ascending
        foreach(len; sort!((a,b) { return b > a; }) (inOpts.trimSeqLenHash.keys))
            HISTO.writeln(len,"\t",inOpts.trimSeqLenHash[len]);    
        return;
    }

    if(inOpts.sortLenDesc) {
        // sort descending
        foreach(len; sort!((a,b) { return a > b; }) (inOpts.trimSeqLenHash.keys))
            HISTO.writeln(len,"\t",inOpts.trimSeqLenHash[len]);    
        return;
    }
    
    // unsorted
    foreach(len, freq; inOpts.trimSeqLenHash) 
        HISTO.writeln(len,"\t",freq);

}
/**********************************************************************/
void writeStatsFile(std.stdio.File STATS,       in inputOptions inOpts,
                    in seqLenStats myStats,     in ulong programTime,
                    in string[] args) {

    STATS.write  ("\nCMD: ");

    foreach(arg; args) STATS.write(arg," ");
    STATS.writeln();

    STATS.writeln();
    STATS.writeln("                 \tRAW\tSPLIT-TRIMMED");
    STATS.writeln("                 \t===\t=============");
    STATS.write  ("      # of Reads:\t",inOpts.rawSeqNum,     "\t",inOpts.trimSeqNum,"\t(");
    STATS.writef ("%.2f", cast(double)inOpts.trimSeqNum/inOpts.rawSeqNum*100);
    STATS.writeln(" %)");

    STATS.write  ("      # of Bases:\t",inOpts.totalRawSeqLen,"\t",inOpts.totalTrimSeqLen,"\t(");
    STATS.writef ("%.2f", cast(double)inOpts.totalTrimSeqLen/inOpts.totalRawSeqLen*100);
    STATS.writeln(" %)");
    
    STATS.write  ("Mean Read Length:\t");
    double rawMeanReadLen = cast(double) inOpts.totalRawSeqLen/inOpts.rawSeqNum;
    STATS.writef ("%.2f", rawMeanReadLen);
    STATS.write  ("\t");
    STATS.writef ("%.2f", myStats.trimAvg);
    STATS.write  ("\t(");
    STATS.writef ("%.2f", cast(double)myStats.trimAvg/rawMeanReadLen*100);
    STATS.writeln(" %)");
    STATS.write  ("           stdev:\t0.00\t");
    STATS.writef ("%.2f", myStats.trimStdev);
    STATS.writeln();
    STATS.writeln("---------------------------------------------");
    STATS.writeln("   Trashed Bases:\t",inOpts.trashedBases);
    STATS.writeln(" Total Time (ms):\t",programTime);


/*


    STATS.writeln( "\nBefore Trimming");
    STATS.writeln( "Reads #:\t", rawSeqNum);
    STATS.writeln( "Total bases:\t", totalRawSeqLen);
    STATS.write(   "Mean Read Length:\t");
    STATS.writefln("%.2f", cast(float)totalRawSeqLen/rawSeqNum);
    STATS.writeln( "\nAfter Trimming");
    STATS.write(   "Reads #:\t", trimSeqNum,"\t(");
    STATS.writef(  "%.2f", cast(float)trimSeqNum/rawSeqNum*100);
    STATS.writeln( " %)");
    STATS.write(   "Total bases:\t",totalTrimSeqLen,"\t(");
    STATS.writef(  "%.2f", cast(float)totalTrimSeqLen/totalRawSeqLen*100);
    STATS.writeln( "%)");
    
    STATS.write(   "Mean Read Length:\t");
    if(trimSeqNum > 0) {
        STATS.writef(  "%.2f", myStats.trimAvg);
        STATS.write(   "\t+/-\t");
        STATS.writefln("%.2f", myStats.trimStdev);
    }
    else {
        STATS.writeln("0.00");
    }
    
    // Edit:
    STATS.writeln();
    STATS.writeln("Total Trim Time (ms):\t", programTime);
*/


/*  Perl Equiv:
    
    if (@paired_files){
      printf $fh ("  Paired Reads #:\t\%d (%.2f %%)\n",$paired_seq_num, $paired_seq_num/$trimmed_num*100);
      printf $fh ("  Paired total bases:\t\%d (%.2f %%)\n",$total_paired_bases,$total_paired_bases/$total_trimmed_seq_len*100);
      printf $fh ("  Unpaired Reads #:\t\%d (%.2f %%)\n", $trimmed_num - $paired_seq_num, ($trimmed_num - $paired_seq_num)/$trimmed_num*100);
      printf $fh ("  Unpaired total bases:\t\%d (%.2f %%)\n", $total_trimmed_seq_len - $total_paired_bases , ($total_trimmed_seq_len - $total_paired_bases)/$total_trimmed_seq_len*100);
    }
    
    printf $fh ("\nDiscarded reads #:\t\%d (%.2f %%)\n", $total_num - $trimmed_num , ($total_num - $trimmed_num)/$total_num*100);
    printf $fh ("Trimmed bases:\t\%d (%.2f %%)\n", $total_raw_seq_len - $total_trimmed_seq_len, ($total_raw_seq_len - $total_trimmed_seq_len)/$total_raw_seq_len*100);
*/    
    
}
/******************************************************************************/
//void splitFilter(in string oligo, in ushort[] quals, std.stdio.File OUT) {
/*
void splitFilter(in string oligo, in ushort[] quals) {

    ushort i = 0;
    ushort startPos = 0;
    ushort stopPos  = 0;

    while (i < oligo.length) {
        if(quals[i] >= minQ) {
            startPos = i;
            i++;
            while((quals[i] >= minQ) && (i < oligo.length)) {
                i++;
            }
           
            stopPos = i;
            
            // Is fragment of sufficient length to attempt to split on N's?
            if((stopPos-startPos) >= minL) {
                // Split on N's: Save those that are of sufficent length
                foreach(frag; std.regex.split(oligo[startPos..stopPos], regex("N+"))) {
                    //if(frag.length >= minL) OUT.writeln(frag);
                    if(frag.length >= minL) writeln(frag);
                }
            }
        }
        i++;
    }
}
*/
/******************************************************************************/
void parseFASTQ(in ulong idx, in ulong offsetStart, ref inputOptions inOpts, in ulong stagger, std.stdio.File OUT) { //in string infilename, in ulong fileSize) {
    
    std.stdio.File FASTQ = std.stdio.File(inOpts.inFastq, "r");
    
    // Determine offset of last byte to read in
    ulong chunkEnd = (idx == inOpts.partitionOffsets.length-1)
                   ? (inOpts.inFastqFileSize)
                   : (inOpts.partitionOffsets[idx+1]-1);

    FASTQ.seek(offsetStart, SEEK_SET);
    string h1, s, h2, q;
    ulong currOffset = offsetStart;

    writeln("IDX = ",idx,": reading from ",currOffset," to ",chunkEnd);


    ulong           localRawSeqNum       = 0;
    ulong           localTotalRawSeqLen  = 0;
    ulong           localTrimSeqNum      = 0;
    ulong           localTotalTrimSeqLen = 0;
    ulong           localTrashedBases    = 0;
    ulong[ushort]   localTrimSeqLenHash;

    ulong           localRawSeqCounter   = 0;
    string          trimmedFqStack;

    // Staggering the flushing of trimmed reads to disk so as to prevent direct 
    // competition between threads for the single filehandle.
    do {
        h1 = FASTQ.readln().chomp();
        s  = FASTQ.readln().chomp();
        ++localRawSeqNum;
        localTotalRawSeqLen += s.length;
        h2 = FASTQ.readln();
        q  = FASTQ.readln().chomp();
        
        if(s.length < inOpts.baseLength) continue;
        
        trimEntry(inOpts, h1, s, h2, q, localTrimSeqNum, localTotalTrimSeqLen, localTrimSeqLenHash, localTrashedBases, trimmedFqStack);              // <-------------------- HERE!

    } while ((localRawSeqNum < stagger) && (FASTQ.tell() < chunkEnd));    
    
    writeln("Staggering at ",stagger," reads; processed ",localRawSeqNum," reads");
    flushEntries(localRawSeqCounter, trimmedFqStack, OUT);


    /*  
        FIX: A single emtpy line will hoze this entire program. Implement an
             additional FOR or WHILE loop to screen 4 valid, non-empty lines.
    */

    while (FASTQ.tell() < chunkEnd) {
    //do {
    
        /////////////////////////////////
        // FASTQ entry every 4 lines
        /////////////////////////////////
        h1 = FASTQ.readln().chomp();                    // 1of4: seq header
        //writeln("HEADER:      \"", h1,"\"");
    
        s = FASTQ.readln().chomp();                     // 2of4: sequence
        ++localRawSeqNum;
        localTotalRawSeqLen += s.length;                // total length of all input sequences
        //writeln("SEQUENCE:    \"", s,"\"");

        h2 = FASTQ.readln();                            // 3of4: qual header (no need to chomp)
        //std.stdio.write("DESCRIPTION: ", h2);

        q = FASTQ.readln().chomp();                     // 4of4: qualities
        //writeln("QUALITY:     \"", q,"\"");

        if(s.length < inOpts.baseLength) continue;      // filter out short seq lengths

        /////////////////////////////////
        // TRIM READS AND PIPE TO DISK
        /////////////////////////////////
        // Need to update the global values of: inOpts.trimSeqNum, inOpts.totalTrimSeqLen, inOpts.trimSeqLenHash
        // Output is synchronized to the OUT filehandle so that buffer is flushed after localTrimSeqNum == X
        trimEntry(inOpts, h1, s, h2, q, localTrimSeqNum, localTotalTrimSeqLen, localTrimSeqLenHash, localTrashedBases, trimmedFqStack);              // <-------------------- HERE!

        (++localRawSeqCounter == inOpts.flush) && flushEntries(localRawSeqCounter, trimmedFqStack, OUT);

    } //while (FASTQ.tell() < chunkEnd);

    if(localRawSeqCounter != 0) flushEntries(localRawSeqCounter, trimmedFqStack, OUT);

    writeln("IDX(",idx,") counted ",localRawSeqNum," reads");

    synchronized {            
        inOpts.rawSeqNum       += localRawSeqNum;                           // declare
        inOpts.totalRawSeqLen  += localTotalRawSeqLen;                      // declare
        inOpts.trimSeqNum      += localTrimSeqNum;                          // declare
        inOpts.totalTrimSeqLen += localTotalTrimSeqLen;                     // declare
        inOpts.trashedBases    += localTrashedBases;                        // declare
        foreach (ref elem; localTrimSeqLenHash.keys) {
            inOpts.trimSeqLenHash[cast(ushort)elem] += localTrimSeqLenHash[cast(ushort)elem];       // declare
        }
    }

}
/******************************************************************************/
void flushEntries(ref ulong localRawSeqCounter, ref string trimmedFqStack, std.stdio.File OUT) {

    synchronized {
        //writeln("flushEntries() Called! currCounter: ", localRawSeqCounter);

        // Write trimmedFqStack to disk
        OUT.seek(0, SEEK_END);
        //OUT.writeln(trimmedFqStack);  // If using this, must take off last "\n" from stack
        OUT.rawWrite(trimmedFqStack);

        //foreach (ref trimmedEntry; trimmedFqStack)
        //    OUT.writeln(trimmedEntry);

        trimmedFqStack = [];
        localRawSeqCounter = 0;
    }
}
/******************************************************************************/
void trimEntry(inputOptions inOpts, 
               in string h1,                        in string s, 
               in string h2,                        in string q, 
               ref ulong trimSeqNum,                ref ulong totalTrimSeqLen, 
               ref ulong[ushort] trimSeqLenHash,    ref ulong trashedBases,
               ref string trimmedFqStack) {

    ushort i = 0;
    ushort qualStart = 0;
    ushort qualStop  = 0;
    bool   changingEncoding = (inOpts.asciiEncoding == inOpts.outEncoding)
                            ? false
                            : true;

    //string trimmedEntry;

    ulong fragCounter = 0;

    while (i < s.length) {

        // Decode quality and check its validity        
        if((cast(ushort)q[i] - inOpts.asciiEncoding) >= inOpts.userMinQ) {

            qualStart = i++;

            // LOOP: Decode quality and check its validity
            while(((cast(int)q[i] - inOpts.asciiEncoding) >= inOpts.userMinQ) && (i < s.length)) {
                i++;
            }
           
            qualStop = i;

            // Is fragment of sufficient length to attempt to split on N's?
            if(qualStop - qualStart >= inOpts.baseLength) {

                // Split on N's: Save those that are of sufficent length
                ushort n      = qualStart;
                ushort nStart = qualStart;
                ushort nStop  = qualStop;

                while(n < qualStop) {
                    if(s[n] != 'N') {       // Entering valid stretch
                        nStart = n++;
                        while((s[n] != 'N') && (n < qualStop)) {
                            n++;
                        }
                        nStop = n;

                        ushort fragLen = cast(ushort) (nStop - nStart);

                        //writeln("---------> OUTSIDE LOOP1");

                        // Ensure ATGC-stretch is of sufficient length            
                        if(fragLen >= inOpts.baseLength) {
                            
                            //writeln("---------> INSIDE LOOP1");

                            ////////////////////////////
                            // Append to FASTQ TRIM file
                            ////////////////////////////

                            enum ctr = ctRegex!(`^(\S+)\s?(.*)`);        // generate native machine code for regex at compile time
                            auto m = match(h1, ctr);                    // match and store
                            string h1leader  = m.captures[1];           // extract match #1
                            string h1trailer = m.captures[2];           // extract match #2

                            // Force conversion of h1leader from the NEW ILLUMINA FORMAT to the OLD FORMAT
                            //if(inOpts.newIlluminaFastqFormat) {
                            //    h1leader ~= "/1";
                            //}

                            //writeln("---------> OUTSIDE fixL LOOP");

                            // OPT: Fix-trim read
                            if(inOpts.userFixL) {
    
                                //writeln("---------> INSIDE fixL LOOP");

                                debug writeln("FRAG: ", s[nStart..nStop]);
    
                                float div = cast(float) fragLen / inOpts.userFixL;  // No. of cycles (real)
                                ushort rounds = cast(ushort) div;       // No. of cycles (int)

                                //writeln("    floatRounds = ",div,"\t","ushortRounds = ",rounds);

                                // Need to drop one cycle if keeping overhangs
                                if(inOpts.recycle)
                                    rounds--;

                                //////////////////////////////////////////////////
                                // Trim to fixed length & Append Coords to Header
                                //////////////////////////////////////////////////
                                //writeln("ROUNDS=0..",rounds);

                                //writeln("---------> OUTSIDE rounds LOOP");

                                foreach(j; 0..rounds) {
                                    ushort jStart = cast(ushort) (cast(int)nStart + (j * cast(int)inOpts.userFixL));
                                    ushort jStop  = cast(ushort) (jStart + inOpts.userFixL); // = nStart+((j+1)*inOpts.userFixL);

                                    //writeln("---------> INSIDE rounds LOOP");

                                    // We display "jStart+1" because sequences are 1-offset outside of D
                                    //writeln("....... BEFORE1: ",fragCounter, ".......");
                                    //write(h1leader,"_",fragCounter,"/1 ",h1trailer," ",jStart+1,":",jStop,"\n",s[jStart..jStop],"\n"); //,h2);
                                    trimmedFqStack ~= text(h1leader,"_",fragCounter,"/1 ",h1trailer," ",jStart+1,":",jStop,"\n",s[jStart..jStop],"\n",h2);
                                    if(!changingEncoding) {
                                        trimmedFqStack ~= q[jStart..jStop] ~ "\n";
                                        //writeln(q[jStart..jStop]);
                                    }
                                    else {
                                        foreach(k, ref elem; q[jStart..jStop]) {
                                            trimmedFqStack ~= cast(char) (cast(int)(elem) + inOpts.outOffset);
                                            //write(cast(char) (cast(int)(elem) + inOpts.outOffset));
                                        }
                                        trimmedFqStack ~= "\n";
                                        //writeln();
                                    }
                                    
                                    ++fragCounter;
                                    //writeln("....... AFTER1: ",fragCounter, ".......");

                                }

                                ////////////////////////////
                                // Track trimmed read data
                                ////////////////////////////
                                trimSeqNum += rounds;                       // Add one frag to read list
                                totalTrimSeqLen += rounds*inOpts.userFixL;  // Add length of frag
                                trimSeqLenHash[inOpts.userFixL] += rounds;


                                //writeln("---------> OUTSIDE recycle LOOP");
    
                                ///////////////////////////////////////////////////
                                // OPTIONAL: Process last frag containing overhang                                        
                                ///////////////////////////////////////////////////
                                ushort overhang = cast(ushort) ((div - cast(float)rounds) * cast(float)inOpts.userFixL);
                                if(inOpts.recycle) {

                                    //writeln("---------> INSIDE recycle LOOP");


                                    // Note the difference in range variable names: from jStart to nStop (*NOT* jStart to jStop)
                                    ushort jStart = cast(ushort)(nStart+(rounds*inOpts.userFixL));

                                    // We display "jStart+1" because sequences are 1-offset outside of D
                                    //writeln("....... BEFORE2: ",fragCounter, ".......");
                                    //write(h1leader,"_",fragCounter,"/1 ",h1trailer," ",jStart+1,":",nStop,"\n",s[jStart..nStop],"\n",h2);
                                    trimmedFqStack ~= text(h1leader,"_",fragCounter,"/1 ",h1trailer," ",jStart+1,":",nStop,"\n",s[jStart..nStop],"\n",h2);
                                    if(!changingEncoding) {
                                        trimmedFqStack ~= q[jStart..nStop] ~ "\n";
                                        //writeln(q[jStart..nStop]);
                                    }
                                    else {
                                        foreach(ref elem; q[jStart..nStop]) {
                                            trimmedFqStack ~= cast(char) (cast(int)(elem) + inOpts.outOffset);
                                            //write(cast(char) (cast(int)(elem) + inOpts.outOffset));
                                        }
                                        trimmedFqStack ~= "\n";
                                        //writeln();
                                    }

                                    ++fragCounter;
                                    //writeln("....... AFTER2: ",fragCounter, ".......");

                                    ////////////////////////////
                                    // Track trimmed read data
                                    ////////////////////////////
                                    trimSeqNum++;                       // Add one more fragment!
                                    totalTrimSeqLen += overhang;        // Add length of overhanging frag
                                    trimSeqLenHash[overhang]++;         // Track last read length
                                }
                                else {
                                    ////////////////////////////
                                    // Track trimmed read data
                                    ////////////////////////////
                                    trashedBases += overhang;    // Record no. of trashed bases
                                }
                            }
                            // OPT: Use standard minL instead of fixL
                            else {
                                ////////////////////////////
                                // Append coords to header
                                ////////////////////////////
                                // Append qualStart..qualStop to header to make split frags unique

                                // We display "jStart+1" because sequences are 1-offset outside of D
                                //write(h1leader,"_",fragCounter,"/1 ",h1trailer," ",nStart+1,":",nStop,"\n",s[nStart..nStop],"\n"); //,h2);
                                trimmedFqStack ~= text(h1leader,"_",fragCounter,"/1 ",h1trailer," ",nStart+1,":",nStop,"\n",s[nStart..nStop],"\n",h2);
                                if(!changingEncoding) {
                                    trimmedFqStack ~= q[nStart..nStop] ~ "\n";
                                    //write(q[nStart..nStop]);
                                }
                                else {

                                    ////////////////////////////
                                    // Transcode Qualities
                                    ////////////////////////////
                                    // Transcode entire fragment and select range from this
                                    char[] qConv = new char[](fragLen);
                                    if(changingEncoding) {
                                        foreach (j, ref elem; q[nStart..nStop])
                                            qConv[j] = cast(char)(cast(int)(elem) + inOpts.outOffset);
                                    }
                                    trimmedFqStack ~= qConv ~ "\n";
                                    //write(qConv);
                                }

                                ++fragCounter;

                                ////////////////////////////
                                // Track trimmed read data
                                ////////////////////////////
                                trimSeqNum++;
                                totalTrimSeqLen += (fragLen);
                                trimSeqLenHash[cast(ushort)(fragLen)]++;
                            }
                            
                        } // if(fragLen >= inOpts.baseLength)
                    }
                    n++;
                } //while(n)
            } //if(minL)
        } //if(minQ)
        i++;
    } //while(i)

    //readln();
}
/******************************************************************************/
/*
    Returns: byte offset of next FASTQ entry, or -1 signifying no valid entry found.
*/
ulong getFastqEntryOffset (std.stdio.File FASTQ, in ulong offset) {

    ulong nextOffset = -1;

    char[] buf;
    ulong entry_offset     = 0;
    bool entry_found       = false;
    bool validated_entry   = false;
    ulong validated_offset = 0;
    
    // Max no. of lines to read equals 7.
    // LINE1: fragment of FASTQ header  <--- absence of ^@ means invalid start site
    // LINE2: ^(ATGC)+$                 <--- absence of ^@ means invalid start site
    // LINE3: ^+(.*)$                   <--- absence of ^@ means invalid start site
    // LINE4: <QUALITY ENCODING>        <--- absence of ^@ means invalid start site
    // LINE5: ^@(\S+)$                  <--- entry_found = 1, entry_offset = FASTQ.tell()
    // LINE6: ^(ATGC)+$                 <--- absence of ^@ means invalid start site
    // LINE7: ^+(.*)$                   <--- validates entry start site above
    
    // If there happens to be two ^@ in a row, then the only valid situation 
    // this could happen in is if the first ^@ belongs to the quality encoding line.
    // That means the second ^@ may be the official start.
    // 
    // LINE1: fragment of FASTQ header  <--- absence of ^@ means invalid start site
    // LINE2: ^(ATGC)+$                 <--- absence of ^@ means invalid start site
    // LINE3: ^+(.*)$                   <--- absence of ^@ means invalid start site
    // LINE4: ^@(\S+)$                  <--- entry_found = 1, entry_offset = FASTQ.tell()
    // LINE5: ^@(\S+)$                  <--- entry_found = 1, entry_offset = FASTQ.tell()
    // LINE6: ^(ATGC)+$                 <--- absence of ^@ means invalid start site
    // LINE7: ^+(.*)$                   <--- validated_entry = 1, validated_offset = entry_offset

    // Move file pointer to given offset
    FASTQ.seek(offset, SEEK_SET);

    OUTER: while(
          !validated_offset                         // entry position not found
      &&  FASTQ.readln(buf)                         // EOF not encountered
         ) {

        if(buf[0] == '@') {
            entry_found  = true;
            entry_offset = FASTQ.tell()-buf.length;
            debug writeln("----> entry found at offset ",entry_offset);
            continue OUTER;
        }
        
        if(entry_found && buf[0] == '+') {
            validated_entry = true;
            validated_offset = entry_offset;
            debug writeln("----> entry validated at offset ", validated_offset);
            break OUTER;
        }

    }

    if(validated_entry == true) {
        nextOffset = entry_offset;

        debug writeln("----> seeking to offset ", nextOffset);
        // DEBUG:
        FASTQ.seek(nextOffset, SEEK_SET);
        char[] entry;
        FASTQ.readln(entry);
        std.stdio.write("----> ENTRY HEADER:", entry);
    }


    return nextOffset;
}
/******************************************************************************/
void profileFastqFileOffsetVitals (std.stdio.File FASTQ, ref inputOptions inOpts) {

    // Estimate size of a single entry (4 lines) from the first entry
    char[] buf;
    ulong bufLen;
    
    inOpts.inFastqFileSize = getSize(inOpts.inFastq);

    foreach(i; 0..4) {
        FASTQ.readln(buf);
        bufLen += buf.length;
    }

    inOpts.fastqEntrySize = bufLen;
    debug writeln("A single FASTQ entry takes up ",inOpts.fastqEntrySize," bytes, on average");

    inOpts.fastqEntries           = cast(ulong)inOpts.inFastqFileSize/inOpts.fastqEntrySize;
    inOpts.fastqEntriesPerThread  = cast(ulong)inOpts.fastqEntries/inOpts.threads;
    inOpts.fastqEntriesPerThread |= 1;
    inOpts.fastqExactEntries      = cast(float)inOpts.inFastqFileSize/inOpts.fastqEntrySize;
    inOpts.effThreadCount         = 0;

    debug {
        writeln("FileSize  = ", inOpts.inFastqFileSize," bytes");
        writeln("Buffer    = ",inOpts.fastqEntrySize," bytes");
        writeln("# Entries = ", inOpts.fastqEntries," (approximate)\t",inOpts.fastqExactEntries," (exact)");
        writeln("Threads   = ", inOpts.threads);
        writeln("# Entries/thread = ", inOpts.fastqEntriesPerThread," (approximate)");
    }

    // Given the average position of each potential FASTQ entry offset as a 
    // starting point, we extract the exact position of the nearest FASTQ entry
    // and specify its end offset as well.
    foreach (t; 0..inOpts.threads) {

        /*
        start = (entriesPerThread * t)*bufLen
        end   = (entriesPerThread * (t+1)*bufLen - 1
              = entriesPerThread * (t*bufLen + bufLen) - 1
              = (entriesPerThread * t * bufLen) + (entriesPerThread * bufLen) - 1
              = start + entriesPerThread*bufLen - 1
        */
        ulong seekStart = (inOpts.fastqEntriesPerThread)*(t)*(inOpts.fastqEntrySize);
        ulong seekEnd   = (t == inOpts.threads-1) 
                        ? (inOpts.inFastqFileSize)
                        : seekStart + (inOpts.fastqEntriesPerThread)*(inOpts.fastqEntrySize) - 1;

        ulong nextOffset = getFastqEntryOffset(FASTQ, seekStart);
        if(nextOffset == -1) break;

        debug writeln("[",t,"] MEAN Byte Offset = ", seekStart,"\tEND: ", seekEnd);

        inOpts.partitionOffsets ~= nextOffset;     // Add next valid offset
        ++inOpts.effThreadCount;
        
        debug {
            writeln("----> chunk begins at byte offset ", nextOffset);
            writeln();
        }
        
    }

    // Output effective thread count (since FASTQ file might not be partitionable 
    // into the provided number of threads
    writeln("Threads: ",inOpts.effThreadCount," (effective)\t",inOpts.threads," (requested)");

    // Last entry will never be the total file size, since that means no entries
    // could follow.
    debug { 
        foreach(b, ulong borderStart; inOpts.partitionOffsets) {
            ulong chunkEnd = (b == inOpts.partitionOffsets.length-1)
                           ? (inOpts.inFastqFileSize)
                           : (inOpts.partitionOffsets[b+1]-1);
            FASTQ.seek(borderStart, SEEK_SET);
            char[] hdr;
            FASTQ.readln(hdr);
            writeln(b, ":", borderStart,"-",chunkEnd,":\"",hdr,"\"");
        }
    }

    FASTQ.rewind;

}
/******************************************************************************/
ulong determineFlushLimit(in uint effThreadCount) {

    // Standard: 500k reads per 4 threads
    return cast(ulong) (125_000*effThreadCount);

    //return 25_000;
    //return 100_000;
    //return 200_000;
    //return 1_000_000;
    //return 500_000;
    //return 50_000;
}
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
/******************************************************************************/
void main(string[] args) {
    
    StopWatch swProgram;
    swProgram.start();

    string[] args_copy = args.dup;

    ////////////////
    // CONSTANTS
    ////////////////
    ubyte              asciiEncoding    = 33;
    ushort             userMinL         = 0;    // min should be 30
    ushort             userFixL         = 0;    // min should be 30
    bool               recycle          = false;
    uint               threads          = 1;
    ubyte              userMinQ         = 10;
    ushort             userMaxN         = 50;
    float              lowComplexity    = 0.85;
    float              minAveQ          = 20.0;
    bool               sortLenAsc       = false;
    bool               sortLenDesc      = false;
    bool               showUsage        = false;
    bool               showHelp         = false;
    bool               verbose          = false;
    ulong              trashedBases     = 0;
    string             outPath;
    ubyte              outEncoding      = asciiEncoding;
    byte               outOffset;
    string             prefix;
    string             inFastq;
    string             statsFilename;
    string             histoFilename;
    string             progName         = "splitrim";
    string             progVersion      = "0.1j3e6";
    bool               showVersion      = false;
    

    getopt(
        args, std.getopt.config.caseSensitive,
        "inFile",       &inFastq,
        "prefix",       &prefix,
        "ascii",        &asciiEncoding,     // 33 or 64
        "outPath",      &outPath,           // string
        "minL",         &userMinL,          // > 0, < stringlength
        "fixL",         &userFixL,          // > 0, < stringlength
        "recycle",      &recycle,           // true or false
        "t|threads",    &threads,           // 1
        "minQ",         &userMinQ,          // 0-41
        "maxN",         &userMaxN,          // > 0, < stringlength
        "minAveQ",      &minAveQ,           // 0.0-41.0
        "outEncoding",  &outEncoding,       // 33 or 64
        "sortLenAsc",   &sortLenAsc,        //
        "sortLenDesc",  &sortLenDesc,       //
        "usage",        &showUsage,
        "h|help",       &showHelp,
        "v|verbose",    &verbose,
        "V|version",    &showVersion,       //
        "stats",        &statsFilename,     // General statistics
        "histo",        &histoFilename,     // Post-trim read length histogram
    );

    ////////////////////////////////////////////
    // Retain all user input parameters 
    ////////////////////////////////////////////
    inOpts.progName      = progName; //args[0];
    inOpts.progVersion   = progVersion;
    inOpts.prefix        = prefix;
    inOpts.asciiEncoding = asciiEncoding;
    inOpts.outPath       = outPath;         // o_path
    inOpts.userMinL      = userMinL;        // opt_min_L
    inOpts.userFixL      = userFixL;        // 
    inOpts.recycle       = recycle;         //
    inOpts.threads       = threads;         //
    inOpts.trashedBases  = trashedBases;    // No. of bases trashed from using --fixL without --recycle
    inOpts.userMinQ      = userMinQ;        // opt_q
    inOpts.userMaxN      = userMaxN;        // N_num_cutoff
    inOpts.minAveQ       = minAveQ;         // average read quality (after trim)
    inOpts.outEncoding   = outEncoding;     // out_offset
    inOpts.outOffset     = outOffset;
    inOpts.sortLenAsc    = sortLenAsc;
    inOpts.sortLenDesc   = sortLenDesc;
    inOpts.newIlluminaFastqFormat = false;
    inOpts.verbose       = verbose;
    if(showVersion) {
        writeln(inOpts.progName," ", inOpts.progVersion);
        exit(0);
    }
    if(showUsage || showHelp || (outPath.length == 0) || (inFastq.length == 0)) usage();
    inOpts.inFastq         = inFastq;         // $_[1]
    inOpts.statsFilename   = statsFilename;
    
    inOpts.baseLength = (inOpts.userMinL > 0 )  ? inOpts.userMinL 
                      : (inOpts.userFixL > 0 )  ? inOpts.userFixL
                      : (0);
                      
    ////////////////////////////////////////////
    // Enforce some basic input requirements
    ////////////////////////////////////////////
    processOptions1(inOpts);

    //inOpts.inFastqBasename = baseName(stripExtension(inFastq));
    inOpts.inFastqBasename = inOpts.prefix;
    inOpts.outFastq        = buildPath(inOpts.outPath, inOpts.inFastqBasename~"_splitrim.fastq");

    //////////////////////////
    // Statistics File
    //////////////////////////
    inOpts.outStatsFilename = (inOpts.statsFilename.length == 0)
                            ? buildPath(inOpts.outPath, inOpts.inFastqBasename~"_splitrim.stats.txt")
                            : buildPath(inOpts.outPath, inOpts.statsFilename);

    // if a path + filename was provided, use that path & filename
    // if only a filename was provided, use inOpts.outPath & given filename
    // if nothing is provided, use inOpts.outPath and default filename
    inOpts.outHistoFilename = (histoFilename.length == 0)
                            ? buildPath(inOpts.outPath, inOpts.inFastqBasename~"_splitrim.histo.txt")
                            : buildPath(inOpts.outPath, histoFilename);
                            
    // ===========================
    // FASTRIM EDITED (2012-11-09):
    // ===========================
    //std.stdio.File IN, OUT, LENGTH, STATS;
    std.stdio.File IN, OUT, STATS, HISTO;

    ////////////////////////////////////////////
    // Enforce addtnl basic input requirements
    ////////////////////////////////////////////
    processOptions2(inOpts);

    ////////////////////////////////////////////
    // Get filehandles for all file I/O
    ////////////////////////////////////////////
    initFiles(inOpts, IN, OUT, STATS, HISTO);                // Assign filehandles (read, write, write)

    ////////////////////////////////////////////
    // Peek at FASTQ file: is it valid?
    ////////////////////////////////////////////
    verifyFastqFormat(IN, inOpts.newIlluminaFastqFormat);

    //ulong[ushort]       trimSeqLenHash;
    //ulong               rawSeqNum, trimSeqNum, totalRawSeqLen, totalTrimSeqLen = 0;

    ////////////////////////////////////////////
    // Profile FASTQ file for offsets
    ////////////////////////////////////////////
    profileFastqFileOffsetVitals(IN, inOpts);
    inOpts.flush         = determineFlushLimit(inOpts.effThreadCount);

    ////////////////////////////////////////////
    // Determine staggered flush output value
    ////////////////////////////////////////////
    ulong   stagger = cast(uint) (inOpts.flush/inOpts.effThreadCount);
    ulong[] staggeredFlushCount;
    for(int i = 0; i < inOpts.effThreadCount; i++) {
        staggeredFlushCount ~= (i+1)*stagger;
    }

    StopWatch swTrim;
    swTrim.start();
    defaultPoolThreads(inOpts.effThreadCount);           // sets no. of threads avail
    auto workers = new TaskPool;
    foreach(i, ref offsetStart; workers.parallel(inOpts.partitionOffsets)) {
        parseFASTQ(i, offsetStart, inOpts, staggeredFlushCount[i], OUT); //inOpts.inFastq, inOpts.inFastqFileSize);
    }    
    workers.finish();
    swTrim.stop();
    writeln("GLOBAL READ COUNT = ", inOpts.rawSeqNum);
    writeln("Trim Time: ",swTrim.peek().msecs, " ms");

    /////////////////////////////////
    // Calc read statistics
    /////////////////////////////////
    seqLenStats myStats;
    
    // EDIT: NEED TO UPDATE FOR trimSeqLenHash COMPATABILITY
    bool validSeqs = calcStats(myStats, inOpts.trimSeqLenHash);    // returns true if calc is OK
    if(!validSeqs) {
        writeln("\n**Warning**: fault detected in statistics calculation!");
        writeln("    -> Check either (1) the input ASCII encoding or (2) the # of reads passing filters");
    }

    /////////////////////////////////
    // Close up shop
    /////////////////////////////////
    IN.close();
    OUT.close();

    writeln;
    // edit
    auto programTime = swProgram.peek().msecs;
    writeln("PROGRAM ELAPSED TIME: ", swProgram.peek().msecs, " ms");
    writeStatsFile(STATS, inOpts, myStats, programTime, args_copy);
    writeHistoFile(HISTO, inOpts);

    swProgram.stop();
    
    return;
}
/***************************************************************************************/
void usage() {

    writeln;
    writeln(inOpts.progName, " ", inOpts.progVersion);
    writeln;
    writeln("Usage: ", inOpts.progName, " [INPUT] [OUTPUT] [OTHERS]");
    write(
        "\n",       
        "Input:\n",
        "==============\n",
        "REQUIRED:\n",
        "   --inFile=        string   Name of the FASTQ file containing all the single-end reads\n",
        "\n",
        "   --minL=          int      Minimum length for a trimmed read to be considered valid [default: ",inOpts.userMinL,"]\n",
        "or --fixL=          int      Fixed length to which each trimmed read will be cut down to [default: disabled]\n",
        "OPTIONAL:\n",
        "   --recycle        bool     When --fixL is specified and a read length is not a multiple of \"fixL\", this option will append any\n",
        "                             remaining bases (up to the maximum fixL-1 bases) to the last fragment of length \"fixL\" [default: false]\n",
        "   --ascii=         int      ASCII encoding (33 or 64) [default: ",inOpts.asciiEncoding,"]\n",
        "   --minQ=          int      Minimum quality for a read to be considered valid (0-",qBest,") [default: ",inOpts.userMinQ,"]\n",
//        "   --maxN=          int      Maximum number of N-bases for a read to be considered valid [default: ",inOpts.userMaxN,"]\n",
//        "   --minAveQ=       float    Minimum average read quality (post trim) (0.0- [=DISABLED=]");
//    writef("%.1f",cast(float)qBest);
//    writeln(") [default: ",inOpts.minAveQ,"]\n",
        "   -t, --threads=   uint     <disabled> no. of threads to use [1]\n",
        "\n",
        "Output:\n",
        "==============\n",
        "REQUIRED:\n",
        "   --outPath=       string   Location output files will be placed\n",
        "OPTIONAL:\n",
        "   --prefix=        string   Prefix of output files\n",
        "   --outEncoding=   int      ASCII encoding of the output (33 or 64) [default: mirrors input]\n",
        "   --stats=         string   Basic read statistics output [default: uses basename from --inFile]\n",
        "   --histo=         string   Post-trim read length histogram [default: uses basename from --inFile]\n",
        "   --sortLenAsc     bool     Sort read length frequency table in ascending order [default: unordered]\n",
        "   --sortLenDesc    bool     Sort read length frequency table in descending order [default: unordered]\n",
        "\n",
        "Others:\n",
        "==============\n",
        "   -h, --help       bool     display HELP\n",
        "   --usage          bool     display this message\n",
        "   -v, --verbose    bool     verbosity level\n",
        "   -V, --version    bool     print program version and exit\n",
        "\n\n",
        "Ex1: Split single-end reads at bases < Q20 and retain all read fragments >= 30-bp, transcoding qualities from ASCII-33 to ASCII-64\n",
        "     ",inOpts.progName," --inFile=/user/518893/seq.fastq --ascii=33 --minL=30 --minQ=20 --prefix=MYSEQ --outPath=/user/518893/output --outEncoding=64 --sortLenAsc\n",
        "\n",
        "Ex2: Split single-end reads at bases < Q20, split remaining frags into 30-bp frags, trashing the last fragment if < 30-bp, and transcode qualities from ASCII-33 to ASCII-64\n",
        "     ",inOpts.progName," --inFile=/user/518893/seq.fastq --ascii=33 --fixL=30 --minQ=20 --prefix=MYSEQ --outPath=/user/518893/output --outEncoding=64 --sortLenAsc\n",
        "\n",
        "Ex3: Split single-end reads at bases < Q20, split remaining frags into 30-bp frags, appending the last fragment (if < 30-bp) to previous fragment, and transcode qualities from ASCII-33 to ASCII-64\n",
        "     ",inOpts.progName," --inFile=/user/518893/seq.fastq --ascii=33 --fixL=30 --recycle --minQ=20 --prefix=MYSEQ --outPath=/user/518893/output --outEncoding=64 --sortLenAsc\n",
        "\n",
    );
    writeln;
    exit(0);
    
}
/******************************************************************************/

    // Case 1: Catch a single 5'-low-quality base
    // Case 2: Catch contiguous 5'-low-quality bases
    // Case 3: Catch a single 3'-low-quality base
    // Case 4: Catch contiguous 3'-low-quality bases
    // Case 5: Catch both 5'- and 3'-low quality bases
    // Case 6: Catch one string of contiguous low-quality bases in the middle
    // Case 7: Catch multiple strings of contiguous low-quality bases in the middle 
    // ----------------------------------------------------------------------------
    // Case 8: Catch a single, high quality, 5' N
    // Case 9: Catch contiguous, high quality, 5' N's
    // Case 10: Catch a single, high quality, 3' N
    // Case 11: Catch contiguous, high quality, 3' N's
    // Case 12: Catch a single string of contiguous, high-qualityt N's in the middle
    // Case 13: Catch multiple strings of contiguous, high-qualilty N's in the middle
    // ----------------------------------------------------------------------------
    // Case 14: Catch contiguous, low-quality bases at the ends (5' and 3'), along
    //          with contiguous, high-quality strings of N's in the middle
    // Case 15: Catch contiguous, high-quality N's at the ends (5' and 3'), along
    //          with contiguous, low-quality bases in the middle. 
    
/*

GLOBAL READ COUNT = 8966559
PROGRAM ELAPSED TIME: 173512 ms

        t   flush       ms          TrimTime(s) read count  SystemRAM(GB)   Cores
      ---   ------      ---------   ---------   -------------   -------
        2    25000      173512                  8966559     16              4
        2   100000      179741                  8966559     16              4
        2  1000000      255656                  8966559     16              4

        2    25000      170471                  8966559     16              4   (used && with readCounter)
        4    25000      182330                  8966559     16              4   (used && with readCounter)        
        
        1   100000      162215                  8966559     16              4   (used && with readCounter)
        2   100000      173448
        4   100000      193537
        8   100000      290494
        
        4   100000      204606                                                  (Lock-Free; hitting same file)
        4   100000      206233      187                                         (L-F)
        4   100000      199519      181                                         (L-F; FASTQ.ftell())
        
        **larger thread count with lower flush values may increase thread wait times!

*/



