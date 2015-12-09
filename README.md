# TEorthology
This script was written to help assessing the evolutionary history of a TE family (=> presence or absence of the copies based on homology). 

It runs on any number of TE annotations (repeat masker outputs) and can interogate any number of genomes (it is threaded).

Briefly, all insertions not nested in another (masked) repeat are extracted with 5’ and 3’ flanking sequences (100 nt by default), 
and queried against all genomes of interest using BLASTN with a default E cutoff of < 1x10-50.
Hits are then filtered based on the presence of both 5’ and 3’ flanking sequence, which is used as a proxy to determine if matches in 
other genomes represented syntenic/orthologous elements. Both the total number and number of potentially orthologous matches are recorded.


    USAGE (v3.6)
    This script will assess lineage for various repeats, from Repeat Masker output(s) - and output a summary table.
    Orthology will be approximated by presence of 50% of the flankings of the element in the best hit.	
    Note that outputs will be located in the same DIRECTORY of the <dir_with_all_files> (so you need writing access)
    
    MANDATORY ARGUMENT:	
    -dir (STRING)
        => directory that should contain:
           - at least 1 repeat masker output [there can be several]
           - Corresponding genomes (exact same ones that are masked) 
               /!\\ Same name is required before extensions = .out and .fa; note that .out.XXXX is OK if you need more info in the file name
               NOTE: sym links will work (ln -s path/file new_name) so no need to actually copy the genome files
        For ex. minimal list of files in this directory = species1.out, species1.fa
        Outputs will be located in this folder too.
	  
    OPTIONAL ARGUMENTS (flagged by brackets [...] around each of them)
     -tree (TSTRING)
        This can be a file .dng OR directly the relationships on command line using \"\"
        To order species in the output, based on phylogenetic relationships. Follow Newick format = (A,B,(C,D));
        ex: -tree \"(((hg19,panTro4),rheMac3),(canFam3,felCat5));\"
            whill lead to order in output being: hg19,panTro4,rheMac3,canFam3,felCat5 + names with parenthesis will be printed on line above
            /!\\ IMPORTANT: 
               - names have to exactly match genome file name, anything before .fa
               - parenthesis and ; at the end are required
    -presplit (BOOL)
    	if RMoutput files have been split already => do not do it again
    	Note that filtering from -filter is done at that step; unless it was -filter name,XX do not use -presplit if -filter changed
    -blast (STRING)
        path to blast bin. 
        Default = /home/software/ncbi-blast-2.2.28+/bin
    -e <evalue>
        set minimum evalue for blast search. 
        Default = 10e-50; same species blast, evalue is set to 10e-120
        Only the best hit is printed in the blast output, but lowering evalue will increase the speed.
        However, for distant species such as 60 My, this should not be decreased to much.
        The best way to set this evalue is to manually check what evalue you would expect for a blast of non coding regions between species of interest.
    -filter <type,name>
        run the script on only a subset of repeats. Not case sensitive.
        This is advised instead of grepping the Repeat Masker output, because elements will be skipped when nested (e.g. flankings are repeated as well).
        The type can be: name, class or family and it will be EXACT MATCH unless -contain is chosen as well
        ex: name,nhAT1_ML => only fragments corresponding to the repeat named exactly nhAT1_ML will be looked at
            class,DNA => all repeats with class named exactly DNA (as in ...#DNA/hAT or ...#DNA/Tc1)
            family,hAT => all repeats with family named exactly hAT (so NOT ...#DNA/hAT-Charlie for example)
        Importantly, unless you use -filter name,XX -presplit SHOULD NOT BE USED (can't filter on class or family)
    -contain (BOOL)
        to check if the \"name\" determined with -filter is included in the value in Repeat Masker output, instead of exact match
        ex: name,HERVK => all fragments containing HERVK in their name
            family,hAT => all repeats with family containing hAT (...#DNA/hAT, ...#DNA/hAT-Charlie, etc)
        No effect if -list is used and not -filter
    -list <file>
         to filter on a list of TEs
         file = path of a file with a list of repeat names, one per line (if full names are given, only what is before # or any space will be considered)
    -min_frg (INT)
        filter for a minimum number of fragments (potentially useful if -filter class is chosen and library is not really clean)
        note that this filtering is done before reconstruction of interrupted repeats (X is really the number of fragments and not of elements)
        ex: -min_frg 3   
    -min_len (INT)
        filter for a minimum length of the fragment (in nt) to avoid considering very small pieces in orthology assessment
        ex: -min_len 80
    -flank (INT)
        To change the length of flankings extracted with the copy (fragment) of the element. 
        At least some flanking is necessary, to be able to probe for potential orthology.
        Default = 100.
    -TEs (STRING)
        optional input file, to allow correction of class and/or family of repeats before filtering
        with TE information as follow, 3 first columns mandatory: 
        Rname \\t Rclass \\t Rfam \\t Rclass/Rfam \\t etc (only first 3 columns will be used)
        an easy way to get it is to run my other script parseRM.pl.
        -> copy the first columns of the output all-repeats.tab, modify them accordingly and then copy this in a text file
    -out (STRING)
        To rename the <out> part of the output directory [Default = orthology]
        Output = <where_project_is>/<out>.<project>
                 Whith project = name of the directory set with -dir
    -append (BOOL)
    	existing output files won't be deleted; newly extracted sequences will be added
    	and blasts will be rerun only for sequences with no blast output
    	(chose this if a genome or a Repeat masker output was added in the directory)
    	(also useful if the job did not go through; somehow threads die with no errors sometimes)
    	Note that use of -append will automatically activate -presplit
    -ponly (BOOL)
    	Nothing will be run besides parsing existing blast outputs
    	You should still provide all options from the original run
    -cpus (INT)
        chose the amount of CPUs used to run the script.
        Default = 1
        More than 1 is strongly advised though, it can be very long.
    -v (BOOL)
        verbose mode, make the script talks to you
        print version if only option
    -chlog (BOOL)
    	print change log (updates)
    -h|help (BOOL) 
        Print this help
