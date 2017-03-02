#!/usr/bin/perl -w
##########################################################
# Author  :  Aurelie Kapusta
# version :  see below + see changelog
# email   :  4urelie.k@gmail.com
# PURPOSE :  Evaluate if a TE is shared between species of interest. 
#			 Check Usage for the many options. Assemblies as well as repeat masker outputs are required.
##########################################################
BEGIN{
   #what to do on: kill -s ALRM <pid> so I can check where it is if it stalls
   $SIG{ALRM}  = sub {print STDERR "SIGALRM received\n"; print STDERR Carp::longmess; print "\n";};
   #what to do on ^C
   $SIG{INT}  = sub {print STDERR "SIGINT received\n"; print STDERR "\n\n".Carp::longmess; exit;};
   #add a folder in INC
   #unshift(@INC, "~/bin/BioPerl-1.6.901");
}
#always load forks before anything else
use forks;
use forks::shared;
#load the rest
use warnings;
use strict;
use Carp;
use Getopt::Long;
use Bio::DB::Fasta;
use Bio::SearchIO; 
use Bio::SeqIO;

my $version = "3.7";

my $changelog = "
#	- v1.0 = May 2014
#	- v1.1 = 23 Jun 2014
#		use blast option to get only best hit and avoid huge outputs + gain time (-num_alignment 1)
#		bug fixes
#		filter out non TE stuff
#	- v2.0 = 03 Jul 2014
#		PARALLELIZATION [ required (lots of) changes ]
#			=> mix of script parseRM_ExtractSeqs_P.pl + TE_Orthology.pl
#		filtering using -filter is now done after filtering out nested elements
#	- v2.1 = 08 Jul 2014
#		Corrected cleaning previous outputs
#		Return in thread loop that shouldn't be there -> changed as a next FILE
#	- v2.2 = 29 Jul 2014
#		extract_one_seq was modified to avoid regexp against whole db header list, it was a long useless step + some variables were not use, remains from a cc of other script
#		missing \"return\" at the end of blast sub
#		Bug fixed: var se remained in extract_sequences (semaphore variable)
#		mv empty output folders instead of rm so I can see which ones since STDERR doesn't print everything for some reason
#		\$v in thread [still trying to find solutions for stuff to print]
#	- v2.3 = 06 Aug 2014
#		script works and run but somehow only one thread was running... Add defined in the loop to get tasks + print before the return of the loop just in case
#		mv empty output folders was done to wrong folder
#		added a DATA folder to avoid having too many folders in the output folder (easier to find final output)
#	- v2.4 = 12 Aug 2014
#		Several bugs were fixed 
#			- filter on family was not working
#			- problem headers in extract_one_seq from the cc of other script
#	- v2.5 = 13 Aug 2014
#		Threads were returning even though there were still files to process => changed the loop in thread subroutine
#	- v3.0 = 15 Jan 2015
#		Blasts were too long. Solution = use the full power of blast and output in tabular + already ask for only best hits.
#       Also added a step to get and store lengths of sequences of genomes (the Bio:DB length tool is not working well for some reason).
#       It takes longer but once it's done info can be extracted back
#	- v3.1 = 20 Jan 2015
#		Still too long; lowered evalue for same species blasts to 10e-120
#	- v3.2 = 04 Feb 2015
#		Added -contain for filtering + actually fixed the filtering
#	- v3.3 = 24-25 Mar 2015
#		-append now checks all files => no overwrite of anything
#       Bug fix in -presplit
#	- v3.4 = 26/30 Mar 2015
#		-out added
#       Added -list <X> option
#       Changes in -presplit behavior, so that -presplit can be used with -filter name,XX or -list
#          (so that several jobs with different filters can be run at the same time, as long as all repeats where split)
#          (Still not possible for -filter class or family)
#	- v3.5 = 30 Apr 2015
#		Tiny bugs fixed
#	- v3.6 = 14 Sep 2015
#		-list option was not set up properly (not included in the filtering while splitting RMout files)
#	- v3.7 = 23 Feb 2017
#       parse blast with ls => error from the use of * by a user => use find command instead

# TO DO
# Better and faster if -append would detect directly the .out that have been processed or not...
\n";

my $usage = "\nUsage [v$version]: 
	perl <scriptname.pl> -dir <dir_with_all_files> [-tree <tree>] [-presplit] [-blast <path_to_ncbi-blast>] [-e <evalue>] [-filter <type,name>] [-contain] [-min_frg <X>] [-min_len <X>] [-flank <X>] [-TEs <TEclass>] [-out] [-append <X>] [-ponly] [-cpus <X>] [-v]
	
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
        Print this help\n\n";      
        
#keep STDOUT and STDERR from buffering [could be issue I have with nohup]
select((select(STDERR), $|=1)[0]); #make STDERR buffer flush immediately
select((select(STDOUT), $|=1)[0]); #make STDOUT buffer flush immediately

################################################################################
# Get arguments/options, check some of them, print details of the run if verbose chosen
################################################################################
my $blastloc = "/home/software/ncbi-blast-2.2.28+/bin";
my $evalue = "10e-50";
my $flank = "100";
my $cpus = 1;
my $outname = "orthology";
my ($tree,$filter,$flist,$min_frg,$min_len,$TEclass) = ("na","na","na","na","na","na");
my ($dir,$presplit,$f_regexp,$append,$tinyponey,$help,$v,$chlog);
GetOptions ('dir=s' => \$dir, 'tree=s' => \$tree, 'blast=s' => \$blastloc, 'e=s' => \$evalue, 'filter=s' => \$filter, 'list=s' => \$flist, 'contain' => \$f_regexp, 'min_frg=s' => \$min_frg, 'min_len=s' => \$min_len, 'flank=s' => \$flank, 'TEs=s' => \$TEclass, 'out=s' => \$outname, 'append' => \$append, 'presplit' => \$presplit, 'ponly' => \$tinyponey, 'cpus=s' => \$cpus, 'chlog' => \$chlog, 'h' => \$help, 'help' => \$help, 'v' => \$v);

($f_regexp)?($f_regexp = "y"):($f_regexp="n");

#
# CHECKING STEPS
#
#check step to see if mandatory argument is provided + if help/changelog
die "\n version $version\n\n" if ((! $dir) && (! $help) && (! $chlog) && ($v));
die $changelog if ($chlog);
die $usage if ((! $dir) || ($help));

print STDERR "\n --- Script TE_Orthology started (v$version)\n" if ($v);

#avoid / at the end of paths + check blast location
$dir = $1 if ($dir =~ /^(.*)\/$/);
$blastloc = $1 if ($blastloc =~ /^(.*)\/$/);
print STDERR "      - Directory containing input files = $dir\n" if ($v);
print STDERR "      - Blast software location = $blastloc\n" if ($v);
die "\n      ERROR (main): please provide accurate blast bin directory location\n\n" unless (-e "$blastloc/blastn");
print STDERR "         -> max. evalue for blasts = $evalue (lowered to 10e-120 for same species blasts)\n" if ($v);

#check $tree value
if ($tree ne "na") {
	if ($tree =~ /\.dng$/) { #if file is provided, load it
		open (my $fh, "<", $tree) or die "\n      ERROR (main): could not open $tree $!\n\n";
		$tree = <$fh>;
		chomp $tree;
	}
	if ($tree !~ /^\(.*;$/) { #now test format of it, e.g. if it starts with ( and ends with ; and die if not
		die "\n      ERROR (main): please provide valid -tree option (not a .dng file and/or not valid format)\n$usage";
	}
	print STDERR "      - Tree info = $tree\n" if ($v);
}
if ($v) {
	print STDERR "      - TE class file provided = $TEclass\n" unless ($TEclass eq "na");
	print STDERR "      - Max number of CPUs used = $cpus\n";
	print STDERR "      - Repeats with same repeat masker block ID will be reconstructed when extracted\n";
	print STDERR "      - Extraction of sequences will be done:\n";
	print STDERR "         with additional $flank nt in 5' and 3' when possible\n";
	print STDERR "            (number of nt added will be added in sequences descriptions)\n";
	print STDERR "         only when repeat has more than $min_frg fragments\n" unless ($min_frg eq "na");
	print STDERR "         only when reconstructed length is > $min_len nt\n" unless ($min_len eq "na");
	print STDERR "         if repeat name matches names provided in $flist\n" unless ($flist eq "na");
	print STDERR "         if not filtered out based on -filter option = $filter\n" unless ($filter eq "na");
	print STDERR "            -> using exact match\n" if (($filter ne "na") && ($f_regexp eq "n"));
	print STDERR "            -> using regular expression\n" if (($filter ne "na") && ($f_regexp eq "y"));
}
#if relevant, check if filter has the ,
die "\n      ERROR (main): check -filter option (use -h if you need to see the usage)\n" if (($filter ne "na") && ($filter !~ /,/));
#Also warn user in case presplit is chosen
my ($ft,$fn) = split(",",$filter);
print STDERR "            -> Note that since -presplit is used, .out files that are going to be processed are from the previous run\n" if (($v) && ($presplit) && ($flist eq "na") && ($ft ne "name"));
#Presplit log
$presplit = 1 if ($append);
print STDERR "      - option -presplit activated since -append chosen\n" if ($append);
if ($v) {
	print STDERR "      - Split RM outputs will be fetched from $dir/*.split\n" if ($presplit);
	($append)?(print STDERR "      - option -append chosen, no cleaning previous outputs\n"):
			  (print STDERR "      - option -append not chosen, previous outputs with same name will be cleaned\n");
}


##########################################################################################################
# "PREP" steps
##########################################################################################################
print STDERR "\n --- Now running prep steps:\n" if ($v);

#Get list of species and ordered
my @Genomes = `ls $dir/*.fa` or die "\n      ERROR (main, prep): can't list files .fa in $dir $!\n\n";
print STDERR "      - getting species list\n" if ($v);
my $sp_ordered = get_species_list(\@Genomes,$tree,$v);
die "\n      ERROR (main, prep): number of species provided in -tree =/= number of genome files provided\n -> check input file or command line\n$usage\n" unless ($sp_ordered);
print STDERR "        ...done\n" if ($v);


#Output directory, in same location as $dir
my $loc = get_path($dir);
my $project = filename($dir);
my $outpath;
($project eq ".")?($outpath = $loc."/".$outname):($outpath = $loc."/".$outname.".".$project);

#Now check if only parsing should be done; some info will have to be fetched though
if ($tinyponey) {
	just_parse($outpath,$flank,$dir,$sp_ordered,$tree,$v);
	print STDERR "\n --- Script TE_Orthology done\n\n" if ($v);
	exit;
}

#clean previous files unless -append or -presplit + give them values
unless ($append) {
	print STDERR "      - Clean previous output directories\n" if ($v);
	if (-e $outpath) {
		print STDERR "        rm -Rf $outpath\n" if ($v);
		`rm -Rf $outpath`; #or die "\n      ERROR (main, prep): can not rm -Rf $outpath $!\n\n";
	}
	my @checksplit = `ls $dir | grep split`;
	my $checksplitnb = @checksplit;
	if ($checksplitnb > 0) {
		my @checktab = `ls $dir/*.split | grep tab`;
		my $checktabtnb = @checktab;
		if ($presplit) {
			if ($checktabtnb > 0) {
				print STDERR "        rm -f $dir/*.split/*tab\n" if ($v);
				`rm -f $dir/*.split/*tab`; #or die "\n      ERROR (main, prep): can not rm -f tab files in $dir/*.split directories $!\n\n";
			}	
		} else {
			print STDERR "        rm -Rf $dir/*.split\n" if ($v);
			`rm -Rf $dir/*.split`; #or die "\n      ERROR (main, prep): can not rm -Rf $dir/*.split $!\n\n";
		}	
	}
	$append = "no";
} else {
	$append = "yes";
}
#Creating output directory
print STDERR "      - Creating output directory\n" if ($v);
if (-e $outpath) {
	mkdir ($outpath."/DATA_Fa_no-seq", 0755) unless (-e $outpath."/DATA_Fa_no-seq");
	mkdir ($outpath."/DATA", 0755) unless (-e $outpath."/DATA");
	print STDERR "        -> $outpath existed, not created (-append chosen)\n" if ($v);
} else {
	mkdir ($outpath, 0755) or die "\n      ERROR (main, prep): can not mkdir $outpath $!\n\n";
	mkdir ($outpath."/DATA_Fa_no-seq", 0755) or die "\n      ERROR (main, prep): can not mkdir $outpath/DATA_Fa_no-seq $!\n\n";
	mkdir ($outpath."/DATA", 0755) or die "\n      ERROR (main, prep): can not mkdir $outpath/DATA $!\n\n";
	print STDERR "        -> $outpath + subdirectories created\n" if ($v);
}

#Get TE infos if provided
my $TE->{"na"} = "na"; #will be replaced with TE info if TEclass provided
unless ($TEclass eq "na") {
	print STDERR "      - getting TE infos from $TEclass\n" if ($v);
	$TE = get_TEs_infos($TEclass);
	print STDERR "        ...done\n" if ($v);
}	

#get list of genomes + do the makeblastdb files now + get lengths of sequences + make index files
print STDERR "      - Doing operations on genome files required for the rest of the script\n        (make Bio::DB index files + creating blast db = makeblastdb + getting lengths of sequences)\n" if ($v);
my $genlen = ();
foreach my $g (@Genomes) {
	chomp ($g);
	# index the genomes if necessary
	my $index = "$g.index";
	if (-e $index) {
		print STDERR "         -> $index file exists, skipped\n" if ($v);
	} else {
		my $db = Bio::DB::Fasta->new($g, -reindex => 1) or confess "\n      ERROR (Main): Failed to open create Bio::DB object (to index) from $g $!\n\n";
		print STDERR "         -> $index file created\n" if ($v);
		undef ($db);
	}	
	#Now deal with makeblastdb
	if (! -e "$g.nhr") {
		system "$blastloc/makeblastdb -in $g -dbtype nucl -out $g -logfile $g.makeblastdb.log";
		print STDERR "         -> makeblastdb $g done\n" if ($v);
	} else {
		print STDERR "         -> $g makeblastdb files exist, skipped\n" if ($v);
	}
	#Now get lengths
	$genlen = get_lengths_noBio($g,$genlen,$v);
}		

#list of repeats from -list if relevant
my @flist = ();
if ($flist ne "na") {
	open (my $lfh, "<", $flist) or confess "\n      ERROR (main, prep): can't open to read $flist $!\n\n";
	LIST: while(<$lfh>) {
		chomp(my $r = $_);
		next LIST if ($r !~ /\w/);
		$r =~ s/^>//;
		$r = $1 if ($r =~ /^(.+?)#/);
		$r = $1 if ($r =~ /^(.+?)\s/);
		push(@flist,$r);
	}
	close $lfh;
} else {
	$flist[0] = "na";
}

##########################################################################################################
# MAIN
##########################################################################################################
#Initialize and open Thread stuff
my @RMout_list :shared; 
my $RMout_list_nb :shared;
my @RMout_done :shared;
my @posi_list :shared; 
my @posi_done :shared; 
my %frgs_all :shared;
my %frgs_nr :shared;

#Split RMout files per repeat and get list of RMout files
if ($presplit) {
	print STDERR "     - No splitting RM output files (-presplit chosen) => gathering previous split RMout in $dir/*split\n" if ($v);
	my @RMout_files = `find $dir/*split -type f -name "*.out*##*" | grep -v ".tab"` or die "\n      ERROR (main): can't find files *.out*##* in $dir/*split $!\n\n";
	my $list = check_for_genome(\@RMout_files,$dir,$v); #check that corresponding genome exists before doing anything
	@RMout_list = @{$list};	
	$presplit = 1;
} else {
	print STDERR "     - Splitting RM output files and getting list of RM output files to process...\n" if ($v);
	my @RMout_files = `ls $dir/*.out*` or die "\n      ERROR (main): can't list files .out in $dir $!\n\n"; #no need to avoid .tab here since files would have been deleted
	my $split_RMouts = split_RM_files($dir,\@RMout_files,$filter,\@flist,$TE,$TEclass,$f_regexp,$v); #checking that corresponding genome exists is done in sub (gain of time, less files)
	@RMout_list = @{$split_RMouts};
	print STDERR "     ...Splitting done\n" if ($v);
	$presplit = 0;
}

print STDERR "\n --- Main steps now\n" if ($v);
$RMout_list_nb = @RMout_list;
print STDERR "     => Starting $cpus thread(s) to process a total of $RMout_list_nb RM output files...\n" if ($v);

#start threads
for(my $i = 1; $i < $cpus; $i++){
    threads->create({ 'context' => 'scalar' }, \&bring_it_on, \@RMout_list, \$RMout_list_nb, \@RMout_done, \%frgs_all, \%frgs_nr, \@posi_list, \@posi_done, \$flank, \$min_frg, \$min_len, \$filter, \$f_regexp, \@flist, \$dir, \$outpath, $genlen, \$append, \$presplit, \$blastloc, \$v);
}

#run threads
bring_it_on(\@RMout_list, \$RMout_list_nb, \@RMout_done, \%frgs_all, \%frgs_nr, \@posi_list, \@posi_done, \$flank, \$min_frg, \$min_len, \$filter, \$f_regexp, \@flist, \$dir, \$outpath, $genlen, \$append, \$presplit, \$blastloc, \$v);

#clean threads
print STDERR "\n --- Cleaning the $cpus threads\n" if ($v);
my $pflag = 0;
my $eflag = 0;
foreach my $thr (threads->list){
    my ($flags) = $thr->join();
    $pflag+=$flags->[0];
    $eflag+=$flags->[1];
}
my $totRMdone = 0;
foreach my $done (@RMout_done) {
	$totRMdone++;
}
my $totP = 0;
foreach my $p (@posi_list) {
	$totP++;
}
my $totPdone = 0;
foreach my $p (@posi_done) {
	$totPdone++;
}
#Just check that files were processed and same number of position files generated
print STDERR "     => $totRMdone RMout files processed\n" if ($v);
print STDERR "        $totP position files generated\n" if ($v);
print STDERR "     => total of $totPdone position files processed to extract sequences\n" if ($v);
print STDERR "     !! Some files were not processed ($pflag)\n" unless ($pflag == 0);
print STDERR "     !! Some sequences had no sequences extracted ($eflag)\n        (check $outpath/DATA_Fa_no-seq)\n" unless ($eflag == 0);

#Now parse blast
print STDERR "\n --- Parse blast outputs (check for presence of flankings in the best hit)\n" if ($v);
my ($ortho,$tothit) = parse_blast($outpath,$flank,$v);

#Now print
print STDERR "\n --- Printing results (orthology)\n" if ($v);
get_ortho_output($dir,$outpath,$sp_ordered,$tree,$ortho,$tothit,\%frgs_all,\%frgs_nr,$v);

print STDERR "\n --- Script TE_Orthology done\n" if ($v);
print STDERR "     --> see files in $outpath\n\n" if ($v);
exit;

##########################################################################################################
# SUBROUTINES
##########################################################################################################
#----------------------------------------------------------------------------
# MAIN SUBROUTINE:
# Threaded actions, loop on .out files split by repeat name
# bring_it_on(\@RMout_list, \$RMout_list_nb, \@RMout_done, \%frgs_all, \%frgs_nr, \@posi_list, \@posi_done, \$flank, \$min_frg, \$min_len, \$filter, \$f_regexp, \@flist, \$dir, \$outpath, $genlen, \$append, \$presplit, \$v);
#----------------------------------------------------------------------------
sub bring_it_on {
    my ($RMout_list,$RMout_list_nb,$RMout_done,$frgs_all,$frgs_nr,$posi_list,$posi_done,$flank,$min_frg,$min_len,$filter,$f_regexp,$flist,$dir,$outpath,$genlen,$append,$presplit,$blastloc,$v) = @_; 
	
	my @flags = (0,0);
#	RMFILE: while(defined(my $RMout = shift @$RMout_list) || (! $$finished)){
    RMFILE: while($$RMout_list_nb > 0){
    	my $RMout = shift @$RMout_list;
# 		if(!$RMout){
# 			sleep 1;
# 			next;
# 		}		
		next RMFILE unless ($RMout);
		chomp ($RMout);		
				
		#Check if this file should be parsed based on its name [which is why for now class and family filtering not possible with presplit]
		#Only relevant if -presplit was chosen since otherwise filtering was done there
		my $return = "yes";
		$return = check_Rname($RMout,$filter,$flist,$f_regexp) if ($$presplit == 1);
		if ($return eq "no") {
			$$RMout_list_nb--;
			next RMFILE;
		} else {
			print STDERR "     STARTING: $RMout (thr ".threads->tid().") [$$RMout_list_nb files to process]...\n" if ($$v);
			$$RMout_list_nb--;
		} 			
			
		#OK went through => rest of the pipeline
		#get positions if needed (if -append then existence of file will be checked and frg number read from the file instead
		my $posi;
		($posi_list,$posi) = get_posi_from_RMout($RMout,$frgs_all,$frgs_nr,$posi_list,$append,$v);
		push(@{$RMout_done},$RMout);
		#check posifiles for frg numbers
		my $posiname = filename($posi);
		if (($$min_frg ne "na") && ($frgs_nr->{$posiname} < $$min_frg)) {			
			$flags[0] += 1;
			next RMFILE;
		}		
		#now extract	
		my $RMname = filename($RMout);
		my $faloc = $$outpath."/DATA/".$RMname.".extract";
		`mkdir $faloc` unless (-e $faloc);
		extract_sequences($dir,$faloc,$posi,$flank,$min_len,$frgs_all,$genlen,$append,$v);
		push(@{$posi_done},$posi);		
		#check that fasta file is not empty; if it is, move folder [as a "log"] and go to the next RMoutput
		my @fafiles = `ls $faloc`;
		@fafiles = grep(/\.fa$/, @fafiles);
		unless (@fafiles) {		
			print STDERR "     ...$RMout skipped ($faloc moved), didn't provide any sequence to extract (e.g. min_len or min_frg)\n" if ($$v);
			`mv $faloc $$outpath/DATA_Fa_no-seq/$RMname.extract`;
			$flags[1] += 1;
			next RMFILE;
		}
	
		#blast now
		blast($dir,$faloc,$blastloc,$evalue,$append,$v);	
		
		#done for this RMoutput file
		print STDERR "     ...DONE: $RMout (thr ".threads->tid().") [$$RMout_list_nb files left]\n\n" if ($$v);		
		
		undef(@fafiles);
	}
	my $local_RMout_list_n = @$RMout_list;
	print STDERR "\n     ====> thread ".threads->tid()." returning\n"  if ($$v);
    print STDERR "           => size of list of files still to process = $local_RMout_list_n [should be 0]\n\n" if ($$v);
	return (\@flags);
}


#----------------------------------------------------------------------------
# get a filename from a full path
# my $name = filename($filename);
#----------------------------------------------------------------------------
sub filename {
	my($file) = shift;
	$file =~ s/.*\/(.*)$/$1/;
	return $file;
}

#----------------------------------------------------------------------------
# from a filename or a directory keep only its path - independently of the current dir
# my $path = path($filename);
#----------------------------------------------------------------------------
sub get_path {
	my($file) = shift;
	($file =~ /\//)?($file =~ s/(.*)\/.*$/$1/):($file = ".");
	return $file;
}

#----------------------------------------------------------------------------
# get TE infos
# $TE = get_TEs_infos($TEclass)
#----------------------------------------------------------------------------
sub get_TEs_infos {
	my $input = shift; #file name
	my %TEs = ();
	open(my $input_fh, "<", $input) or confess"\n      ERROR (sub get_TEs_infos): could not open $input $!\n\n";
	LINE: while(<$input_fh>) {
		chomp (my $line = $_);
		next LINE if ($line !~ /\w/);
		my @TEs = split('\t', $line); 
		my $lcRname = lc ($TEs[0]);
		$TEs{$lcRname} = \@TEs;
	}	
	return \%TEs;
}

#----------------------------------------------------------------------------
# get species list, ordered
# my $sp_ordered = get_species_list(\@Genomes,$tree,$v);
# called by main
#----------------------------------------------------------------------------
sub get_species_list {
	my ($G,$tree,$v) = @_;
	
	#get list of species using list of genomes => used for print of the output
	my @species = ();
	foreach my $fa (@{ $G } ) {
		chomp $fa;
		my $sp = filename($fa);
		$sp = $1 if ($sp =~ /^(.+?)\.fa/);
		push(@species,$sp);
	}

	#Now order species in a specific way if $tree info
	my @sp_ordered = ();
	if ($tree ne "na") {
		my @sp_tree = split(",",$tree) ;	
		#check that @sp_tree contains all species that are in the directory
		#just compare their dimension if needed
		return if ($#species != $#sp_tree);
		#if OK, carry on
		print STDERR "        ordering species based on tree info provided\n" if ($v);
		@sp_ordered = @sp_tree;
	} else {
		print STDERR "        ordering species alphabetically\n" if ($v);
		@sp_ordered = sort { $a cmp $b } @species; #if no tree info just do alphabetically
	}
	return \@sp_ordered;
}

#----------------------------------------------------------------------------
# Get lengths of all sequences and store that by sequence ID. Note that if some are not unique, it just replaces by last length.
# This sub does not use BioPerl - avoid having to index the genome
# $genlen = get_lengths_noBio($fa,$genlen,$v);
# called by main
#----------------------------------------------------------------------------
sub get_lengths_noBio {
	my ($fa,$genlen,$v) = @_;
	my $lengthfile = "$fa.lengths";		
	if (-e $lengthfile) {
		print STDERR "         -> sequences lengths have been previously calculated ($lengthfile exists) => extracting\n" if ($v);
		#extract lengths now
		open (my $lengths_fh, "<", $lengthfile) or confess "   \nERROR (sub get_lengths_noBio): could not open $lengthfile $!\n";
		while (<$lengths_fh>) {
			chomp (my $line = $_);
			my ($id,$len) = split(/\s+/,$line);
			$genlen->{$fa}{$id}=$len;
		}	
		close ($lengths_fh);
	} else {
		#looping through fasta file
		print STDERR "         -> sequences lengths are being extracted from $fa\n" if ($v);
		my $id = "";
		my $l = 0;
		my $c = 0;
		open (my $fa_fh, "<", $fa) or confess "   \nERROR (sub get_lengths_noBio): could not open $fa $!\n";
		open (my $len_fh, ">", $lengthfile) or warn "   \nERROR (sub get_lengths_noBio): could not create $lengthfile, but lengths will be calculated $!\n";
		while (<$fa_fh>) {
			chomp (my $line = $_);
			if (substr($line,0,1) eq ">") {
				#first get and print unless first header
				unless ($c == 0) {
					print $len_fh "$id\t$l\n";
					$genlen->{$fa}{$id}=$l;
				}
				$c=1;
				#store header and reinitialize length
				my @id = split (/\s+/,$line);
				$id = $id[0];
				$id =~ s/>//;
				$l = 0;
			} else {
				#get length; could be more than one line so increment
				$l+=length($line);
			}
		}
		#get and print len last sequence
		print $len_fh "$id\t$l\n";
		$genlen->{$fa}{$id}=$l;
		
		close ($fa_fh);
		close ($len_fh);
	}	
	return ($genlen);
}

#----------------------------------------------------------------------------
# check for genome
# my $list = check_for_genome(\@RMout_files,$dir,$v)
#----------------------------------------------------------------------------
sub check_for_genome {
	my ($files,$dir,$v) = @_;
	my @list = ();
	foreach my $file (@$files) {
		chomp $file;
		my $fa = filename($file);
		$fa = $1.".fa" if $fa =~ /^(.*)\.out/;
		unless (-e "$dir/$fa") {
			print STDERR "     ERROR (sub check_for_genome): $fa doesn't exist - check genome file name (see usage)\n" if ($v);
			print STDERR "       -> $file skipped for analysis\n" if ($v);
		} else {
			push(@list,$file);
		}
	}
	return (\@list);
}

#----------------------------------------------------------------------------
# split RMoutput files and filter them - so that it's faster to load in arrays + threading
# my $split_RMouts = split_RM_files($dir,\@RMout_files,\@flist,$filter,$TE,$TEclass,$f_regexp,$v); 
#----------------------------------------------------------------------------
sub split_RM_files {
	my ($dir,$RMouts,$filter,$flist,$TE,$TEclass,$f_regexp,$v) = @_;	
	my @list = ();

	#get filter if relevant
	my ($f_type,$f_name) = split(",",$filter) unless ($filter eq "na");
	my %check = ();	
	#check that corresponding genome exists, if not exclude from list
	my $RMs = check_for_genome($RMouts,$dir,$v);
	
	FILE: foreach my $file (@{$RMs}) {
		chomp $file;
		next FILE unless (-e $file); #for some reason there was an extra empty value in the list. Better check on file existence.
		print STDERR "       -> $file...\n" if ($v);

		#get genome name
		my $fa = $file;
		$fa = $1.".fa" if $fa =~ /^(.*)\.out/;
		
		#make outputdir
		my $filename = filename($file);
		my $splitpath = $dir."/".$filename.".split";
		unless (-e $splitpath) {
			system "mkdir $splitpath"; #or die "\n      ERROR (sub split_RM_files): can not mkdir $splitpath $!\n\n";
		}	
		
		#put RMout in array to be able to loop on it (because need to access to previous and next $block)
		my @RMout = get_RMout_array($file);		
		#now loop on array
		my %check = ();
		RMLINE: for (my $i = 0; $i <= $#RMout; $i++){
			my ($Rname,$classfam,$block) 
			   = ($RMout[$i]->[9],$RMout[$i]->[10],$RMout[$i]->[14]);	
			next RMLINE unless ($block); #apparently, that happens	
			my ($Rclass,$Rfam,$Rclassfam) = get_Rclass_Rfam($classfam);
						 			
 			#filter if fragment is nested (as detected by block ID).
			unless (($i == 0) || ($i == $#RMout)) {
				my $Ublock = $RMout[$i-1]->[14];
				my $Dblock = $RMout[$i+1]->[14];
				next RMLINE if (($Dblock) && ($Ublock) && ($Ublock eq $Dblock) && ($block ne $Ublock));
			}	
 			
			#Deal with class and family, if replacement required
			if ($TEclass ne "na") {
				if ($TE->{lc($Rname)}) {
					$Rclass = $TE->{lc($Rname)}->[1];
					$Rfam = $TE->{lc($Rname)}->[2];
					$Rclassfam = $TE->{lc($Rname)}->[3];	
				} else {
					print STDERR "       !! $Rname not found in $TEclass => RM output class / fam used instead\n" if (($TEclass) && ($v) && (! $check{lc($Rname)}{"class"}));
					$check{lc($Rname)}{"class"}=1;
				}	
			} 
			#filter out non TE stuff
			next RMLINE if (($Rclass eq "nonTE") || ($Rclass eq "Simple_repeat") || ($Rclass eq "Low_complexity") 
			               || ($Rclass eq "Satellite") || ($Rclass =~ /RNA$/) || ($Rclass =~ /omeric$/));
		
		
			#filter out stuff from -filter if relevant
			my ($lcRname,$lcRclass,$lcRfam,$lcf_name) = (lc($Rname),lc($Rclass),lc($Rfam),lc($f_name));
			if ($filter ne "na") {				
				if ($f_regexp eq "y") {
					#check if what's det as the filter is included in the names
					next RMLINE unless ((($f_type eq "name") && ($lcRname =~ /$lcf_name/))
					               || (($f_type eq "class") && ($lcRclass =~ /$lcf_name/)) 
					               || (($f_type eq "family") && ($lcRfam =~ /$lcf_name/)));
				} else {
					next RMLINE unless ((($f_type eq "name") && ($lcf_name eq $lcRname))
					               || (($f_type eq "class") && ($lcf_name eq $lcRclass)) 
					               || (($f_type eq "family") && ($lcf_name eq $lcRfam)));
				}	
			}
			
			#filter out stuff from -list if relevant
			if ($flist->[0] ne "na") {	
				FLIST: foreach my $name (@{$flist}) {
					my $lcname = lc($name);	
					next RMLINE if (($f_regexp eq "y") && ($lcname !~ /$lcRname/));
					next RMLINE if (($f_regexp eq "n") && ($lcname ne $lcRname));
					last FLIST if (($f_regexp eq "y") && ($lcname =~ /$lcRname/));
					last FLIST if (($f_regexp eq "n") && ($lcname eq $lcRname));
				}
			}		
												
			#now print in split files what went through
			my $split = $splitpath."/".$filename."##$Rname";
			push(@list,$split) unless (-e $split);
			my $line = "$RMout[$i]->[0]\t$RMout[$i]->[1]\t$RMout[$i]->[2]\t$RMout[$i]->[3]\t$RMout[$i]->[4]\t$RMout[$i]->[5]\t$RMout[$i]->[6]\t$RMout[$i]->[7]\t$RMout[$i]->[8]\t$RMout[$i]->[9]\t$RMout[$i]->[10]\t$RMout[$i]->[11]\t$RMout[$i]->[12]\t$RMout[$i]->[13]\t$RMout[$i]->[14]\n";
			open (my $fh_split,">>",$split) or confess "\n      ERROR (Sub split_RM_files): Failed to open $split $!\n\n";
			print $fh_split "$line";
			close($fh_split);
		}
		print STDERR "       ...$file done\n" if ($v);
	}
	return \@list;
}

#----------------------------------------------------------------------------
# load RM output in array of array
# my @RMout = get_RMout_array($RMout);
#----------------------------------------------------------------------------
sub get_RMout_array {
	my $rm = shift;
	my @array = ();
	open my $rm_fh, "<", "$rm" or confess "\n      ERROR (Sub get_RMout_array): can't open $rm $!\n\n";
	LINE: while(<$rm_fh>) {
		chomp (my $line = $_);
		next LINE if (($line =~ /position|[sS]core|Gname/) || ($line !~ /\w/)); #skip headers and white lines
		$line =~ s/^\s+//; #remove spaces in beginning of lines
		my @line = split(/\s+/,$line);
		push(@array, \@line);
	}	
	return @array;
}

#----------------------------------------------------------------------------
# Get Rclassfam from RMout
# my ($Rclass,$Rfam,$Rclassfam) = get_Rclass_Rfam($classfam);
#----------------------------------------------------------------------------
sub get_Rclass_Rfam {
	my($classfam) = shift;
	my ($Rclass,$Rfam);
	if ($classfam =~ /\//) {
		#deal if several / in name for some reason: replace only first /
		$classfam =~ s/\//%/;
		($Rclass,$Rfam) = split("%", $classfam);
	} else {
		$Rfam = $classfam;
		$Rfam=~ s/^(.*)\..*$/$1/;
		$Rclass = $classfam;
		$Rclass =~ s/^.*\.(.*)$/$1/;
	}
	my $Rclassfam = "$Rclass/$Rfam";
	return ($Rclass,$Rfam,$Rclassfam);
}

#----------------------------------------------------------------------------
# Check the name of the repeat => if file should be processed; 
# my $return = check_Rname($RMout,$filter,$flist,$f_regexp) if ($$presplit == 1);
#----------------------------------------------------------------------------
sub check_Rname {
	my($RMout,$filter,$flist,$f_regexp) = @_;
	#deal with filter
	my ($f_type,$f_name) = split(",",$$filter) unless ($$filter eq "na");
	my $lcf_name = lc($f_name);
		
	#get repeat name from file name
	$RMout =~ s/^.+?##(.+?)$/$1/;
	my $lcRname = lc($RMout);
		
	#Now check if goes through filter
	#1) -filter	
	if (($$filter ne "na") && ($f_type eq "name"))  {
		return ("yes") if (($$f_regexp eq "y") && ($lcRname =~ /$lcf_name/)); #do not skip if it contains 
		return ("yes") if (($$f_regexp eq "n") && ($lcRname eq $lcf_name)); #do not skip if it matches
	}	
	#2) -list	
	if ($flist->[0] ne "na") {	
		FLIST: foreach my $name (@{$flist}) {
			my $lcname = lc($name);
			if ((($$f_regexp eq "y") && ($lcname =~ /$lcRname/)) || (($$f_regexp eq "n") && ($lcname eq $lcRname))) {
				return ("yes"); #do not skip if it contains 
				last FLIST;
			} 
		}
	}
	return ("no"); #if no match then 1 means skip it
}

#----------------------------------------------------------------------------
# Get position files to extract repeats
# ($posi_list,$posi) = get_posi_from_RMout($RMout,$frgs_all,$frgs_nr,$posi_list,$append,$v);
#----------------------------------------------------------------------------
sub get_posi_from_RMout {
	my ($RMout,$frgs_all,$frgs_nr,$posi_list,$append,$v) = @_; 
	
	print STDERR "       ..in progress: get positions / fragment numbers from $RMout (thr ".threads->tid().")\n" if ($$v);
	
	#put RMout in array to be able to loop on it (because need to access to previous and next $block)
	my @RMout = get_RMout_array($RMout);
	
	#now loop on array
	my %blockcheck = ();
	my %c = ();
	my $posi = $RMout.".tab";
	push(@{$posi_list},$posi);
	open (my $fh, ">>","$posi") or confess "\n      ERROR (Sub get_posi_from_RMout): Failed to create/open to write $posi $!\n\n" unless (($$append eq "yes") && (-e $posi));
	#Note that all filters were previously done in the split step	
	for (my $i = 0; $i <= $#RMout; $i++){
		my ($div,$Gname,$Gstart,$Gend,$strand,$Rname,$classfam,$block) 
		   = ($RMout[$i]->[1],$RMout[$i]->[4],$RMout[$i]->[5],$RMout[$i]->[6],$RMout[$i]->[8],$RMout[$i]->[9],$RMout[$i]->[10],$RMout[$i]->[14]);		
 		my ($Rclass,$Rfam,$Rclassfam) = get_Rclass_Rfam($classfam);		
 		
		#count all frgs
		my $posiname = filename($posi);
		($frgs_all->{$posiname})?($frgs_all->{$posiname}+=1):($frgs_all->{$posiname}=1);
		
		#count number of fragments, corrected by blockID. Shared, so can't be complex data...
		my $fullblock = $Gname."_".$block;
		unless ($blockcheck{$Rname.$fullblock}) {
			($frgs_nr->{$posiname})?($frgs_nr->{$posiname}+=1):($frgs_nr->{$posiname}=1);
		}	
		$blockcheck{$Rname.$fullblock}=1; #keep in mind so that counts of fragments are corrected
		
		#now print coordinates for each copy - each file will be accessed by only one thread so no need to use Semaphore	
		my $cc = $Rname."_".$block."_".$Gname;
		($c{$cc})?($c{$cc}++):($c{$cc}=1);
		print $fh "$Gname\t$Gstart\t$Gend\t$strand\t$Rclass\t$Rfam\t$Rclassfam\t$div\t$block\t$c{$cc}\n" unless (($$append eq "yes") && (-e $posi));
	}	
	close($fh) unless (($$append eq "yes") && (-e $posi));
	print STDERR "       ..done: get positions / fragment numbers from $RMout (thr ".threads->tid().")\n" if ($$v);
	return ($posi_list,$posi);
}

#----------------------------------------------------------------------------
# Get fragment numbers from posi files already generated
# my ($frgs_all,$frgs_nr) = get_frgs_fromposi($dir,$v);
#----------------------------------------------------------------------------
sub get_frgs_fromposi {
	my ($dir,$v) = @_; 	
	print STDERR " --- in progress: get fragment numbers from posi files\n" if ($v);
	my @posi = `find $dir/*split -type f -name "*.out*##*" | grep ".tab"` or die "\n      ERROR (main): can't find posi files *.out*##*.tab in $dir/*split $!\n\n";
	my $frgs_all;
	my $frgs_nr;
	my %blockcheck = ();
	PF: foreach my $p (@posi) {
		chomp ($p);
		next PF unless ($p);
		if (! -e $p) {
			print STDERR "     ERROR: posi file $p does not exist?\n" if ($v);
			next PF;
		}
		my $pname = filename($p);
		open (my $fh,"<",$p) or confess "\n      ERROR (Sub get_posi_from_RMout): Failed to open to read $p $!\n\n"; 
		while(<$fh>) {
			my @l = split("\t",$_);
			my $fullblock = $l[0]."_".$l[8]; #Gname and block
			($frgs_all->{$pname})?($frgs_all->{$pname}+=1):($frgs_all->{$pname}=1); #each line is a fragment
			unless ($blockcheck{$fullblock}) {
				($frgs_nr->{$pname})?($frgs_nr->{$pname}+=1):($frgs_nr->{$pname}=1);
			}	
			$blockcheck{$fullblock}=1;
		}		
		close $fh;
	}
	print STDERR "     ...done\n" if ($v);
	return($frgs_all,$frgs_nr);
}

#----------------------------------------------------------------------------
# Sub to loop on big hash to extract sequences
# extract_sequences($dir,$faloc,$posi,$flank,$min_len,$frgs_all,$genlen,$append,$v);	
#----------------------------------------------------------------------------
sub extract_sequences {
	my ($dir,$faloc,$posi,$flank,$min_len,$frgs_all,$genlen,$append,$v) = @_;
	
	print STDERR "       ..in progress: extract sequences from $posi (thr ".threads->tid().")\n" if ($$v);
		
	#get genome (existence already checked)
	my $sp = filename($posi);
	$sp = $1 if $sp =~ /^(.*)\.out/;
	my $fa = $sp.".fa";
	$fa = $$dir."/".$fa;
	
	#check if sequences should be extracted or not
	my $checkseq = 0;
	if ($$append eq "yes") {
		my @exist = `ls $faloc`; #or confess "\n      ERROR (sub extract_sequences): can't list files in $faloc $!\n\n";
		if ($#exist > 0) {
			EXIST: foreach my $exists (@exist) {
				chomp $exists;
				$exists = $1.".fa" if ($exists =~ /^(.+?)\./);
				my $faname = filename($fa);
				$checkseq=1 if ($exists eq $faname);
				last EXIST if ($checkseq == 1);
			}
		}
	}	
	return if ($checkseq == 1);
		
	# connect to the fasta file
	my $db = Bio::DB::Fasta->new($fa, -reindex => 0) or confess "\n      ERROR (Sub extract_sequences): Failed to open Bio::DB::Fasta object from $fa $!\n\n";
		
	#read posifile to extract sequences - one fasta file per element
	# NB: my $posi = $RMout."##".$Rname.".tab"; # $RMout = RMout first provided; otherwise, with $RMout from thread it's $RMout.tab
	my $posiname = filename($posi);
	my ($RMout,$Rname) = ($1,$2) if $posiname =~ /^(.+?\.out.*?)##(.+?)\.tab$/;
	#deal with output fasta file with extracted sequences - 1 file per repeat for each RMout/genome
	my $extracted = $faloc."/".$sp.".".$Rname.".extract.fa";
	
	#Now deal with extraction
	my $check = ();
	my %prevcoords = ();
	my @seqs = ();
	my $s = 0;
	my $counts = 0;
	open (my $posi_fh, "<", $posi) or confess "\n      ERROR (Sub extract_sequences): Failed to open $posi $!\n\n";
	REP: while (<$posi_fh>) {
		chomp (my $l = $_);
		my ($Gname,$Gstart,$Gend,$strand,$Rclass,$Rfam,$Rclassfam,$div,$block,$c) = split ("\t",$l);
		$counts++;		
		#If first copy of any given block is processed that means print what is stored + store
		if ($c == 1) {
			#extract previous unless it is first line and nothing is stored; temp values because $seqs and $div may not be defined.
			my ($seqs,$check) = extract_one_seq(\%prevcoords,$sp,$Rname,$db,$flank,$min_len,$genlen,$fa,$check,$v) if ($prevcoords{1});
			$seqs[$s]=$seqs if ($seqs);
			#extraction done so increment the array coordinate in $seqs
			$s++ if ($seqs);
			#empty hashes and store this current line as the new previous
			%prevcoords = ();
		} 
		#Now store lines, in hash that has been reset or not (eg $c = 1 others)
		$prevcoords{$c} = "$Gname\t$Gstart\t$Gend\t$strand\t$Rclassfam\t$div";
		
		#Check for last line whatever the number
		my ($seqs,$check) = extract_one_seq(\%prevcoords,$sp,$Rname,$db,$flank,$min_len,$genlen,$fa,$check,$v) if ($counts == $frgs_all->{$posiname});
		$seqs[$s]=$seqs if ($seqs);
		
		#print $seqs if more than 500 to avoid big memory usage
		my $seqnb = @seqs;
		if ($seqnb > 499) {
			open (my $fh, ">>","$extracted") or confess "\n      ERROR (Sub extract_sequences): Failed to create/write in $extracted $!\n\n";		
			for (my $i = 0; $i <= $#seqs; $i++) {
				print $fh ">".$seqs[$i]."\n" if ($seqs[$i]);
			} 
			close($fh);
			#reinitialize
			@seqs = ();
			$s = 0;
		}
	}
	close($posi_fh);
	
	#now print whatever is left in @seqs and was not printed, unless specified not to
	my $seqnb = @seqs;
	if ($seqnb > 0) {
		open (my $fh, ">>","$extracted") or confess "\n      ERROR (Sub extract_sequences): Failed to create/write in $extracted $!\n\n";		
		for (my $i = 0; $i <= $#seqs; $i++) {
			print $fh ">".$seqs[$i]."\n" if ($seqs[$i]);
		} 
		close($fh);
	}
	print STDERR "       ..done: extract all sequences from $posi (thr ".threads->tid().")\n" if ($$v);
	undef($db);
	undef(@seqs);
	undef(%prevcoords);
	return;	
}			
			
#----------------------------------------------------------------------------
# Extract a sequence using coordinates 
# ($seqs,$check) = extract_one_seq(\%prevcoords,$sp,$Rname,$db,$flank,$min_len,$genlen,$fa,$check,$v);
#----------------------------------------------------------------------------
sub extract_one_seq {					
	my ($coords,$sp,$Rname,$db,$flank,$min_len,$genlen,$fa,$check,$v) = @_;
	
	#get list of the ID of the genome file
	my @dbIDs = $db->get_all_ids();
	
	my $Gs;
	my $Ge;
	my $len;
	my $nb = 0;	
	my $div_corr;
	my $pondiv;
	#First round to get start and end, if more than one fragment per block, for the name of the sequence
	foreach my $c (sort keys %{ $coords }) {
		my ($Gname,$Gstart,$Gend,$strand,$Rclassfam,$div) = split ("\t",$coords->{$c});
		my $currlen = $Gend - $Gstart + 1;
		#Correct coords for new name + get real length
		$div_corr = $div;
		if ($c == 1) {
			$Gs = $Gstart;
			$Ge = $Gend; #to intialize
			$len = $currlen;
			$pondiv = $div*$currlen;
		} else {
			$Ge = $Gend; #to correct - Gs remains the same
			$len+=$currlen;
			$pondiv+=$div*$currlen;
		}
		$nb++;
	}
	#correct %div if $nb >1
	$div_corr = $pondiv/$len if ($nb > 1);
	
	#no extraction if too small
	if ($$min_len ne "na") {
		return if ($len < $$min_len);
	}
	
	#Now loop to extract sequence(s)
	my $newId = "";
	my $descr = "";
	my $seq = "";
	my $fl5;
	my $fl3;
	foreach my $c (keys %{ $coords }) {
		my ($Gname,$Gstart,$Gend,$strand,$Rclassfam,$div) = split ("\t",$coords->{$c});
		
		#name of the sequence extracted will be only gb when there is gi and gb in name (otherwise too long)
		my $Gnamemod = $Gname;
		$Gnamemod =~ s/^gi\|.*\|gb/gb/;
		
		#prep new sequence name
		$newId = $Rname."__".$nb."_".$Gnamemod."--".$Gs."-".$Ge; #avoid the : for later potential uses of Bio::DB::Fasta	
		$descr = $Rclassfam."__frg=".$nb."__len=".$len."_st=".$strand."_div=".$div_corr."_sp=".$sp;
		
		#Check if $Gname is in $db and get $Glen at the same time or exit the sub
		my $Glen;
		if ($genlen->{$fa}{$Gname}) {
			$Glen = $genlen->{$fa}{$Gname}; 
		} else {
			print STDERR "        ERROR (Sub extract_one_seq): $Gname was not found in $fa\n" if ((! $check->{$fa}{$Gname}) && ($v));
			$check->{$fa}{$Gname}=1;
			return ($seq,$check);
		}	
		
		#get coords with flanks, but only if outside fragments 
		if ($c==1) { #first fragment => get flank at start
			$Gstart = $Gstart-$$flank;
			$Gstart = 1 if ($Gstart < 1); #correct to avoid negative start
		} 
		if ($c == $nb) { #last fragment => get flank at end
			$Gend = $Gend+$$flank;
			$Gend = $Glen if ($Gend > $Glen); #correct to avoid getting out of scaffold
		}
		
		# now extract target sequence, with or without flankings
		my $subSeq = $db->seq($Gname,$Gstart,$Gend);
		if ($c == 1){
			$seq = $subSeq; #get the sequence of the first fragment
		} else {
			$seq = $seq.$subSeq; # add other piece(s)
		}
		
		# get actual flankings if relevant
		if ($c == 1) {
			$fl5 = $Gs - $Gstart; #real TE start minus start of flanking => real 5' flanking size
			$fl3 = $Gend - $Ge; #Do to avoid errors Use of uninitialized value $fl3 - Value is corrected when it gets to last frg
		}
		if ($c == $nb) {
			$fl3 = $Gend - $Ge; #new end after fl determination minus actual end => real 3' flanking size
		}
	}	
	#correct desc to add flanking length
	$descr = $descr."_5fl=".$fl5."_3fl=".$fl3;	
	#get full sequence now
	my $fullseq = "$newId\t$descr\n$seq";	
	#now return it
	return ($fullseq,$check);	
}

#----------------------------------------------------------------------------
# Blast all sequences in a file against all genomes in $dir
# blast($dir,$faloc,$blastloc,$evalue,$append,$v);
# called by main
#----------------------------------------------------------------------------
sub blast {					
	my ($dir,$faloc,$blastloc,$evalue,$append,$v) = @_;
	print STDERR "       ..in progress: blastn extracted sequences located in $faloc (thr ".threads->tid().")\n" if ($$v);

	my @fafiles = `ls $faloc/*.fa` or confess "\n      ERROR (sub blast): can't list files .fa in $faloc $!\n\n";
	my @genomes = `ls $$dir/*.fa` or confess "\n      ERROR (sub blast): can't list files .fa in $dir $!\n\n";
	my $skipped = 0;
	my $newe = $evalue;
	foreach my $g (@genomes) {
		chomp ($g);
		chomp (my $fa = $fafiles[0]);
		my $gname = filename($g);
		$gname = $1 if ($gname =~ /^(.+)\.fa$/);
		
		#lower evalue for same species blast
		my $faname = filename($fa);
		$newe = "10e-120" if ($faname =~ /^$gname/); #Is same start
		
		#now blast unless append and output exists
		my $out = "$fa.blast.$gname.out";
		if (($$append eq "yes") && (-e $out)) {
			print STDERR "          - skipped (because -append and output exists): blast $fa against $g\n" if ($v);		
		} else {
			print STDERR "          - $fa against $g (thr ".threads->tid().")\n" if ($v);
# 			print STDERR "         Blast cmd line = $blastloc/blastn -query $fa -out $fa.blast.$gname.out -db $g -evalue $evalue -num_alignments 1\n" if ($$v);
		`$blastloc/blastn -query $fa -out $out -db $g -evalue $newe -max_target_seqs 1 -outfmt "6 qseqid qstart qend qlen qcovs qcovhsp sseqid sstart send sstrand slen evalue bitscore score length pident nident mismatch ppos positive gapopen gaps"`;
		}	
	}	
	print STDERR "       ..done: blast all sequences from $faloc against provided genomes (thr ".threads->tid().")\n" if ($$v);
	return;
}

#----------------------------------------------------------------------------
# Get Qdescs from fa file that was used to blast
# $Qdesc = get_Qdesc_from_extracted_seqs($faloc,$v);
#----------------------------------------------------------------------------
sub get_Qdesc_from_extracted_seqs {					
	my ($loc) = @_;
	my $Qdesc;
	#Qdescr - because blast output does not keep that info
	my @seqs = `ls $loc/*.fa` or confess "\n      ERROR (sub get_Qdesc_from_extracted_seqs): can't list files .fa in $loc $!\n\n";
	foreach my $seq (@seqs) {
		chomp $seq;
		my $fa = Bio::SeqIO->new(-file => $seq, -format => "fasta") or confess "\n      ERROR (sub get_Qdesc_from_extracted_seqs): Failed to open $seq $!\n\n";
		while( my $s = $fa->next_seq() ) {
			my $id = $s->display_id;
			my $desc = $s->desc;
			$Qdesc->{$id}=$desc;
		}
	}	
	return ($Qdesc);	
}

#----------------------------------------------------------------------------
# parse blast output
# my ($ortho,$tothit) = parse_blast($faloc,$flank,$v);
# called by main
#----------------------------------------------------------------------------
sub parse_blast {
	my ($outpath,$flank,$v) = @_;
	$flank = $flank/100*50; #if flank of 100 nt, 50 nt are required to decide it is orthologous
	my %ortho = ();
	my %tothit = ();
#	my @blastout = `ls $outpath/DATA/*.extract/*.blast.*.out`;
	my @blastout = `find $outpath/DATA/*.extract/ -name *.blast.*.out`; #Check to see if this fix error sh: /bin/ls: Argument list too long
	foreach my $Bout (@blastout) {
		chomp $Bout;
# 		print STDERR "     - parsing $Bout\n" if ($v);
		my $Bfile = filename($Bout);
		my $loc = get_path($Bout);
		#Fetch Query descriptions from fa files
		#print STDERR "     - Getting sequence descriptions by opening fa files from $loc\n" if ($v);
		my $Qdesc = get_Qdesc_from_extracted_seqs($loc);
		#Now parse
		my ($Q_sp,$T_sp) = ($1,$2) if ($Bfile =~ /^(.+?)\..*\.blast\.(.+)\.out$/);
		
		#Structure of the tabular file:
		#qseqid qstart qend qlen qcovs qcovhsp sseqid sstart send sstrand slen evalue bitscore score length pident nident mismatch ppos positive gapopen gaps
# 
#    	    qseqid means Query Seq-id
#    	    qstart means Start of alignment in query
#    	      qend means End of alignment in query   	     	    
#    	      qlen means Query sequence length
# 
#          qcovs means Query Coverage Per Subject
#    	   qcovhsp means Query Coverage Per HSP
#    	   
#    	    sseqid means Subject Seq-id
#    	    sstart means Start of alignment in subject
#    	      send means End of alignment in subject
#    	      slen means Subject sequence length
#    	   sstrand means Subject Strand
# 
#    	    evalue means Expect value
#    	  bitscore means Bit score
#    	     score means Raw score
#    	    length means Alignment length
#    	    pident means Percentage of identical matches
#    	    nident means Number of identical matches
#    	  mismatch means Number of mismatches
#    	      ppos means Percentage of positive-scoring matches
#    	  positive means Number of positive-scoring matches
#    	   gapopen means Number of gap openings
#    	      gaps means Total number of gaps
		
		my $Boutp = "$Bout.ortho.tab";
		open(my $Boutp_fh,">",$Boutp) or confess "\n      ERROR (Sub parse_blast): Failed to open to write $Boutp $!";	
		print $Boutp_fh "#query_id\tstart_in_query\tend_in_query\tquery_len\tquery_coverage_per_hit\tquery_coverage_per_hsp\thit_id\tstart_in_hit\tend_in_hit\thit_strand\thit_len\tevalue\tbit_score\traw_score\talignment_length\t%_identical_matches\tnumber_identical_matches\tnumber_mismatch\t%_positive-scoring_matches\tnumber_positive-scoring_matches\tnumber_gap_opening\tnumber_total_gaps\n\n";
		open(my $Bout_fh,"<",$Bout) or confess "\n      ERROR (Sub parse_blast): Failed to open to read $Bout $!";	
		while(<$Bout_fh>) {
			chomp (my $line = $_);			
			my ($Qname,$start,$end,$Rlen,$qcovs,$qcovhsp,$sseqid,$sstart,$send,$sstrand,$slen,$evalue,$bitscore,$score,$length,$pident,$nident,$mismatch,$ppos,$positive,$gapopen,$gaps) = split (/\t/,$line);			
			my $Rname = $1 if ($Qname =~ /^(.+?)__/);
			my ($Rclassfam,$Rinfo) = split(/__/,$Qdesc->{$Qname});
			my $Rfullname = $Rname."#".$Rclassfam;
			#check boundaries to assess if potentially orthologous, with a 10nt security
			if (($start < $flank) && ($end > ($Rlen - $flank))) {
				($ortho{$Rfullname}{$Q_sp}{$T_sp})?($ortho{$Rfullname}{$Q_sp}{$T_sp}++):($ortho{$Rfullname}{$Q_sp}{$T_sp}=1);
				print $Boutp_fh "$line\n";
			}
			($tothit{$Rfullname}{$Q_sp}{$T_sp})?($tothit{$Rfullname}{$Q_sp}{$T_sp}++):($tothit{$Rfullname}{$Q_sp}{$T_sp}=1);
		}
		close ($Bout_fh);
		close ($Boutp_fh);
	}
	return (\%ortho,\%tothit);
}	

#----------------------------------------------------------------------------
# get orthology output
# get_ortho_output($dir,$outpath,$sp_ordered,$tree,$ortho,$tothit,\%frgs_all,\%frgs_nr,$v);
# called by main
#----------------------------------------------------------------------------
sub get_ortho_output {
	my ($dir,$outpath,$species,$tree,$ortho,$tothit,$frgs_all,$frgs_nr,$v) = @_;
	
	my $out = "$outpath/_TE_Orthology.results.tab";
	open(my $fh,">",$out) or confess "\n      ERROR (Sub get_ortho_output): Failed to create $out $!";
	
	#Prep the output file headers, with tree info if any
	print $fh "#\"Total_hits\" means total of hits with evalue less than the one set through -e\n";
	print $fh "#\"ortho_hits\" means POTENTIALLY orthologous elements = hits that contain at least 20% of the extracted flanking sequence\n";
	print $fh "#Acess to individual copies through files in $outpath/DATA/*.extract/*.out.ortho.tab\n\n";
	print $fh "#\t[Query_species]\t.\tTarget_species";
	unless ($tree eq "na") {
		print $fh "#(with tree info):\n\t\t\t";
		foreach my $sp (@{ $species }) {
			print $fh "$sp\t,\t";
		}
	}
	print $fh "\n#\t\t\t";
	foreach my $sp (@{ $species }) {
		$sp =~ s/\(|\)//g if ($tree ne "na"); #replace parentheses if it was tree info
		print $fh "$sp\t$sp\t";
	}
	print $fh "\n#Rname\t[name]\t[Nb_extracted_seq]\t";
	foreach my $sp (@{ $species }) {
		print $fh "Total_hits\tOrtho_hits\t";
	}
	print $fh "\n";
	
	#now print output
	foreach my $Rfullname (sort keys %{ $ortho }) {
		foreach my $query (sort keys %{ $ortho->{$Rfullname} }) {
			my $Rname = $1 if ($Rfullname =~ /^(.+?)#/);
			my $key = `ls $dir/*split/$query.out*##$Rname.tab`;
			chomp $key;
			$key = filename($key);
			print $fh "$Rfullname\t$query\t".$frgs_nr->{$key}."\t";
			#now follow same order of species pre-established
			foreach my $target (@{ $species }) {
				$target =~ s/\(|\)|;//g if ($tree ne "na"); #replace parentheses if it was tree info
				$ortho->{$Rfullname}{$query}{$target}=0 unless ($ortho->{$Rfullname}{$query}{$target}); #to avoid blanks in print
				$tothit->{$Rfullname}{$query}{$target}=0 unless ($tothit->{$Rfullname}{$query}{$target}); #to avoid blanks in print
				print $fh $tothit->{$Rfullname}{$query}{$target}."\t".$ortho->{$Rfullname}{$query}{$target}."\t";
			}
			print $fh "\n";
		}
	}
}

#----------------------------------------------------------------------------
# get orthology output
# just_parse($outpath,$flank,$dir,$sp_ordered,$tree,$v);
# called by main
#----------------------------------------------------------------------------
sub just_parse {
	my ($outpath,$flank,$dir,$sp_ordered,$tree,$v) = @_;
	print STDERR "\n --- -ponly chosen, so only action will be to parse some previously obtained outputs\n" if ($v);	
	my ($frgs_all,$frgs_nr) = get_frgs_fromposi($dir,$v);
	
	
	#Now parse blast
	print STDERR "\n --- Parse blast outputs (check for presence of flankings in the best hit)\n" if ($v);
	my ($ortho,$tothit) = parse_blast($outpath,$flank,$v);

	#Now print
	print STDERR "\n --- Printing results (orthology)\n" if ($v);
	get_ortho_output($dir,$outpath,$sp_ordered,$tree,$ortho,$tothit,\%frgs_all,\%frgs_nr,$v);

	return;
}