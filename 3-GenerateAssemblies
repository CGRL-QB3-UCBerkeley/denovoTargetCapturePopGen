#!/usr/bin/perl

#Ke Bi (kebi@berkeley.edu)

use strict;
use warnings;
use Getopt::Std;
use Getopt::Long;
use File::Basename;



&main;
exit;

sub main {
        &usage if (@ARGV<1);
        my $command = shift(@ARGV);
        my %fun = (abyss=>\&abyss,trinity=>\&trinity, soapdenovo=>\&soapdenovo);
        die("Unknown command \"$command\"\n") if (!defined($fun{$command}));
        &{$fun{$command}};
      }


sub usage {
  die(qq/
Usage: 3-GenerateAssemblies <command> [<arguments>]\n

Command: 

trinity: Assembly using trinity (for RNAseq)

abyss: Assembly using ABySS (for exon capture)

\n/);
}

sub trinity {
  die(qq/
Usage 3-generateAssemblies_LOCAL.pl trinity [options]

Options: -a library folder
         -c JM memory [32]
         -d min_kmer_cov [2]
         -e cpus [6]

Assuming the unpaired reads (XXX_u_final.txt) are 
merged with the left end reads (XXX_1_final.txt), 
with a new name of "XXX_1_final.txt".
  
\n/) if (!@ARGV);
  
  my %opts = (a=>undef,  c=>32, d=>2, e=>6);
  getopts('a:c:d:e', \%opts);
  
  
  

  my $dir;
  
  if ($opts{a} =~ m/\/$/ ){
    $dir = $opts{a}; 
  }
  else {
    $dir = $opts{a} . "/";
  }
  
   
  
  my @files = <$dir*_1_final.txt>;
  foreach my $file1 (@files) {
    my $file2 = $file1;
    $file2 =~ s/_1_/_2_/;
    my $lib = $1 if  basename($file1) =~ m/(\S+)_1_final.txt/;
    my $resultDir = $dir.$lib;
    mkdir $resultDir unless -d $resultDir;
    
    system("Trinity.pl --seqType fq -JM $opts{c}G --left $file1 --group_pairs_distance 999 --right $file2 --min_kmer_cov $opts{d} --CPU $opts{e} --output $resultDir");
    
    #system ("mv $resultDir'Trinity.fasta'  $resultDir$lib'.fasta' ");
    
    
  }  
}


sub abyss {
  
  die(qq/
Usage 3-GenerateAssemblies abyss [Options]

external dependencies: 
ABySS (compiled with OpenMPI and Google sparsehash)


Options: 
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-reads    DIR             Directory with all sequence reads
-lib      CHAR ...        Particular libraries to process? 
                          (e.g. AAA BBB CCC). If -lib is not 
                          used then process all libraries in
                          The folder (-reads) 
-mpi      DIR             Full path to openmpi mpirun   
-out      CHAR            Directory where results will go
-kmer     INT ...         A list of kmers. At least privide 
                          one kmer (e.g. 21 43 67 81)
-kcov     INT ...         A list of non-zero kmer-coverage. 
                          Activate default without using kcov
                          (e.g. 5 10 20). Suggest use default.
-np       INT             number of processors used for assembly
         
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


ABySS Multi-kmer assembly for exon capture data. 

Reference:
Singhal S. 2012. Molecular Ecology Resources 
Bi et al. 2012. BMC genomics
  
\n/) if (!@ARGV);
  
 
  my ($reads, $out,$mpi,$kmer, $np, $kcov, $lib) = (undef, undef, undef, undef, undef,undef,undef);
  GetOptions('kmer=s@{1,}' => \$kmer,'reads=s@{1,1}' => \$reads,'out=s@{1,1}' => \$out,'mpi=s@{1,1}' => \$mpi,'np=s@{1,1}' => \$np,'kcov=s@{,}' => \$kcov, 'lib=s@{,}' => \$lib );
   
  
  my $dir;
  
  if (@{$reads}[0] =~ m/\/$/ ){
    $dir = @{$reads}[0]; 
    }
  else {
    $dir = @{$reads}[0] . "/";
  }
  
  my $resDir;
  
   if (@{$out}[0] =~ m/\/$/ ){
    $resDir = @{$out}[0]; 
    }
  else {
    $resDir = @{$out}[0] . "/";
  }
  mkdir $resDir unless -d $resDir;
 
  my @files = ();
  if (!$lib) {
    @files= <$dir*_1_final.txt>;
  }
  
  if ($lib) { 
    foreach (@{$lib}) {
      my $file = $dir . $_ .'_1_final.txt';	
      push (@files, $file);      
    }
  } 
  
  
  foreach my $file1 (@files) {
    my $file2 = $file1;
    my $fileu = $file1;
    $file2 =~ s/_1_/_2_/;
    $fileu =~ s/_1_/_u_/;
    my $lib = $1 if  basename($file1) =~ m/(\S+)_1_final.txt/;
    my $resultDir = $resDir.$lib. "/";
    mkdir $resultDir unless -d $resultDir;
    foreach my $k (@{$kmer}) {
      if ($kcov){
	foreach my $cov (@{$kcov}) {
	  system("abyss-pe mpirun=@{$mpi}[0] k=$k c=$cov e=$cov np=@{$np}[0] E=0 n=10 in=\'$file1 $file2\' se=$fileu name=$resultDir$lib\"_k\"$k\"_cov\"$cov  ");
	  system("rm coverage.hist");
	}
      }
      if (!$kcov) {
	  system("abyss-pe mpirun=@{$mpi}[0] k=$k np=@{$np}[0] E=0 n=10 in=\'$file1 $file2\' se=$fileu name=$resultDir$lib\"_k\"$k\"_cov_default\"  ");
	  system("rm coverage.hist");
      }
    }
  } 
}
