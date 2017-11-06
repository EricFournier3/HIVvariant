GP120VariantComputer

Developped at the Laboratoire de santé publique du Québec, Canada, Qc

Current contributor:

	Eric Fournier


1. About

This script use as input deep sequencing reads from the HIV gp120 envelpe generated from chronic/recent patient isolates and their corresponding consensus dna sequence. For each specimen, reads are mapped to their respective consensus gp120 sequence in order to measure the level of
heterogeneity in different pre-established region. Final output is a statistic file compatible with R. Comparison between regions and chronic versus recent isolates can be made based on nucleotidic diversity, complexity and normalized Shannon entropy.

2. Running GP120VariantComputer

2.1 Dependencies
        -Linux operating system
	-Python 2.7
        -The following Python modules; Bio, yaml, argparse, textwrap

2.2 Command line options

	*To display help run 
         	python GP120VariantComputer.py --help

        *Required options are
		--basedir <String>
                 	Required. Path in which the working directory VariantAnalysis will be created

		--global-align-in <String>
			Required. Path to the global alignment with HXB2 reference in fasta format

		--yaml-spec-in <String>
			Required. Path to specimens list file in yaml format. See ExampleSpecYaml.yaml
			for example.

		--smalt-bin <String>
			Required. Path to Smalt executable

		--fastq-1 <String>
			Required. Path to fastq files

		--fastq-2 <String>
			Optional. Path to extra fastq files

2.3 Input

	User must specify path to the following input files :

			- A global dna gp120 alignement of specimens sequence with HXG2 as reference in fasta format.

			- A specimen list in yaml format. See ExampleSpecYaml.yaml as example.

			- Paired end fastq files (uncompressed or gz compressed). Two different path are accepted.  

2.4 Output

	The following directory scaffold is automatically created when the script is launched.

		VariantAnalysis
			|
			|_ Output        => The final statistic file StatResults.txt in generated in this directoy.
			|_ ErrorLog      => Log files is saved in this directory when exception is generated during diversity computation.
			|_ SmaltOutput   => Directory for mapping files generated by Smalt.
			|_ TempDiversity => Intermediary sequences alignement in fasta alignment are stored here. They are used for
                        |                   for diversity calculation by ComputeDiversity.R
			|_ Gp120Fasta    => Used to store g120 sequences from each specimen in fasta format. Those sequences are extracted
                                            from the global alignement input.

		
|-| | \/

	|-| | \/       |-| | \/

                |-| | \/        
