# coding=utf-8
__author__ = 'ericfournier'


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import yaml
import os
import re
import math
import subprocess
import argparse
import textwrap
import threading
import time

#Option parsing
parser = argparse.ArgumentParser(prog="GP120VariantComputer",formatter_class=argparse.RawDescriptionHelpFormatter,description=textwrap.dedent('''\
               |-| | \/
                     |-| | \/
                            |-| | \/
                            |-| | \/
                     |-| | \/
               |-| | \/
            '''))


#TODO ajouter / dans les path des base dir , confirmer le besoin dans le code

parser.add_argument('--basedir',metavar='[Required : Path to base directory]',required=True,nargs=1)

parser.add_argument('--global-align-in',type=str,metavar='[Required : Path to global alignement]',required=True,nargs=1,dest='globalalignin')

parser.add_argument('--yaml-spec-in',type=str,metavar='[Required : Path to yaml with specimens id]',required=True,nargs=1,dest='yamlspecin')

parser.add_argument('--smalt-bin',nargs=1,type=str,metavar='[Required : Path to Smalt executable]',required=True,dest='smaltbin')

parser.add_argument('--fastq-1',nargs=1,type=str,metavar='Required : Base directory for first set of fastq]',required=True,dest='fastq1')

#Optional
parser.add_argument('--fastq-2',nargs=1,type=str,metavar='[Optional : Base directory for second optional set of fastq]',required=False,dest='fastq2')


args=parser.parse_args()


#TODO SUPPRIMER LES PATH COMMENTé mais faire un check in avant

class PathManager():
    """
    Manage all path in this script
    """

    def __init__(self):

        #Path to base directory
        self.basedir = args.basedir[0]
        print 'basedir ' + self.basedir

        #Path to the global alignment file
        self.global_align_in = args.globalalignin[0]

        #Path to the yaml file with specimens list
        self.yaml_spec_in = args.yamlspecin[0]

        #Path to the smalt executable
        self.smalt_bin = args.smaltbin[0]

        #Path to fastq files
        if not args.fastq1[0].endswith('/'):
            self.fastq_1 = args.fastq1[0] + '/'
        else:
            self.fastq_1 = args.fastq1[0]

        print 'fastq_1 ' + self.fastq_1

        #Path second set set of fastq files if applicable
        try:
            if not args.fastq2[0].endswith('/'):
               self.fastq_2 = args.fastq2[0] + '/'
            else:
                self.fastq_2 = args.fastq2[0]

        except:
            self.fastq_2 = None
        print 'fastq_2 ' + str(self.fastq_2)


    def BuildWorkDir(self):
        """
        Build the working directory with all needed sub directories
        and files
        :return:
        """

        try:
            if not self.basedir.endswith('/'):
                self.basedir = self.basedir + '/'

            #Path to working directory
            self.work_dir = self.basedir + 'VariantAnalysis/'

            #Path to input directory
            self.input_dir = self.work_dir + 'Input/'

            #Path to output directory
            self.output_dir = self.work_dir + 'Output/'

            #Path to error directory files
            self.error_dir = self.work_dir + 'ErrorLog/'

            #Path to smalt output
            self.smalt_out_dir = self.work_dir + 'SmaltOutput/'

            #Path to fasta alignement pour diversity calculation
            self.diversity_align_dir = self.work_dir + 'TempDiversity/'

            #Path to gp120 assemblies fasta files
            self.gp120_seq_dir = self.work_dir + 'Gp120Fasta/'

            os.mkdir(self.work_dir)
            os.mkdir(self.input_dir)
            os.mkdir(self.output_dir)
            os.mkdir(self.error_dir)
            os.mkdir(self.smalt_out_dir)
            os.mkdir(self.diversity_align_dir)
            os.mkdir(self.gp120_seq_dir)

            print "Working directory created"

        except:
            print "Exception in build working directory"

#Create the working directory
my_path_manager = PathManager()
my_path_manager.BuildWorkDir()



class StatComputer():
    """
    Class used to compute mapping statistic.
    """

    def __init__(self,map_dir,ref_dir,dna_regions,region_length_inref):

        #Length of each region in HXB2
        self.region_length_inref=region_length_inref

        #Specimen list
        self.spec_list=[]

        #Base directory for smalt output
        self.map_dir=map_dir

        #Base directory of consensus sequence for each specimen
        self.ref_dir=ref_dir

        #Coordinate of each region in HXB2 in global alignement
        self.dna_regions=dna_regions

        '''
        For each specimen, coordinate of each region in mapping file
        Key => specimen id
        Value => another dictionary  Key => region
                                     Value => [start position, end position]
                                     
        '''
        self.dna_regions_inmap_file={}

        '''
        In each region, mapping statistic for each specimen
        Key => region
        Value => another dictionnary  Key => specimen id
                                      Value => this list [Complexity,Shannon,Diversity,Total number of sequences, Number of unique sequences]
                 
        '''
        self.all_stat={}

        '''
        For each specimen, number of gap in each region in the global alignment
        Key => specimen id
        Value => another dictionnary Key => Region
                                     Value => Number of gap
                                     
        see function CreateGapCountDict(self,spec_rec)
        
        '''
        self.dna_regions_gap={}

        #Number of specimens for each region
        self.nb_spec_per_region={'v1':0,'v2':0,'c2_1':0,'c2_2':0,'c2_3':0,'v3':0,'c3_1':0,'c3_2':0,'v4':0,'c4':0,'v5':0,'c5':0,'NHR':0}

    def InitAllStatDict(self):
        """
        Initialization of self.all_stat dictionnary

        Key => region
        Value => another dictionnary  Key => specimen
                                      Value => list of final stats
        :return:
        """

        for region in  self.dna_regions.keys():
            self.all_stat[region]={}

            for spec in self.spec_list:
                self.all_stat[region][spec]=[]

    def CreateGapCountDict(self,spec_rec):
        """
        For each specimen, compute for each region the number
        of gap in the golbal alignement
        :param spec_rec:
        :return:
        """

        '''
        A dictionnary for this specimen
        Key => Region
        Value => Number of gap in this region
        '''
        self.dna_regions_gap[spec_rec.id]={}

        #Total number of gap for this specimen in the global alignement
        tot_gap=0

        for region in self.dna_regions:
            self.dna_regions_gap[spec_rec.id][region]=self.CountGap(self.dna_regions[region],spec_rec.seq)
            tot_gap+=self.dna_regions_gap[spec_rec.id][region]

    def CountGap(self,location,dna_seq):
        """
        For the specimen, compute the number
        of gap in  the region

        :param location:
        :param dna_seq:
        :return:
        """
        return dna_seq[location[0]-1:location[1]].count('-')

    def CreateRegionsInMapFile(self,dna_region_list,spec):
        """
        Find coordinate for each region in a mapping file,

        :param dna_region_list:
        :param spec:
        :return:
        """

        #For each specimen, coordinate of each region in the mapping file
        self.dna_regions_inmap_file[spec]={}

        #Init the start position
        start_pos=1

        for region in dna_region_list:

            #End position of the region in mapping file
            stop_pos=self.FindStopPosInMapFile(self.dna_regions[region],self.dna_regions_gap[spec][region],start_pos)

            #Coordinate of the region in mapping file
            self.dna_regions_inmap_file[spec][region]=[start_pos,stop_pos]

            #Increment start position for next region
            start_pos=stop_pos+1

    def FindStopPosInMapFile(self,align_location,gap_nb,start_pos):
        """
        For one specimen, find the end position of the region
        in the mapping file
        :param align_location:
        :param gap_nb:
        :param start_pos:
        :return:
        """

        #Length of the region in the mapping file
        region_length=(align_location[1]-align_location[0]+1)-gap_nb
        #End position of the region
        stop_pos=region_length+start_pos-1

        return  stop_pos

    def AnalyseReadsInRegions(self,region):
        """
        For a specific specimen, extract unique sequences list in this region
        :param region:
        :return:
        """

        #a position of first match with read in the subject (reference sequence) in sam file
        start_pos_readOnRef=None

        #a position of last match with read in the subject (reference sequence) in sam file
        stop_pos_readOnRef = None

        #a mapping quality score in sam file
        qual=None

        #a cigar string in sam file
        cigar=None

        #a read sequence without soft and hard clipped nucleotide
        read_seq=None

        #a read sequence inside a region
        read_seq_onregion=None

        '''
        Boolean to continue the analysis for a specific specimen in this region. See next
        '''
        proceed=None

        #Line number in mapping file
        line_nb=0

        #For every specimen
        for spec in self.spec_list:

            proceed=True

            '''
            Dictionnary a unique read sequences covering this region for this specimen.
            Key => the unique sequence
            Value => absolute frequency
            '''
            uniq_seq_list={}

            #Length of the region in the mapping file
            spec_region_length=self.dna_regions_inmap_file[spec][region][1] - self.dna_regions_inmap_file[spec][region][0] + 1

            '''
            Keep extremity regions in the analysis only if their length are
            at least 80% of the corresponding region in HXB2 
            '''
            if region in ['v1','v2','NHR']:
                if float(spec_region_length)/self.region_length_inref[region] < 0.8:
                    proceed=False

            if proceed == True:

                '''
                Increment the number of specimens
                for which mapping statistic will
                be computed in this region
                '''
                self.nb_spec_per_region[region]+=1

                #Open map file and skip first three lines
                map_file=open(self.map_dir+spec+'.sam')
                map_file.readline()
                map_file.readline()
                map_file.readline()

                for line in map_file:
                    #line parsing
                    format_line=re.search(r'^\S+\t\S+\t\S+\t(\S+)\t(\S+)\t(\S+)\t\S+\t\S+\t\S+\t(\S+)',line)

                    #position of first match in subject (reference sequence) is in column 4
                    start_pos_readOnRef=int(format_line.group(1))

                    #mapping quality in column 5
                    qual=int(format_line.group(2))

                    #cigar string in column 6
                    cigar=format_line.group(3)

                    #read sequence without soft and hard clipping
                    read_seq=self.RemoveClippingSeq(format_line.group(4),cigar)

                    '''
                    Don't consider read alignment with low mapping quality 
                    '''
                    if cigar!='*' and qual >25:

                        #position of last match with read in subject
                        stop_pos_readOnRef=self.FindEndPos(cigar,start_pos_readOnRef)

                        '''
                        a read is kept only if his sequence cover completely the region
                        '''
                        if start_pos_readOnRef <= self.dna_regions_inmap_file[spec][region][0] and stop_pos_readOnRef >= self.dna_regions_inmap_file[spec][region][1]:

                            #Position of last match with read in subject
                            end_pos_readOnRef=self.FindEndPos(cigar,start_pos_readOnRef)

                            #read sequence inside the region
                            read_seq_onregion=self.FindReadSeqOnRegion(read_seq,start_pos_readOnRef,region,spec)

                            #add a new unique sequence
                            if read_seq_onregion not in uniq_seq_list:
                                uniq_seq_list[read_seq_onregion]=1
                            #increment its frequency
                            else:
                                uniq_seq_list[read_seq_onregion]+=1

                            line_nb+=1

                map_file.close()

            #Update statistic for this specimen
            self.UpdateStatDict(spec,uniq_seq_list,region)

    def SaveStatOutput(self):
        """
        Save Complexity, Shannon and Diversity for each
        specimen
        :return:
        """

        type=''

        '''
        Final mapping statistic file : 
        For each specimen, in each region we save Complexity, Shannon, Diversity, number of sequences 
        and number of haplotypes 
        '''
        #stat_file=open('/home/ericfournier/ProjetPersonnel/ProjetsNGS/HIV/SlidingWindow/Results/stat_results.txt','w')
        stat_file = open(my_path_manager.output_dir + 'StatResults.txt','w')
        # Header
        stat_file.write('Spec\tType\tRegion\tComplexity\tShannon\tDiversity\tNb_seq\tNb_haplotype\n')

        for region in self.all_stat:
            for spec in self.all_stat[region]:
                if str(spec).startswith('CH'):
                    #CHRONIC specimens
                    type='CHRONIC'
                else:
                    #RECENT specimens
                    type='RECENT'

                stat_file.write(spec+'\t'+type+'\t'+region+'\t')

                #Save all statistic for the current specimen in the current region
                for stat in self.all_stat[region][spec][0:-1]:
                    stat_file.write(str(stat)+'\t')
                stat_file.write(str(self.all_stat[region][spec][-1])+'\n')

        stat_file.close()

        '''
        Save the number of specimens contributing to each region
        '''
        #nbspec_per_region_file=open('/home/ericfournier/ProjetPersonnel/ProjetsNGS/HIV/SlidingWindow/Results/nb_spec_per_region.txt','w')
        nbspec_per_region_file = open(my_path_manager.output_dir + 'NbSpecPerRegion.txt','w')

        nbspec_per_region_file.write('Regions\tNb_specimen\n')

        for region in self.nb_spec_per_region:
            nbspec_per_region_file.write(region+'\t'+str(self.nb_spec_per_region[region])+'\n')
        nbspec_per_region_file.close()

    def UpdateStatDict(self,spec,seq_list,region):
        """
        Add Complexity, Shannon and Diversity values for this specimen
        in self.all_stat

        :param spec:
        :param seq_list:
        :param region:
        :return:
        """

        #Delete under-represented sequences in the region
        seq_list=self.FilterSeqList(seq_list)

        #List of kept sequences
        key_to_keep=seq_list.keys()

        #Total number of kept sequence
        seq_nb=sum(seq_list.values())

        #Satistic to compute in this region
        complexity=None
        shannon=None
        diversity=None

        #Only one unique sequence in this region
        if len(seq_list)==1:
            complexity=0
            shannon=0
            diversity=0
        #No sequence
        elif len(seq_list)==0:
            complexity='NA'
            shannon='NA'
            diversity='NA'
        else:
            #Compute statistics
            complexity=self.ComputeComplexity(seq_list)
            shannon=self.ComputeShannon(seq_list)
            diversity=self.ComputeDiversity(seq_list,spec,region)

        self.all_stat[region][spec].extend([complexity,shannon,diversity,sum(seq_list.values()),len(seq_list)])

    def ComputeDiversity(self,seq_list,spec,region):
        """
        Compute diversity from a multi fasta file built
        from the unique sequences list;
        Need R package seqinr, msa and pegas

        :param seq_list:
        :param spec:
        :param region:
        :return:
        """

        #Create a multi fasta file from unique sequences list
        self.CreateFastaFromUniqReads(seq_list)

        diversity=0

        #Path to the multi fasta file
        #fasta_read="/home/ericfournier/ProjetPersonnel/ProjetsNGS/HIV/SlidingWindow/ReadAlign/temp_read.fasta"
        fasta_read = my_path_manager.diversity_align_dir + 'temp_read.fasta'

        #Compute diversity with the in-house R script ComputeDiversity.R
        #r_output=subprocess.check_output("Rscript /home/ericfournier/ProjetPersonnel/ProjetsNGS/HIV/SlidingWindow/RSCRIPT/ComputeDiversity.R "+fasta_read, shell=True)
        r_output = subprocess.check_output("Rscript " + "ComputeDiversity.R" + " " + fasta_read, shell=True)

        #Extract the diversity value from the R output
        try:
            diversity=re.search(r'\[1\] (\d+\.*\d*)',r_output).group(1)

        except:
            #with open("/home/ericfournier/ProjetPersonnel/ProjetsNGS/HIV/SlidingWindow/ERROR/error.txt",'a') as writef:
            with open(my_path_manager.error_dir+'DiversityErrorLog','a') as writef:
                writef.write('No diversity for '+spec+'\n')
            writef.close()

            #Keep a copy of the alignment if exception
            #os.system('cp /home/ericfournier/ProjetPersonnel/ProjetsNGS/HIV/SlidingWindow/ReadAlign/temp_read.fasta /home/ericfournier/ProjetPersonnel/ProjetsNGS/HIV/SlidingWindow/ReadAlign/'+'temp_read_'+spec+'_'+region+'.fasta')
            os.system('cp ' + fasta_read + " " + my_path_manager.error_dir + 'temp_read_' + spec + '_' + region + '.fasta')
            diversity='NA'

        os.system('rm '+fasta_read)

        return diversity

    def CreateFastaFromUniqReads(self,seq_list):
        """
        Create a multi-records fasta file from unique sequences list

        :param seq_list:
        :return:
        """

        #List of unique sequences in fasta format
        rec_list=[]

        #Build rec_list
        for key in seq_list:
            rec_list.append(SeqRecord(Seq(key),id=str(seq_list[key]),description='NA'))

        #Create the corresponding multi fasta file
        #SeqIO.write(rec_list,"/home/ericfournier/ProjetPersonnel/ProjetsNGS/HIV/SlidingWindow/ReadAlign/temp_read.fasta","fasta")
        SeqIO.write(rec_list,my_path_manager.diversity_align_dir + 'temp_read.fasta', "fasta")

    def ComputeShannon(self,seq_list):
        """
        Normalized Shannon entropy

        = -[1/log(Total number of sequences)] * SUM[i=1 to Number of unique sequence]Prob(i)log(Prob(i))

        :param seq_list:
        :return:
        """

        # Number of unique sequences
        nb_uniq_seq=len(seq_list)

        # Total number of sequences
        nb_total_seq=sum(seq_list.values())

        #Compute SUM[i=1 to Number of unique sequence]Prob(i)log(Prob(i))
        shannon_sum=0
        for key in seq_list:
            prob=float(seq_list[key])/nb_total_seq
            shannon_sum+=prob*math.log(prob,10)

        #Compute normalized Shannon entropy
        shannon=round(-(1/math.log(nb_total_seq,10))*shannon_sum,10)

        return shannon

    def ComputeComplexity(self,seq_list):
        """
        Comlexity is simply the number of unique sequences
        divided by the total number of sequences

        :param seq_list:
        :return:
        """

        #Number of unique sequences
        nb_uniq_seq=len(seq_list)

        #Total number of sequences
        nb_total_seq=sum(seq_list.values())

        try:
            complexity=round(float(nb_uniq_seq)/nb_total_seq,10)

        except:
            complexity='NA'

        return  complexity

    def FilterSeqList(self,seq_list):
        """
        Clean the unique sequences list in a specific region for a given specimen.
        A unique sequence is deleted if its abundance is not at least 0.5%
        :param seq_list:
        :return:
        """

        '''
        Total number of reads kept in this region
        '''
        seq_nb=sum(seq_list.values())

        #Unique sequences list to delete
        key_to_del=[]

        '''
        Delete unique sequences that do not represent at leat 0.5% of total sequences
        '''
        for key in seq_list:
            if  not ((float(seq_list[key])/seq_nb)*100)  >= 0.5:
                key_to_del.append(key)

        for key in key_to_del:
            del(seq_list[key])

        return seq_list

    def RemoveClippingSeq(self,read_seq,cigar):
        """
        Remove hard and soft clipping nucleotide from the read
        sequence in sam file
        :param read_seq:
        :param cigar:
        :return:
        """

        #number of nucleotides to remove in 3'
        stop_clip=0

        #number of nucleotides to remove in 5'
        start_clip=0

        try:
            #soft clipping present in 5' and 3'
            start_clip=int(re.search(r'^(\d+)S',cigar).group(1))
            stop_clip=int(re.search(r'(\d+)S$',cigar).group(1))

        except:
            try:
                #soft clipping only in 5'
                start_clip=int(re.search(r'^(\d+)S',cigar).group(1))

            except:
                try:
                    #soft clipping only in 3'
                    stop_clip=int(re.search(r'(\d+)S$',cigar).group(1))

                except:
                    pass

        '''
        Note: there is no hard clipping in the smalt
        output files 
        '''

        '''
        Slice the read accoringly
        '''
        if stop_clip==0:
            read_seq_clip=read_seq[start_clip:]
        else:
            read_seq_clip=read_seq[start_clip:-stop_clip]

        return read_seq_clip

    def FindReadSeqOnRegion(self,read_seq,start_pos,region,spec):
        """
        Extract the nucleotide sequence of the read
        that is inside the target region
        :param read_seq:
        :param start_pos:
        :param region:
        :param spec:
        :return:
        """

        #start position of the region in map file
        start_reg=self.dna_regions_inmap_file[spec][region][0]

        #end position of the region in map file
        stop_reg=self.dna_regions_inmap_file[spec][region][1]

        #convert read sequence in nucleotid list
        nuc_list=list(read_seq)

        '''
        Position of each nucleotide of the read in the reference
        '''
        read_pos_onref_list=[]

        '''
        Nucleotide sequence of the read inside the target region
        '''
        read_seq_inregion=''

        for cigar_pair in self.cigar_list:

            #simply add the position value in subject if its a match and increment position value
            if cigar_pair[1] in ['M']:
                i=0
                while i< cigar_pair[0]:
                    read_pos_onref_list.append(start_pos)
                    start_pos+=1
                    i+=1

            #if insertion => add the previous position value
            elif cigar_pair[1] in ['I']:
                read_pos_onref_list.append(read_pos_onref_list[-1])

            #if deletion => only increment the position value
            elif cigar_pair[1] in ['D']:
                start_pos+=1
                pass

        #List of tuple : match read nuclotide to its corresponding position in the subject
        nuc_pos_onref_list=zip(nuc_list,read_pos_onref_list)

        #Build the read sequence covering the region. i.e cut in 5' and 3'
        for nuc_pos_pair in nuc_pos_onref_list:

            if nuc_pos_pair[1]>= start_reg and nuc_pos_pair[1] <= stop_reg:
                read_seq_inregion+=nuc_pos_pair[0]
                #print read_seq_inregion

        return  read_seq_inregion

    def FindEndPos(self,cigar,start_pos):
        """
        Find position of last match in subject (reference sequence)
        for alignment with one read in sam file
        :param cigar:
        :param start_pos:
        :return:
        """

        #parsed cigar in tuple list format
        self.cigar_list=self.ParseCigar(cigar)

        stop_pos=start_pos

        '''
        Recall:
        M => MATCH : Can be difference or identity
        D => DELETION : Nucleotide over in reference
        I => INSERTION : Nucleotide over in read
        '''
        for my_pair in self.cigar_list:
            if my_pair[1] in ['M','D']:

                #Increment stop_pos from start_pos only for match and deletion
                stop_pos+=my_pair[0]

        stop_pos+=-1

        return stop_pos


    def ParseCigar(self,cigar):
        """
        Parse a cigar string and output
        a list a tuple having the following format : [(digit,char)]
        Example: cigar 2S80M94S will produce [(2, 'S'), (80, 'M'), (94, 'S')]
        :param cigar:
        :return:
        """

        digit=''

        digit_list=[]
        char_list=[]

        for i in cigar:
            #Build the numerical value
            if i.isdigit():
                digit+=i
            #Build the numerical and the character list
            else:
                digit_list.append(int(digit))
                char_list.append(i)

                digit=''

        #Build the tuple list
        self.cigar_list=zip(digit_list,char_list)

        return self.cigar_list

class Mapper():
    """
    Smalt mapping
    """

    def __init__(self,out_dir):

        #Output directory
        self.out_dir=out_dir

        #Path to smalt executable
        #self.smalt_prog='/home/ericfournier/ProgramInternet/SMALT/smalt-0.7.4/smalt_i686'
        self.smalt_prog = my_path_manager.smalt_bin

        #Path to paired end fastq; first and second set

        #self.fastq_hiv1='/media/ericfournier/Elements/MiSeq/20160316_AlexVIH3/{0}/FASTQ/'
        self.fastq_hiv1 = my_path_manager.fastq_1 + '{0}/FASTQ/'

        #self.fastq_hiv2='/media/ericfournier/Elements/MiSeq/20160321_AlexVIH4/{0}/FASTQ/'
        if my_path_manager.fastq_2 != None:
            self.fastq_hiv2 = my_path_manager.fastq_2 + '{0}/FASTQ/'
        else:
            pass

        #For one specimen fastq paired end path
        self.fastq_pair=[]

        #Smalt parameters
        self.param_k=' -k 20 '
        self.param_s=' -s 13 '
        self.param_f=' -f bam '
        self.param_l=' -l pe '
        self.param_i=' -i 500 '
        self.param_o_part=' -o '+self.out_dir

    def IndexRef(self,fasta_file):
        """
        Index a specimen fasta sequence
        :param fasta_file:
        :return:
        """

        #Prefix for the output file
        index_file_prefix=re.sub(r'\.fasta','',fasta_file)

        #Execute indexing
        os.system(self.smalt_prog +' index '+self.param_k+self.param_s+index_file_prefix+' '+fasta_file)

    def FindFastqPair(self,spec):
        """
        Make the paired end fastq path for this
        specimen

        :param spec:
        :return:
        """

        #For one fastq paired end path
        self.fastq_pair=[]

        #Find the fastq paired end path
        try:
            for fastq in os.listdir(self.fastq_hiv1.format(spec)):
                self.fastq_pair.append(self.fastq_hiv1.format(spec)+'/'+fastq)
        except:
            print "IN EXCEPT"
            for fastq in os.listdir(self.fastq_hiv2.format(spec)):
                self.fastq_pair.append(self.fastq_hiv2.format(spec)+'/'+fastq)

    def Map(self,spec,index_file):
        """
        Execution of the smalt mapping for one specimen
        :param spec:
        :param index_file:
        :return:
        """
        '''
        For test only
        os.system('{0}{1}{2}{3}{4}{5}{6}{7} {8} {9}'.format(self.smalt_prog,' map ',self.param_f,self.param_l,self.param_i,self.param_o_part,spec+'.bam ',index_file,self.fastq_pair[0],self.fastq_pair[1]))
        os.system('samtools sort -f {0}{1}{2} {0}{1}_sort{2}'.format(self.out_dir,spec,'.bam'))
        os.system('samtools calmd -e {0}{1}{2} {3}.fasta > {0}{1}.sam'.format(self.out_dir,spec,'_sort.bam ',index_file))
        os.system('rm {0}{1} {0}{2}'.format(self.out_dir,spec+'.bam',spec+'_sort.bam'))
        #os.system("samtools calmd -e {0}Smalt_{2}_sort.bam {0}{2}_cut.fasta > {1}{2}_sort_md.sam".format(self.in_dir,self.out_bam_dir,spec))
        '''

        #map => bam file
        os.system('{0}{1}{2}{3}{4}{5}{6}{7} {8} {9}'.format(self.smalt_prog,' map ',self.param_f,self.param_l,self.param_i,self.param_o_part,spec+'.bam ',index_file,self.fastq_pair[0],self.fastq_pair[1]))

        #sort the bam file
        os.system('samtools sort -f {0}{1}{2} {0}{1}_sort{2}'.format(self.out_dir,spec,'.bam'))
        #convert to sam
        os.system('samtools view -h {0}{1}{3} > {0}{1}{2}'.format(self.out_dir,spec,'.sam','_sort.bam'))

        #os.system('samtools calmd -e {0}{1}{2} {3}.fasta > {0}{1}.sam'.format(self.out_dir,spec,'_sort.bam ',index_file))
        #remove bam
        os.system('rm {0}{1} {0}{2}'.format(self.out_dir,spec+'.bam',spec+'_sort.bam'))

        #for test only
        #os.system("samtools calmd -e {0}Smalt_{2}_sort.bam {0}{2}_cut.fasta > {1}{2}_sort_md.sam".format(self.in_dir,self.out_bam_dir,spec))

class WindowAnalyser():

    def __init__(self):

        #flag for program status
        self.finish = False

        print "This process has the pid", os.getpid()

        #Display wait message during program execution
        try:
            t = threading.Thread(target=self.WaitMessage)
            t.start()
        except:
            print "Error: unable to start thread"
            #End program
            exit(0)

        #Specimen list
        self.spec_list=[]

        #Ordered region of gp120
        self.gp120_region_ordered_list=['v1','v2','c2_1','c2_2','c2_3','v3','c3_1','c3_2','v4','c4','v5','c5','NHR']

        #Length of each region in reference HXB2
        self.gp120_region_length_inref={}

        #Coordinate of each region in HXB2 in global alignement
        self.gp120_region={'v1':(415,495),'v2':(496,642),'c2_1':(643,742),'c2_2':(743,842),'c2_3':(843,954),'v3':(955,1071),'c3_1':(1072,1173),
        'c3_2':(1174,1276),'v4':(1277,1388),'c4':(1389,1481),'v5':(1482,1547),'c5':(1548,1700),'NHR':(1701,1860)}

        #Base directory of consensus sequence for each specimen
        #self.base_dir_ref="/home/ericfournier/ProjetPersonnel/ProjetsNGS/HIV/SlidingWindow/TempReference/"
        self.base_dir_ref = my_path_manager.gp120_seq_dir

        #Base directory for smalt output
        #self.base_dir_map="/media/ericfournier/Elements/MiSeq/HIV_slidingwindow_map/"
        self.base_dir_map = my_path_manager.smalt_out_dir

        #Specimen list in yaml file
        #self.yam_file=open('/home/ericfournier/ProjetPersonnel/ProjetsNGS/HIV/SlidingWindow/HIV.yaml').read()
        self.yam_file = open(my_path_manager.yaml_spec_in).read()

        #Number of specimen
        self.NbSpec=yaml.load(self.yam_file)[0]

        #Specimen list
        self.SpecList=yaml.load(self.yam_file)[1:self.NbSpec+1]

        #Test only
        #self.global_align_iter=SeqIO.parse(self.base_dir_seq+'GlobalHIVAlign_test.fas','fasta')

        #List of fasta record from global alignement
        #self.global_align_iter=SeqIO.parse(self.base_dir_seq+'GlobalHIVAlign_original.fas','fasta')
        self.global_align_iter = SeqIO.parse(my_path_manager.global_align_in,'fasta')

        #HXB2 record in fasta format
        ref_rec=self.global_align_iter.next()

        self.ComputeRegionLength(ref_rec)

        self.mapper=Mapper(self.base_dir_map)

        self.stat_comput=StatComputer(self.base_dir_map,self.base_dir_ref,self.gp120_region,self.gp120_region_length_inref)


        print " *********  Start Analysis **********"
        self.StartAnalysis()



        print " *********  Save Statistics **********"
        self.SaveStat()

        #program has ended
        self.finish = True

        return

    def WaitMessage(self):
        """
        Printout during program execution
        :param delay:
        :return:
        """
        while not self.finish:
            time.sleep(30)
            print "%s: %s" % ("Please wait", time.ctime(time.time()))


    def ComputeRegionLength(self,ref_rec):
        """
        Compute length of each region without gap in HXB2.
        ref_rec is the sequence of HXB2 in global alignement.
        :param ref_rec:
        :return:
        """
        for region in self.gp120_region:

            #Length of the region with gap
            length_with_gap=len(ref_rec.seq[self.gp120_region[region][0]-1:self.gp120_region[region][1]])

            #Number of gap in the region
            nb_gap=ref_rec.seq[self.gp120_region[region][0]-1:self.gp120_region[region][1]].count('-')

            #Length of the region without gap
            self.gp120_region_length_inref[region]=length_with_gap-nb_gap


    def StartAnalysis(self):
        """
        Initialize the analysis

        :return:
        """

        #For each specimen record in the global alignment
        for my_rec in self.global_align_iter:

            #Add specimen id to specimens list
            self.stat_comput.spec_list.append(my_rec.id)

            self.CreateOneSpecFasta(my_rec)
            self.mapper.IndexRef(self.base_dir_ref+my_rec.id+'.fasta')
            self.mapper.FindFastqPair(my_rec.id)

            print "Mapping step for " + my_rec.id

            self.mapper.Map(my_rec.id,self.base_dir_ref+my_rec.id)

            self.stat_comput.CreateGapCountDict(my_rec)

            self.stat_comput.CreateRegionsInMapFile(self.gp120_region_ordered_list,my_rec.id)

        self.stat_comput.InitAllStatDict()

        '''
        Genetic heterogeneity evaluation in each region, for each
        specimen start here
        '''
        for region in self.gp120_region_ordered_list:
            self.stat_comput.AnalyseReadsInRegions(region)

    def SaveStat(self):
        """
        Save all statistics for each specimen in each region
        :return:
        """

        self.stat_comput.SaveStatOutput()

    def CreateOneSpecFasta(self,rec):
        """
        For a specific specimen, create a single fasta file without gap
        :param rec:
        :return:
        """

        #Extract sequence for this specimen record and delete gap
        rec.seq=Seq(re.sub('-','',str(rec.seq)))

        #Create and save this fasta record. It will be use later as reference for mapping
        SeqIO.write(rec,'{1}{0}.fasta'.format(rec.id,self.base_dir_ref),'fasta')


#Start the program
my_wa=WindowAnalyser()

'''
Command Example
python GP120VariantComputer.py --basedir '/home/eric/ProjetPersonnel/ProjetsNGS/HIV/SlidingWindow/TestGitCode' --global-align-in '/home/eric/ProjetPersonnel/ProjetsNGS/HIV/SlidingWindow/TestGitCode/Alignment/Alignement.fas' --yaml-spec-in '/home/eric/ProjetPersonnel/ProjetsNGS/HIV/SlidingWindow/TestGitCode/SpecList/ExampleSpecYaml.yaml' --smalt-bin '/home/eric/InternetProgram/Smalt/smalt-0.7.4/smalt_i686' --fastq-1 '/media/eric/Elements/MiSeq/20160316_AlexVIH3'


'''