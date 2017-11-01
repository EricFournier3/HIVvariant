__author__ = 'ericfournier'


from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import yaml
import os
import re
import math
import subprocess

class R_GraphMaker():
    def __init__(self):
        pass

class StatComputer():
    def __init__(self,map_dir,stat_dir,ref_dir,dna_regions,region_length_inref):

        self.region_length_inref=region_length_inref
        self.spec_list=[]
        self.map_dir=map_dir
        self.stat_dir=stat_dir
        self.ref_dir=ref_dir
        self.dna_regions=dna_regions
        self.dna_regions_inmap_file={}
        self.all_stat={}
        self.dna_regions_gap={}
        self.nb_spec_per_region={'v1':0,'v2':0,'c2_1':0,'c2_2':0,'c2_3':0,'v3':0,'c3_1':0,'c3_2':0,'v4':0,'c4':0,'v5':0,'c5':0,'NHR':0}

    def InitAllStatDict(self):

        for region in self.dna_regions.keys():
            self.all_stat[region]={}

            for spec in self.spec_list:
                self.all_stat[region][spec]=[]

    def CreateGapCountDict(self,spec_rec):

        self.dna_regions_gap[spec_rec.id]={}

        tot_gap=0

        for region in self.dna_regions:

            self.dna_regions_gap[spec_rec.id][region]=self.CountGap(self.dna_regions[region],spec_rec.seq)
            tot_gap+=self.dna_regions_gap[spec_rec.id][region]

    def CountGap(self,location,dna_seq):
        return dna_seq[location[0]-1:location[1]].count('-')

    def CreateRegionsInMapFile(self,dna_region_list,spec):

        self.dna_regions_inmap_file[spec]={}
        start_pos=1

        for region in dna_region_list:
            stop_pos=self.FindStopPosInMapFile(self.dna_regions[region],self.dna_regions_gap[spec][region],start_pos)
            self.dna_regions_inmap_file[spec][region]=[start_pos,stop_pos]
            start_pos=stop_pos+1

    def FindStopPosInMapFile(self,align_location,gap_nb,start_pos):

        region_length=(align_location[1]-align_location[0]+1)-gap_nb
        stop_pos=region_length+start_pos-1

        return  stop_pos

    def AnalyseReadsInRegions(self,region):

        start_pos_readOnRef=None
        qual=None
        cigar=None
        read_seq=None
        stop_pos_readOnRef=None
        read_seq_onregion=None
        proceed=None

        line_nb=0

        for spec in self.spec_list:

            proceed=True

            uniq_seq_list={}

            print ">>>>>> SPEC IS: ",spec

            spec_region_length=self.dna_regions_inmap_file[spec][region][1] - self.dna_regions_inmap_file[spec][region][0] + 1

            #if (self.dna_regions_inmap_file[spec][region][1] - self.dna_regions_inmap_file[spec][region][0]) > 20:
            #if spec_region_length >= self.region_length_inref[region]-10 and spec_region_length <= self.region_length_inref[region]+10 :
            #if float(spec_region_length)/self.region_length_inref[region] >=0.8:

            if region in ['v1','v2','NHR']:
                if float(spec_region_length)/self.region_length_inref[region] < 0.8:
                    proceed=False

            if proceed==True:

                #print spec_region_length

                self.nb_spec_per_region[region]+=1

                map_file=open(self.map_dir+spec+'.sam')
                map_file.readline()
                map_file.readline()
                map_file.readline()

                for line in map_file:
                    format_line=re.search(r'^\S+\t\S+\t\S+\t(\S+)\t(\S+)\t(\S+)\t\S+\t\S+\t\S+\t(\S+)',line)
                    start_pos_readOnRef=int(format_line.group(1))
                    qual=int(format_line.group(2))
                    cigar=format_line.group(3)
                    read_seq=self.RemoveClippingSeq(format_line.group(4),cigar)

                    if cigar!='*' and qual >25:
                        stop_pos_readOnRef=self.FindEndPos(cigar,start_pos_readOnRef)

                        if start_pos_readOnRef <= self.dna_regions_inmap_file[spec][region][0] and stop_pos_readOnRef >= self.dna_regions_inmap_file[spec][region][1]:

                            end_pos_readOnRef=self.FindEndPos(cigar,start_pos_readOnRef)

                            read_seq_onregion=self.FindReadSeqOnRegion(read_seq,start_pos_readOnRef,region,spec)

                            if read_seq_onregion not in uniq_seq_list:
                                uniq_seq_list[read_seq_onregion]=1
                            else:
                                uniq_seq_list[read_seq_onregion]+=1

                            line_nb+=1

                map_file.close()

            self.UpdateStatDict(spec,uniq_seq_list,region)

    def SaveStatOutput(self):

        type=''

        stat_file=open('/home/ericfournier/ProjetPersonnel/ProjetsNGS/HIV/SlidingWindow/Results/stat_results.txt','w')
        stat_file.write('Spec\tType\tRegion\tComplexity\tShannon\tDiversity\tNb_seq\tNb_haplotype\n')

        #for region in ['v3']:
        for region in self.all_stat:
            for spec in self.all_stat[region]:
                if str(spec).startswith('CH'):
                    type='CHRONIC'
                else:
                    type='RECENT'

                stat_file.write(spec+'\t'+type+'\t'+region+'\t')

                for stat in self.all_stat[region][spec][0:-1]:
                    stat_file.write(str(stat)+'\t')
                stat_file.write(str(self.all_stat[region][spec][-1])+'\n')

        stat_file.close()

        nbspec_per_region_file=open('/home/ericfournier/ProjetPersonnel/ProjetsNGS/HIV/SlidingWindow/Results/nb_spec_per_region.txt','w')

        nbspec_per_region_file.write('Regions\tNb_specimen\n')

        for region in self.nb_spec_per_region:
            nbspec_per_region_file.write(region+'\t'+str(self.nb_spec_per_region[region])+'\n')
        nbspec_per_region_file.close()

    def UpdateStatDict(self,spec,seq_list,region):

        seq_list=self.FilterSeqList(seq_list)

        key_to_keep=seq_list.keys()

        seq_nb=sum(seq_list.values())

        complexity=None
        shannon=None
        diversity=None

        if len(seq_list)==1:
            complexity=0
            shannon=0
            diversity=0
        elif len(seq_list)==0:
            complexity='NA'
            shannon='NA'
            diversity='NA'
        else:

            complexity=self.ComputeComplexity(seq_list)
            shannon=self.ComputeShannon(seq_list)
            diversity=self.ComputeDiversity(seq_list,spec,region)

        self.all_stat[region][spec].extend([complexity,shannon,diversity,sum(seq_list.values()),len(seq_list)])

    def ComputeDiversity(self,seq_list,spec,region):

        self.CreateFastaFromUniqReads(seq_list)

        diversity=0

        fasta_read="/home/ericfournier/ProjetPersonnel/ProjetsNGS/HIV/SlidingWindow/ReadAlign/temp_read.fasta"

        r_output=subprocess.check_output("Rscript /home/ericfournier/ProjetPersonnel/ProjetsNGS/HIV/SlidingWindow/RSCRIPT/ComputeDiversity.R "+fasta_read, shell=True)

        try:
            diversity=re.search(r'\[1\] (\d+\.*\d*)',r_output).group(1)

        except:
            with open("/home/ericfournier/ProjetPersonnel/ProjetsNGS/HIV/SlidingWindow/ERROR/error.txt",'a') as writef:
                writef.write('No diversity for '+spec+'\n')
            writef.close()

            os.system('cp /home/ericfournier/ProjetPersonnel/ProjetsNGS/HIV/SlidingWindow/ReadAlign/temp_read.fasta /home/ericfournier/ProjetPersonnel/ProjetsNGS/HIV/SlidingWindow/ReadAlign/'+'temp_read_'+spec+'_'+region+'.fasta')
            diversity='NA'

        os.system('rm '+fasta_read)

        return diversity

    def CreateFastaFromUniqReads(self,seq_list):

        rec_list=[]

        for key in seq_list:
            rec_list.append(SeqRecord(Seq(key),id=str(seq_list[key]),description='NA'))

        SeqIO.write(rec_list,"/home/ericfournier/ProjetPersonnel/ProjetsNGS/HIV/SlidingWindow/ReadAlign/temp_read.fasta","fasta")

    def ComputeShannon(self,seq_list):

        nb_uniq_seq=len(seq_list)

        nb_total_seq=sum(seq_list.values())

        shannon_sum=0

        for key in seq_list:
            prob=float(seq_list[key])/nb_total_seq
            shannon_sum+=prob*math.log(prob,10)

        shannon=round(-(1/math.log(nb_total_seq,10))*shannon_sum,10)

        return shannon

    def ComputeComplexity(self,seq_list):

        nb_uniq_seq=len(seq_list)
        nb_total_seq=sum(seq_list.values())

        try:
            complexity=round(float(nb_uniq_seq)/nb_total_seq,10)

        except:
            complexity='NA'

        return  complexity

    def FilterSeqList(self,seq_list):
        pass

        seq_nb=sum(seq_list.values())

        key_to_del=[]

        for key in seq_list:
            if  not ((float(seq_list[key])/seq_nb)*100)  >= 0.5:
                key_to_del.append(key)

        for key in key_to_del:
            del(seq_list[key])

        return seq_list

    def RemoveClippingSeq(self,read_seq,cigar):
        stop_clip=0
        star_clip=0

        try:
            star_clip=int(re.search(r'^(\d+)[SH]',cigar).group(1))
            stop_clip=int(re.search(r'(\d+)S$',cigar).group(1))

        except:
            try:
                star_clip=int(re.search(r'^(\d+)[SH]',cigar).group(1))

            except:
                try:
                    stop_clip=int(re.search(r'(\d+)[SH]$',cigar).group(1))

                except:
                    pass

        if stop_clip==0:
            read_seq_clip=read_seq[star_clip:]
        else:
            read_seq_clip=read_seq[star_clip:-stop_clip]

        return read_seq_clip

    def FindReadSeqOnRegion(self,read_seq,start_pos,region,spec):

        start_reg=self.dna_regions_inmap_file[spec][region][0]
        stop_reg=self.dna_regions_inmap_file[spec][region][1]

        nuc_list=list(read_seq)
        read_pos_onref_list=[]

        read_seq_inregion=''

        for cigar_pair in self.cigar_list:
            if cigar_pair[1] in ['M']:
                i=0
                while i< cigar_pair[0]:
                    #print 'i ', i, 'cigar_pair ',cigar_pair[0]
                    read_pos_onref_list.append(start_pos)
                    start_pos+=1
                    i+=1
                    #print "READ POS ", read_pos_onref_list
            elif cigar_pair[1] in ['I']:
                read_pos_onref_list.append(read_pos_onref_list[-1])

            elif cigar_pair[1] in ['D']:
                start_pos+=1
                pass
                #print "testing"

        nuc_pos_onref_list=zip(nuc_list,read_pos_onref_list)

        for nuc_pos_pair in nuc_pos_onref_list:

            if nuc_pos_pair[1]>= start_reg and nuc_pos_pair[1] <= stop_reg:
                read_seq_inregion+=nuc_pos_pair[0]
                #print read_seq_inregion

        return  read_seq_inregion

    def FindEndPos(self,cigar,start_pos):

        self.cigar_list=self.ParseCigar(cigar)

        stop_pos=start_pos

        for my_pair in self.cigar_list:
            if my_pair[1] in ['M','D']:
                stop_pos+=my_pair[0]
                #print 'PAIR ', my_pair[0]

        stop_pos+=-1

        return stop_pos


    def ParseCigar(self,cigar):

        digit=''
        char=''

        digit_list=[]
        char_list=[]

        for i in cigar:
            if i.isdigit():
                digit+=i
            else:
                digit_list.append(int(digit))
                char_list.append(i)

                digit=''

        self.cigar_list=zip(digit_list,char_list)

        return self.cigar_list

class Mapper():

    def __init__(self,out_dir):

        self.out_dir=out_dir

        self.smalt_prog='/home/ericfournier/ProgramInternet/SMALT/smalt-0.7.4/smalt_i686 '
        self.fastq_hiv3='/media/ericfournier/Elements/MiSeq/20160316_AlexVIH3/{0}/FASTQ/'
        self.fastq_hiv4='/media/ericfournier/Elements/MiSeq/20160321_AlexVIH4/{0}/FASTQ/'

        self.fastq_pair=[]

        self.param_k=' -k 20 '
        self.param_s=' -s 13 '
        self.param_f=' -f bam '
        self.param_l=' -l pe '
        self.param_i=' -i 500 '
        self.param_o_part=' -o '+self.out_dir

    def IndexRef(self,fasta_file):

        index_file_prefix=re.sub(r'\.fasta','',fasta_file)
        os.system(self.smalt_prog +'index'+self.param_k+self.param_s+index_file_prefix+' '+fasta_file)

    def FindFastqPair(self,spec):

        self.fastq_pair=[]

        try:
            for fastq in os.listdir(self.fastq_hiv3.format(spec)):
                #print fastq
                self.fastq_pair.append(self.fastq_hiv3.format(spec)+'/'+fastq)
        except:
            for fastq in os.listdir(self.fastq_hiv4.format(spec)):
                #print fastq

                self.fastq_pair.append(self.fastq_hiv4.format(spec)+'/'+fastq)

        print spec, " ",self.fastq_pair

    def Map(self,spec,index_file):
        '''
        os.system('{0}{1}{2}{3}{4}{5}{6}{7} {8} {9}'.format(self.smalt_prog,' map ',self.param_f,self.param_l,self.param_i,self.param_o_part,spec+'.bam ',index_file,self.fastq_pair[0],self.fastq_pair[1]))

        os.system('samtools sort -f {0}{1}{2} {0}{1}_sort{2}'.format(self.out_dir,spec,'.bam'))


        os.system('samtools calmd -e {0}{1}{2} {3}.fasta > {0}{1}.sam'.format(self.out_dir,spec,'_sort.bam ',index_file))
        os.system('rm {0}{1} {0}{2}'.format(self.out_dir,spec+'.bam',spec+'_sort.bam'))
        #os.system("samtools calmd -e {0}Smalt_{2}_sort.bam {0}{2}_cut.fasta > {1}{2}_sort_md.sam".format(self.in_dir,self.out_bam_dir,spec))
        '''

        os.system('{0}{1}{2}{3}{4}{5}{6}{7} {8} {9}'.format(self.smalt_prog,' map ',self.param_f,self.param_l,self.param_i,self.param_o_part,spec+'.bam ',index_file,self.fastq_pair[0],self.fastq_pair[1]))

        os.system('samtools sort -f {0}{1}{2} {0}{1}_sort{2}'.format(self.out_dir,spec,'.bam'))
        os.system('samtools view -h {0}{1}{3} > {0}{1}{2}'.format(self.out_dir,spec,'.sam','_sort.bam'))


        #os.system('samtools calmd -e {0}{1}{2} {3}.fasta > {0}{1}.sam'.format(self.out_dir,spec,'_sort.bam ',index_file))
        os.system('rm {0}{1} {0}{2}'.format(self.out_dir,spec+'.bam',spec+'_sort.bam'))
        #os.system("samtools calmd -e {0}Smalt_{2}_sort.bam {0}{2}_cut.fasta > {1}{2}_sort_md.sam".format(self.in_dir,self.out_bam_dir,spec))

class WindowAnalyser():

    def __init__(self):

        self.spec_list=[]

        self.gp120_region_ordered_list=['v1','v2','c2_1','c2_2','c2_3','v3','c3_1','c3_2','v4','c4','v5','c5','NHR']

        self.gp120_region_length_inref={}

        self.gp120_region={'v1':(415,495),'v2':(496,642),'c2_1':(643,742),'c2_2':(743,842),'c2_3':(843,954),'v3':(955,1071),'c3_1':(1072,1173),
        'c3_2':(1174,1276),'v4':(1277,1388),'c4':(1389,1481),'v5':(1482,1547),'c5':(1548,1700),'NHR':(1701,1860)}

        self.base_dir_seq="/home/ericfournier/ProjetPersonnel/ProjetsNGS/HIV/SlidingWindow/GlobalAlignment/"
        self.base_dir_ref="/home/ericfournier/ProjetPersonnel/ProjetsNGS/HIV/SlidingWindow/TempReference/"
        self.base_dir_map="/media/ericfournier/Elements/MiSeq/HIV_slidingwindow_map/"
        self.base_dir_stat="/home/ericfournier/ProjetPersonnel/ProjetsNGS/HIV/SlidingWindow/StatRes/"

        self.yam_file=open('/home/ericfournier/ProjetPersonnel/ProjetsNGS/HIV/SlidingWindow/HIV.yaml').read()  # read yaml file

        self.NbSpec=yaml.load(self.yam_file)[0]  #number of specimen
        self.SpecList=yaml.load(self.yam_file)[1:self.NbSpec+1]  #list of specimen id

        #self.global_align_iter=SeqIO.parse(self.base_dir_seq+'GlobalHIVAlign_test.fas','fasta')
        self.global_align_iter=SeqIO.parse(self.base_dir_seq+'GlobalHIVAlign_original.fas','fasta')
        ref_rec=self.global_align_iter.next()

        self.ComputeRegionLength(ref_rec)


        self.mapper=Mapper(self.base_dir_map)
        self.stat_comput=StatComputer(self.base_dir_map,self.base_dir_stat,self.base_dir_ref,self.gp120_region,self.gp120_region_length_inref)

        self.StartAnalysis()

        self.SaveStat()

    def ComputeRegionLength(self,ref_rec):
        for region in self.gp120_region:
            length_with_gap=len(ref_rec.seq[self.gp120_region[region][0]-1:self.gp120_region[region][1]])
            nb_gap=ref_rec.seq[self.gp120_region[region][0]-1:self.gp120_region[region][1]].count('-')

            #print nb_gap

            self.gp120_region_length_inref[region]=length_with_gap-nb_gap

        #print self.gp120_region_length_inref

    def StartAnalysis(self):

        for my_rec in self.global_align_iter:
            self.stat_comput.spec_list.append(my_rec.id)
            #print len(my_rec.seq)
            #self.CreateOneSpecFasta(my_rec)
            #self.mapper.IndexRef(self.base_dir_ref+my_rec.id+'.fasta')
            #self.mapper.FindFastqPair(my_rec.id)
            #self.mapper.Map(my_rec.id,self.base_dir_ref+my_rec.id)

            self.stat_comput.CreateGapCountDict(my_rec)
            #print my_rec.id
            self.stat_comput.CreateRegionsInMapFile(self.gp120_region_ordered_list,my_rec.id)

        self.stat_comput.InitAllStatDict()

        #for region in ['v3']:
        for region in self.gp120_region_ordered_list:
            print "REGION IS: ",region
            self.stat_comput.AnalyseReadsInRegions(region)

    def SaveStat(self):
        self.stat_comput.SaveStatOutput()

    def CreateOneSpecFasta(self,rec):

        rec.seq=Seq(re.sub('-','',str(rec.seq)))
        SeqIO.write(rec,'{1}{0}.fasta'.format(rec.id,self.base_dir_ref),'fasta')


my_wa=WindowAnalyser()
