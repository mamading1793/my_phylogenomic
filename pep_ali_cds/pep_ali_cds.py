#!/home/good/anaconda2/bin/python
# -*- coding: utf-8 -*-  
# 
# This program is use a aligned pep file 
# to leading the alignment of a cds file 


from Bio import SeqIO
from Bio import Seq


import argparse

def _parse_args():
    parser = argparse.ArgumentParser(
        usage = '\nThis program is used to use an aligned peptide sequence as a guild, to align its counterpart cds',
        description = 'This program is created by mamading in 2019/03/01. \n Any bug please contact with me, my email zxmlmq@163.com',
        add_help=True)

    parser.add_argument("-p", dest = 'peptide', required = True ,help="The aligned peptide file")
    parser.add_argument("-c", dest = 'cds', required = True, help="The aligned cds file")
    parser.add_argument("-o", dest = 'output', required = True, help = "The output file of the cds")
    #parser.add_argument('-h', '--help', help="Show this help message and exit")
    return vars(parser.parse_args())



def main(args):
    '''
     read seq alignment in the infile_ali file 
     this is the alignment of peptide 
     so we can see what it looks likc 
    '''
    # we build a dic for the aligned seq 

    infile_ali = args['peptide']
    infile_cds = args['cds']
    out_file = args['output']
    ali_dic = {i.name:str(i.seq) for i in SeqIO.parse(infile_ali,'fasta')}

    # read the file of cds 
    cds_dic = {i.name:str(i.seq) for i in SeqIO.parse(infile_cds,'fasta')}
    order = [i.name for i in SeqIO.parse(infile_cds,'fasta')]

    for i in ali_dic:
        if cds_dic.has_key(i):
            pass
        else:
            print 'ERROR: The sequence in aligned pep file do not has a corespondant match in cds file! It is %s' % i
            raise BaseException()
        gapped_seq = ''
        if len(cds_dic[i]) % 3 != 0:
            print 'ERROR: The cds sequence should be a multiple of 3. It is %s' % i
            raise BaseException()
        k = 0 # this para is used to cound whihc postition i have made 
        for j in ali_dic[i]:
            if j == '-':
                gapped_seq += '---'
            else:
                gapped_seq += cds_dic[i][k*3:k*3+3]
                k+=1
        cds_dic[i] = gapped_seq 

    with open(out_file, 'w') as fila:
        for i in order:
            fila.write('>%s\n%s\n' % (i, cds_dic[i]))
    return None


if __name__ =='__main__':
    args = _parse_args()
    main(args)
