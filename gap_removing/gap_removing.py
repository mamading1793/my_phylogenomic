#!/home/good/anaconda2/bin/python
# -*- coding: utf-8 -*-  
# This program is used to remove the gaps site from a file of aligned sequences
# This program should have 3 paras as follow
# -g gap percentage , how ration of the gap precent should be removed
# -c if not select para c, we default it is run by single site, if codon, we used -c to select
# -l the minimal length of the alignment not splited by gap
# -n select the species number, if not select, we assume it equte to the taxon nub


from Bio import Seq
from Bio import AlignIO
import argparse


def _parse_args():
    parser = argparse.ArgumentParser(
        usage='Filter the gaps from a multi-sequence alignment',
        description='''This program is created by mamading in 2019/03/04. 
           Any bug please contact with me, my email is zxmlmq@163.com
        ''',
        add_help=True)
    parser.add_argument("-g", type=float, dest='gap_p', required=False,
                        default=0.5, help="The percentage of the gap, default: 0.5")

    parser.add_argument("-l", type=int, dest='min_l', required=False,
                        default=1, help="The minimal length of the block,default: 1")

    parser.add_argument("-n", type=int, dest='spe_n', required=False, default=None,
                        help="The number of sepcies that should precent in the alignment")

    parser.add_argument("-c", dest='codon', required=False,
                        action='store_true', help="If we split the alignment as codons")

    parser.add_argument("-o", dest='outfile',
                        required=True, help="The output file")

    parser.add_argument("-i", dest='infile', required=True,
                        help="The input alignment file")
    return vars(parser.parse_args())


def main(infile, outfile, gap_p=0.5, spe_n=None, codon=False, min_l=1):
    '''
    infile: the input alignment file 
    outfile: the filtered output file 
    gap_p: the percentage of the gap , upper than this will be filtered 
    spe_n: the number of species selected, if not select, we read the number from the file
    codon: if we filtered as codon 
    min_l: the min length of the block'''

    seq_ali = AlignIO.read(infile, 'fasta')
    # set the para
    length = len(seq_ali[0, :])

    if spe_n == None:
        spe_n = len(seq_ali)
    gap_pn = gap_p * spe_n

    # give me the length of this alignment

    # functions
    ##########################################################
    def find_block(lista):
        '''give me the site lista, 
        return as a list of tuple as (site, block_length) '''
        saved_lista = []
        before_count = 1  # the before count also include itself
        block_lenth = 0
        k = 0  # we need to assign two pin
        j = 0
        while j < len(lista):
            k = j
            block_lenth = before_count
            tmp = 0
            while k + 1 < len(lista):
                if lista[k] + 1 == lista[k + 1]:  # it is continuous
                    block_lenth += 1
                    tmp = 1
                else:
                    break
                k += 1
            if tmp == 1:
                before_count += 1
            else:
                before_count = 1
            saved_lista.append((lista[j], block_lenth))
            j += 1
        return saved_lista
    ##########################################################

    # now we see the seq one by one
    # find the site we need to preserved
    site_lista = []
    if codon == False:
        for i in range(length):
            site = seq_ali[:, i]
            gap_number = site.count('-')  # count how many gaps
            if gap_number < gap_pn:
                site_lista.append(i)
    else:
        for i in range(length / 3):
            site = seq_ali[:, i * 3:i * 3 + 3]
            gap_number = 0
            if len(site[0, :]) != 3:
                print 'We find the sites are not tripled, so they are not codons'
                raise BaseException()
            for j in range(3):
                tmp_number = site[:, j].count('-')
                if tmp_number > gap_number:
                    gap_number = tmp_number
            if gap_number < gap_pn:
                site_lista.append(i)
        # we then expand the site lista by triple
        tmp_lista = []
        for i in site_lista:
            tmp_lista += [i * 3, i * 3 + 1, i * 3 + 2]
        site_lista = tmp_lista

    # filter the site lista with min_l

    site_lista = find_block(site_lista)  # find the length of each block
    new_lista = []
    if codon == True:  # if it is codon, we then triple the min length
        min_l *= 3
    for i in site_lista:
        if i[-1] >= min_l:
            new_lista.append(i)
    site_lista = new_lista

    if site_lista == []:  # if after filter, there is no sequence, we should output a text but with no output file
        print 'After filter there is no alignment remain in %s .' % infile
        raise BaseException()
    # now we make sure the sequence at least has one position
    # we then make a new alignmnet

    new_ali = seq_ali[:, site_lista[0][0]:site_lista[0][0] + 1]
    for i in site_lista[1:]:
        new_ali += seq_ali[:, i[0]:i[0] + 1]

    # at last we write the new alignment file
    AlignIO.write(new_ali, outfile, 'fasta')
    return True


if __name__ == '__main__':

    dic = _parse_args()
    gap_p = dic['gap_p']
    min_l = dic['min_l']
    spe_n = dic['spe_n']
    codon = dic['codon']
    outfile = dic['outfile']
    infile = dic['infile']

    main(infile=infile, outfile=outfile, gap_p=gap_p,
         spe_n=spe_n, codon=codon, min_l=min_l)
