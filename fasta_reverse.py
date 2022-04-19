"""
Generate reverse protein sequence and output to a fasta file
"""
from protein_coverage import fasta_reader, read_description


def fasta_reverse_generator(fasta_file_in, fasta_file_out):
    print ('reverse_algorithm = null, or other protease')
    # read protein sequence into dic
    protein_dict = fasta_reader(fasta_file_in)


    # read description into dic
    ID_descrip_dict = read_description(fasta_file_in)  # dictionary format: {'ID': ('sp'or'tr', description)}

    # write id and reverse sequence into fasta_file_out
    with open(fasta_file_out, 'w', newline='\n') as file_open:
        for id in protein_dict:
            seq = protein_dict[id]
            rev_seq = protein_dict[id][::-1]

            block = range(0,len(rev_seq)+60,60)

            # write forward
            file_open.write('>'+ ID_descrip_dict[id][0]+'|'+id+'|'+ID_descrip_dict[id][1]+'\n')
            for i in range(len(block)-1):
                file_open.write(seq[block[i]:block[i+1]]+'\n')

            # write reverse
            file_open.write('>Rev_'+ ID_descrip_dict[id][0]+'|'+id+'|'+ID_descrip_dict[id][1]+'\n')
            for i in range(len(block)-1):
                file_open.write(rev_seq[block[i]:block[i+1]]+'\n')
    return fasta_file_out


if __name__=='__main__':
    from glob import glob

    fasta_file_in = 'F:/deion_biotin/uniprot-proteome_UP000000589.fasta'

    fasta_reverse_generator(fasta_file_in,fasta_file_in.replace('.fasta','_avidin_luciferase_rev.fasta'))
