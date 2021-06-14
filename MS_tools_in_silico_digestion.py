import multiprocessing
import re
import time
#import cPickle as pkl
#from isotopic_mass_simple import mass_calc


def create_fasta_dict(file_name):   #generates a dictonary from a fasta file
    dict_fasta = {}
    f = open(file_name,"r")     #opens the file

    fasta_file = f.read().split("\n>")      #reads the file and splits into chunks - \n> is before every entry
    for chunk in fasta_file:        # go though chunk by chunk
        lines=chunk.split('\n')         #split by lines - first line contains protein name and other junk
        try:
            uniprot_name=lines[0].split("|")[1]         # most have | around them, but not all. if they do make the uniprot_ID variable
        except:
            uniprot_name=lines[0]       # if there is no | - just use the whole line
        seq = "".join(lines[1:])      # join all but the first line together to give the sequence
        dict_fasta[uniprot_name]= seq       # create dictionary
    f.close()

    return dict_fasta


def peptide_generator(input):
    from MS_tools_parameters import missed  # parameter files contain all the changable user variables and are imported here
    from MS_tools_parameters import expasy_rules  # regex rules are imported for digestion
    from MS_tools_parameters import min_len  # min length of peptides desired
    from MS_tools_parameters import max_len  # MAX length of peptides desired
    from MS_tools_parameters import specificity
    from MS_tools_parameters import custom_rules
    # print ('Specificity: %i' % specificity)
    protein_ID = input[0]  # take the protein ID from the pool - input is a tuple in the form (protein_ID,Sequence)
    seq = input[1] # take the sequence from the pool

    enzyme = 'custom_trypsin_human4h_1'

    list_cuts = [m.end() for m in re.finditer(custom_rules[enzyme], seq)] # Create a list of cuts in the protein - does not include the start and end
    list_cuts.insert(0, 0)
    list_cuts.append(len(seq))  # add start and end to the cuts list
    if specificity == 1: # if semi is selected in parameters.py

        """this generates a list of start and end points for peptides in the sequence in the format (start,end) 
        GOING FORWARDS from tryptic ends (starts at the first cut NOT 0)"""
        list_peptide_coordinates_forwards = [(list_cuts[j], list_cuts[j] + k)              # format (start,end) iterating through k up to the next cut point
                                        for i in range(1, missed + 2)       # (or one after up to the number of missed)
                                        for j in range(1, len(list_cuts) - i)    # j corresponds to the current cut location in the cuts list and iterates up to the end - the number of missed cuts)

                                        for k in range((list_cuts[j + i - 1]) - list_cuts[j] + 1, (list_cuts[j + i] + 1) - list_cuts[j])  # iterates though k from 1 to the next cut point - the range is modified to not repeat the ones from lower values of j
                                        if  min_len <= k <= max_len] # will only calculate lengths between the min and max length

        """this generates a list of start and end points for peptides in the sequence in the format (start,end) 
               GOING BACKWARDS from tryptic ends (starts at the first cut removes iterated values of k and stops at the last cut)"""
        list_peptide_coordinates_backwards = [((list_cuts[j + i] - k), list_cuts[j + i]) for i in range(1, missed + 2) for j in
                                         range(0, (len(list_cuts) - i - 1)) for k in
                                         range(list_cuts[j + i] - list_cuts[j + 1] + 1, list_cuts[j + i] - (list_cuts[j]) + 1) if
                                         k >= min_len and k <= max_len]


        list_peptide_coordinates = list_peptide_coordinates_backwards + list_peptide_coordinates_forwards # combines the Forward and Back coordinates
        dict_peptide = {protein_ID: {seq[i[0]:i[1]] for i in list_peptide_coordinates}} # uses coordinates to generate a dictonary of peptides from the sequence


    elif specificity == 2: # if full is selected in parameters.py


        peptides = [seq[list_cuts[i]:list_cuts[i+j]] for j in range(1,missed+2) for i in range(0, len(list_cuts)-j) if min_len <= (list_cuts[i+j] - list_cuts[i]) <= max_len ]
        #Creates list of coorinates similar to above, but without the iteration with k
        peptide_set = list((set(peptides)))
        dict_peptide = {protein_ID : peptide_set}

    else: # NON tryptic


        list_peptide_coordinates_forwards = [seq[0:i] for i in range(1,len(seq)+1) if i >= min_len and i <= max_len]
        list_peptide_coordinates_backwards = [seq[i:len(seq)]for i in range(1,len(seq)) if (len(seq)-i) >= min_len and (len(seq)-i) <= max_len]

        peptides = list_peptide_coordinates_backwards + list_peptide_coordinates_forwards

        peptide_set = list((set(peptides)))

        dict_peptide = {protein_ID: peptide_set}

    return dict_peptide #returns dictionary to pool




if __name__ == '__main__':
    import MS_tools_parameters
    from protein_coverage import fasta_reader2, fasta_reader
    import pickle as ppp
    start_time = time.time()
    fasta_file = 'D:/data/proteome_fasta/uniprot-proteome_UP000005640_reverse_beta.fasta'
    # print("loading: " + MS_tools_parameters.fasta_filename)
    # dict_fasta = create_fasta_dict(MS_tools_parameters.fasta_filename) # Creates the dictonary of proteins [protein_ID : sequence]
    dict_fasta = fasta_reader(fasta_file)
    print ("done")
    print (time.time() - start_time)
    print ("calculating peptides")
    work = [(i, dict_fasta[i]) for i in dict_fasta] # from the dictonary, create a list of tuples of all proteins to be operated on format (protein_ID,sequence)

    n = len(work) / 500  # create a chunksize number for a desired number of chunks - can be changed 500-1000 seems fast
    print (n)
    with multiprocessing.Pool(multiprocessing.cpu_count()-2) as pool:
        result = pool.map(peptide_generator,work, chunksize=300)
        pool.close()
        pool.join()
     # output is a list of dictionaries
    dict_peptides = {k: v for d in result for k,v in d.items()} #convert to one dictionary
    peptide_set = {pep for v in dict_peptides.values() for pep in v}
    print (len(peptide_set))
    ppp.dump(dict_peptides, open('D:/data/deep_proteome/non_specfic_search/custom_trypsin_human4h_1.p','wb'),protocol=-1)


    """
    list_with_mass = [(mass_calc(j),j,i) for i in dict_peptides for j in dict_peptides[i]]
    print ("done")
    print (time.clock() - start_time)
    print ("sorting peptides")
    list_with_mass.sort()
    print ("done")
    print (time.clock() - start_time)
    print ("creating index")
    dict_index = {}
    j = int(list_with_mass[0][0])
    for i in range(len(list_with_mass)):
        if float(list_with_mass[i][0]) >= float(j):
            #print str(j) + " @ " + str(i)
            dict_index[int(list_with_mass[i][0])] = i
            j = int(list_with_mass[i][0]) + 1

    print ("done")
    print (time.clock() - start_time)
    print ("writing file")
    write_file = open("indexed_file.txt", "w")
    write_file.write( "  ===/INDEX_START/=== \n")
    for i in dict_index:
        write_file.write("masses " + str(i) + " to " + str(i+1) + " starting at line " + str(dict_index[i] + len(dict_index) + 2) + "\n")
    write_file.write("  ===/INDEX_END/=== \n" )

    for i in range(len(list_with_mass)):
        write_file.write(str(list_with_mass[i][0]) +"\t" + list_with_mass[i][1] + "\t" + list_with_mass[i][2] +"\n")



    write_file.close()
    print ("done")
    print (time.clock() - start_time)
    """