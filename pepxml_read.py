from collections import Counter
from tsv_reader import venn_diagram_gen
import re
from tensorflow import keras
from tcn import TCN

# pxm_1 = 'D:/data/mscompress/test/163_3A_1/163-3A_1.pepXML'
# pxm_2 = 'D:/data/mscompress/test/062902021/decomp.pepXML'
#
# from tsv_reader import pepxml_peptide_getter
#
# pxm_original = pepxml_peptide_getter(pxm_1)
# pxm_decom = pepxml_peptide_getter(pxm_2)
# pxm_original_counter_dict = Counter(pxm_original)
# pxm_decom_counter_dict = Counter(pxm_decom)
#
# pxm_original_counter_list = [each+str(i) for each in pxm_original_counter_dict for i in range(pxm_original_counter_dict[each])]
# pxm_decom_counter_list = [each+str(i) for each in pxm_decom_counter_dict for i in range(pxm_decom_counter_dict[each])]
#
# venn_diagram_gen({"original":pxm_original_counter_list,"decomp":pxm_decom_counter_list})
#
# unique = [each for each in pxm_decom if each not in pxm_original]
#
# print (unique)


reg_replace = {r'\<cvParam.*?MS\:1000579.*?\"\/\>',
               r'\<cvParam.*?MS\:1000130.*?\"\/\>',
               r'\<cvParam.*?MS\:1000127.*?\"\/\>',
               r'\<cvParam.*?MS\:1000285.*?\"\/\>',
               r'\<cvParam.*?MS\:1000528.*?\"\/\>',
               r'\<cvParam.*?MS\:1000527.*?\"\/\>',
               r'\<cvParam.*?MS\:1000796.*?\"\/\>',
               r'\<cvParam.*?MS\:1000795.*?\"\/\>',
               r'\<cvParam.*?MS\:1000512.*?\"\/\>',
               r'\<cvParam.*?MS\:1000616.*?\"\/\>',
               r'\<cvParam.*?MS\:1000927.*?\"\/\>',
               r'\<cvParam.*?MS\:1000500.*?\"\/\>',
               r'\<userParam name.*?type.*?\"\/\>',
               r'\<cvParam.*?MS\:1000042.*?\"\/\>',
               r'\<cvParam.*?MS\:1000422.*?\"\/\>',
               r'\<cvParam.*?MS\:100045.*?\"\/\>',
               r'(?<!\>)\r\n'}


ms_accession_necessary = {'MS:1000511':'ms level',
                          'MS:1000016':'scan start time',
                          'MS:1000501':'scan window lower limit',
                          'MS:1000500':'scan window upper limit',
                          'MS:1000523':'64-bit float',
                          'MS:1000574':'zlib compression',
                          'MS:1000514':'m/z array',
                          'MS:1000515':'intensity array',
                          'MS:1000827':'isolation window target m/z',
                          'MS:1000828':'isolation window lower offset',
                          'MS:1000829':'isolation window upper offset',
                          'MS:1000744':'selected ion m/z',
                          'MS:1000041':'charge state',
                          'MS:1000505':'base peak intensity'}

unit_accession_common = {'UO:0000028':'millisecond',
                         'UO:0000031':'minute',
                         'UO:0000266':'electronvolt',
                         'MS:1000040':'m/z',
                         'MS:1000131':'number of detector counts'}




# mzml = 'D:/data/mscompress/163-3A_1.mzML'
# output = 'D:/data/mscompress/163-3A_1_tag_delete.mzML'
#
# with open(mzml,'r') as f_read:
#     f_read = f_read.read()
#     for each in reg_replace:
#         print (each)
#         f_read = re.sub(each,'',f_read)
#
# with open(output,'w') as f_write:
#     f_write.write(f_read)



# unit_acces_dict = {}
# with open('D:/data/mscompress/uo.obo.txt','r') as f_open:
#     file_read = f_open.read()
#     for each in file_read.split('[Term]\n')[1:]:
#         if 'unit_slim' in each:
#             id_name = each.split('id: ')[1].split('\n')[0]
#             name = each.split('name: ')[1].split('\n')[0]
#
#             unit_acces_dict[id_name] = name

# with open('D:/data/mscompress/uo.txt','w',newline='\n') as f_write:
#     for uo in unit_acces_dict:
#         f_write.write(uo+'\t'+unit_acces_dict[uo]+'\n')



model_path = 'D:/data/deep_proteome/bilstm_tcn_suba_06302021_15epoch'
model = keras.models.load_model('D:/data/deep_proteome/bilstm_tcn_tryp_ct_suba_06302021_15epoch',custom_objects={'TCN': TCN})
print (model.summary())
