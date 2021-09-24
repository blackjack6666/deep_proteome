import struct
from xml.etree import cElementTree as ET
import zlib
import base64
import numpy as np
import time
import os

np_dtype_numtype = {'i': np.int32, 'e': np.single, 'f': np.float32, 'q': np.int64, 'd': np.float64}

accession_dict = {"MS:1000519": "32i",
                  "MS:1000520": "16e",
                  "MS:1000521": "32f",
                  "MS:1000522": "64i",
                  "MS:1000523": "64d",
                  "MS:1000574": "zlib",
                  "MS:1000576": "no compression",
                  "MS:1000515": "intensity",
                  "MS:1000514": "mass"}

data_type_dict = {'data_number':'0',
                  'data_encoding': '32f',
                  'data_compression': 'no compression',
                  'data_type': 'mass'}


def base64_decoder(base64_data, number_fmt, compress_method, array_type):
    if base64_data is not None:
        num_type = np_dtype_numtype[number_fmt[-1]]
        decode_base64 = base64.decodebytes(base64_data.encode('ascii'))
        if compress_method == 'zlib':
            decode_base64 = zlib.decompress(decode_base64)
        data = np.frombuffer(decode_base64, dtype=num_type)
        #if array_type == 'intensity':
        #    data = np.log2(np.where(data>0.00001, data, 0.00001)/np.linalg.norm(data))  # performs log only on intensity
    else:
        data = np.array([])
    return data


def base64_lossy_decoder(base64_data, number_fmt, compress_method, array_type):
    if base64_data is not None:
        num_type = np_dtype_numtype[number_fmt[-1]]
        decode_base64 = base64.decodebytes(base64_data.encode('ascii'))
        if compress_method == 'zlib':
            decode_base64 = zlib.decompress(decode_base64)
        data = np.frombuffer(decode_base64, dtype=num_type)
        if array_type == 'intensity':
            data = data/np.linalg.norm(data)*65535
            data = (np.round(np.log2(np.where(data>0.00001, data, 0.00001)),4)*1000).astype(np.ushort)  # performs log only on intensity
        else:
            basemass= int(data[0]*10000)
            data = (np.diff(data)*10000).astype(np.int32)
            data=np.insert(data, 0, basemass)
    else:
        data = np.array([])
    return data


def base64_encoder(number_array: np.ndarray, compress_method: str):
    byte_data = number_array.tobytes()
    if compress_method == 'zlib':
        byte_data = zlib.compress(byte_data)
    return base64.b64encode(byte_data).decode('ascii').replace('\n', '')


def mzml_splitter(mzml_file: str):
    start=time.time()
    # Generate file names
    int_binary_file = mzml_file.replace('.mzML', '.bint')
    mass_binary_file = mzml_file.replace('.mzML', '.bmass')
    mzml_out = mzml_file.replace('.mzML', '.smzml')

    # Open file pointers
    mass_binary_out_fp = open(mass_binary_file, 'wb')
    int_binary_out_fp = open(int_binary_file, 'wb')
    mzml_out_fp = open(mzml_out, 'w', newline='\n', encoding='utf-8')

    # Create empty list for data temporary storage
    data_format={}
    data_position=[]
    # Iterate mzML file and get all data into data_position[]
    for event, elem in iter(ET.iterparse(mzml_file, events=('end',))):
        if elem.tag.endswith('}cvParam'):
            if elem.get('accession').endswith(('MS:1000521', 'MS:1000522', 'MS:1000523', 'MS:1000519', 'MS:1000520')):  # retrieves the datatype based on MS accession
                data_type_dict['data_encoding'] = accession_dict[elem.get('accession')]
            if elem.get('accession').endswith(('MS:1000576', 'MS:1000574')):  # retrieves the compression based on MS accession
                data_type_dict['data_compression'] = accession_dict[elem.get('accession')]
            if elem.get('accession').endswith(('MS:1000515', 'MS:1000514')):  # retrives array_type
                data_type_dict['data_type'] = accession_dict[elem.get('accession')]

        elif elem.tag.endswith('}binary'):
            number_array = base64_lossy_decoder(elem.text, data_type_dict['data_encoding'], data_type_dict['data_compression'], data_type_dict['data_type'])
            data_type_dict['data_number']=str(len(number_array))
            data_position.append(([data_type_dict['data_number'], data_type_dict['data_encoding'], data_type_dict['data_compression'], data_type_dict['data_type']], number_array))


    # Split mzml into smzml, binary mass and binary int
    with open(mzml_file, 'r', encoding='utf-8') as mzml_in_fp:
        i=0
        for line in mzml_in_fp:
            if '<binary>' in line and '</binary>' in line:
                mzml_out_fp.write(line.split('<binary>')[0])
                mzml_out_fp.write("<binary>$%s$</binary>\n" % '$'.join(data_position[i][0]))
                number_array=data_position[i][1]
                number_fmt = np_dtype_numtype[data_position[i][0][1][-1]]
                if data_position[i][0][3]=='mass':
                    mass_binary_out_fp.write(number_array.astype(number_fmt))
                else:
                    int_binary_out_fp.write(number_array.astype(number_fmt))
                i+=1
            else:
                mzml_out_fp.write(line)

    # Close file pointers
    mzml_out_fp.close()
    int_binary_out_fp.close()
    mass_binary_out_fp.close()

    print(time.time()-start)


def mzml_lossy_splitter(mzml_file: str):
    # Generate file names
    int_binary_file = mzml_file.replace('.mzML', '.bint')
    mass_binary_file = mzml_file.replace('.mzML', '.bmass')
    mzml_out = mzml_file.replace('.mzML', '.smzml')

    # Open file pointers
    mass_binary_out_fp = open(mass_binary_file, 'wb')
    int_binary_out_fp = open(int_binary_file, 'wb')
    mzml_out_fp = open(mzml_out, 'w', newline='\n')

    # Create empty list for data temporary storage
    data_position = []

    # Iterate mzML file and get all data into data_position[]
    for event, elem in iter(ET.iterparse(mzml_file, events=('end',))):
        if elem.tag.endswith('}cvParam'):
            if elem.get('accession').endswith(('MS:1000521', 'MS:1000522', 'MS:1000523', 'MS:1000519', 'MS:1000520')):  # retrieves the datatype based on MS accession
                data_type_dict['data_encoding'] = accession_dict[elem.get('accession')]
            if elem.get('accession').endswith(('MS:1000576', 'MS:1000574')):  # retrieves the compression based on MS accession
                data_type_dict['data_compression'] = accession_dict[elem.get('accession')]
            if elem.get('accession').endswith(('MS:1000515', 'MS:1000514')):  # retrives array_type
                data_type_dict['data_type'] = accession_dict[elem.get('accession')]

        elif elem.tag.endswith('}binary'):
            number_array = base64_lossy_decoder(elem.text, data_type_dict['data_encoding'], data_type_dict['data_compression'], data_type_dict['data_type'])
            data_type_dict['data_number']=str(len(number_array))
            data_position.append(([data_type_dict['data_number'], data_type_dict['data_encoding'], data_type_dict['data_compression'], data_type_dict['data_type']], number_array))

    # Split mzml into smzml, binary mass and binary int
    with open(mzml_file, 'r') as mzml_in_fp:
        i=0
        for line in mzml_in_fp:
            if '<binary>' in line and '</binary>' in line:
                mzml_out_fp.write(line.split('<binary>')[0])
                mzml_out_fp.write("<binary>$%s_lossy$</binary>\n" % '$'.join(data_position[i][0]))
                number_array=data_position[i][1]
                if data_position[i][0][3]=='mass':
                    mass_binary_out_fp.write(number_array.astype(np.int32))
                else:
                    int_binary_out_fp.write(number_array.astype(np.ushort))
                i+=1
            else:
                mzml_out_fp.write(line)

    # Close file pointers
    mzml_out_fp.close()
    int_binary_out_fp.close()
    mass_binary_out_fp.close()


def mzml_decoder(smzml_file, bmass_file, bint_file, mzml_file):
    # Create file pointers
    mass_in_fp = open(bmass_file, 'rb')
    int_in_fp = open(bint_file, 'rb')
    smzml_file = open(smzml_file, 'r')

    # Restore mzML file from smzml
    with open(mzml_file, 'w', newline='\n') as f_out:
        for line in smzml_file:
            if '<binary>' in line and '</binary>' in line:
                f_out.write(line.split('<binary>')[0])
                _, data_num, data_fmt, data_compression, data_type, _ = line.split('$')
                if data_type == 'mass':
                    number_array =np.frombuffer(mass_in_fp.read(int(int(data_fmt[:2])*int(data_num)/8)), np_dtype_numtype[data_fmt[-1]])
                else:
                    number_array = np.frombuffer(int_in_fp.read(int(int(data_fmt[:2])*int(data_num)/8)), np_dtype_numtype[data_fmt[-1]])
                f_out.write('<binary>%s</binary>\n'%base64_encoder(number_array, data_compression))
            else:
                f_out.write(line)

    # Close file pointers
    mass_in_fp.close()
    int_in_fp.close()
    smzml_file.close()


if __name__ == '__main__':
    mzml_file = r'K:\test\64bit_zlib_profile\fusion_20200116_qjl_YLD_19.mzML'
    smzml_file = r'J:\mass_cloud\projects_folder\1000111\test\fusion_20200116_qjl_YLD_19.smzml'
    bmass_file = r'J:\mass_cloud\projects_folder\1000111\test\fusion_20200116_qjl_YLD_19.bmass'
    bint_file = r'J:\mass_cloud\projects_folder\1000111\test\fusion_20200116_qjl_YLD_19.bint'
    mzml_out = r'J:\mass_cloud\projects_folder\1000111\test\test2.mzML'

    # Split to lossless files
    mzml_splitter(mzml_file)

    # Split to lossy files
    # mzr.mzml_lossy_splitter(mzml_file)

    # Re-generate mzML file
    # mzr.mzml_decoder(smzml_file, bmass_file, bint_file, mzml_out)
