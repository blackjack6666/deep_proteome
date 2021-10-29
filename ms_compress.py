import mmap
import os
import numpy as np
import time
# import pybase64
import zlib
from xml.etree import cElementTree as ET
import pip

try:
    import pybase64
except ImportError:
    os.system("python3 -m pip install pybase64")
    import pybase64


np_dtype_numtype = {'i': np.int32, 'e': np.single, 'f': np.float32, 'q': np.int64, 'd': np.float64}


def chunks(l, n):
    """Yield successive n-sized chunks from l."""
    for i in range(0, len(l), n):
        yield l[i:i + n]


def decode_pos(data_position, mzml_fp, data_fmt, array_type):
    data_start, data_end = data_position
    mzml_fp.seek(data_start)
    data = mzml_fp.read(data_end - data_start)
    num_data, num_array = base64_decoder(data, data_fmt['data_encoding'], data_fmt['data_compression'], array_type)
    return num_data, num_array


def base64_decoder(base64_data: bytes, number_fmt, compress_method, array_type):
    if base64_data is not None:
        num_type = np_dtype_numtype[number_fmt[-1]]
        decode_base64 = pybase64._pybase64.b64decode_as_bytearray(base64_data)
        if compress_method == 'zlib':
            decode_base64 = zlib.decompress(decode_base64)
        data = np.frombuffer(decode_base64, dtype=num_type)
        # if array_type == 'intensity':
        #    data = np.log2(np.where(data>0.00001, data, 0.00001)/np.linalg.norm(data))  # performs log only on intensity
    else:
        data = 0, np.array([])
    return len(data), data


def base64_encoder(number_array: np.ndarray, compress_method: str):
    byte_data = number_array.tobytes()
    if compress_method == 'zlib':
        byte_data = zlib.compress(byte_data)
    return pybase64.b64encode(byte_data).decode('ascii').replace('\n', '')


accession_dict = {"MS:1000519": "32i",
                  "MS:1000520": "16e",
                  "MS:1000521": "32f",
                  "MS:1000522": "64i",
                  "MS:1000523": "64d",
                  "MS:1000574": "zlib",
                  "MS:1000576": "no compression",
                  "MS:1000515": "intensity",
                  "MS:1000514": "mass"}


def find_string(mzml_read_fp, match_tag_start, match_tag_end, data_format, spec_no):
    start_time = time.time()
    file_name = mzml_read_fp.name
    file_size = os.path.getsize(file_name)
    encoded_start_tag = match_tag_start.encode('utf-8')
    encoded_end_tag = match_tag_end.encode('utf-8')
    len_start_tag = len(match_tag_start)
    len_end_tag = len(match_tag_end)
    data_positions = []
    mass_fmt = ("$%s$" % '$'.join([data_format['mass']['data_encoding'], data_format['mass']['data_compression'], data_format['mass']['data_type']])).encode('utf-8')
    int_fmt = ("$%s$" % '$'.join([data_format['intensity']['data_encoding'], data_format['intensity']['data_compression'], data_format['intensity']['data_type']])).encode('utf-8')
    i = 0
    with open(file_name.replace('.mzML', '.smzml'), 'wb') as fo:
        # memory-map the file, size 0 means whole file
        m = mmap.mmap(mzml_read_fp.fileno(), 0, access=mmap.ACCESS_READ)
        last_end = 0
        while True:
            if i == spec_no * 2:
                print("total time used to position all b64 data:", time.time() - start_time)
                fo.write(m.read(file_size - m.tell()))
                return data_positions
            start = m.find(encoded_start_tag)
            fo.write(m.read(start + len_start_tag - last_end))
            m.seek(start)
            end = m.find(encoded_end_tag)
            data_positions.append((start + len_start_tag, end))
            m.seek(end + len_end_tag)
            if i % 2 == 0:
                fo.write(mass_fmt)
            else:
                fo.write(int_fmt)
            fo.write(b'</binary>')
            last_end = end + len_end_tag
            i += 1


def mzml_splitter(mzml_file: str):
    start = time.time()
    # Generate file names
    int_binary_file = mzml_file.replace('.mzML', '.bint')
    mass_binary_file = mzml_file.replace('.mzML', '.bmass')
    # mzml_out = mzml_file.replace('.mzML', '.smzml')

    # Open file pointers
    mass_binary_out_fp = open(mass_binary_file, 'wb')
    int_binary_out_fp = open(int_binary_file, 'wb')
    mzml_read_fp = open(mzml_file, 'rb')
    # mzml_out_fp = open(mzml_out, 'w', newline='\n', encoding='utf-8')

    data_format = {}

    # Iterate mzML file and get all data into data_position[]
    current_encoding = '32f'
    current_compress = 'no compression'
    spec_no = 0
    for event, elem in iter(ET.iterparse(mzml_file, events=('start',))):
        if elem.tag.endswith('}cvParam'):
            if elem.get('accession').endswith(('MS:1000521', 'MS:1000522', 'MS:1000523', 'MS:1000519', 'MS:1000520')):  # retrieves the datatype based on MS accession
                current_encoding = accession_dict[elem.get('accession')]
            if elem.get('accession').endswith(('MS:1000576', 'MS:1000574')):  # retrieves the compression based on MS accession
                current_compress = accession_dict[elem.get('accession')]
            if elem.get('accession').endswith(('MS:1000515', 'MS:1000514')):  # retrives array_type
                current_type = accession_dict[elem.get('accession')]
                data_format[current_type] = {'data_encoding': current_encoding, 'data_compression': current_compress, 'data_type': current_type}

        elif elem.tag.endswith('}spectrumList'):
            spec_no = int(elem.get('count'))
            print("Total binary data:", spec_no * 2)

        elif len(data_format.keys()) == 2:
            print("Mass and intensity data format", data_format)
            break

    data_position = find_string(mzml_read_fp, '<binary>', '</binary>', data_format, spec_no)
    data_chunks = chunks(data_position, 2)

    mass_num_data_list = []
    int_num_data_list = []

    for each_data in data_chunks:
        mass_pos, int_pos = each_data

        mass_num_data, mass_data_to_write = decode_pos(mass_pos, mzml_read_fp, data_format['mass'], 'mass')
        mass_binary_out_fp.write(mass_data_to_write)
        mass_num_data_list.append(mass_num_data)

        int_num_data, int_data_to_write = decode_pos(int_pos, mzml_read_fp, data_format['intensity'], 'intensity')
        int_binary_out_fp.write(int_data_to_write)
        int_num_data_list.append(int_num_data)

    mass_num_data_list.append(len(mass_num_data_list))
    int_num_data_list.append(len(int_num_data_list))

    mass_binary_out_fp.write(np.array(mass_num_data_list, dtype=np.int32).tobytes())
    int_binary_out_fp.write(np.array(int_num_data_list, dtype=np.int32).tobytes())

    mzml_read_fp.close()
    int_binary_out_fp.close()
    mass_binary_out_fp.close()


def mzml_lossy_splitter(mzml_file: str):
    start = time.time()
    # Generate file names
    int_binary_file = mzml_file.replace('.mzML', '.bint')
    mass_binary_file = mzml_file.replace('.mzML', '.bmass')
    # mzml_out = mzml_file.replace('.mzML', '.smzml')

    # Open file pointers
    mass_binary_out_fp = open(mass_binary_file, 'wb')
    int_binary_out_fp = open(int_binary_file, 'wb')
    mzml_read_fp = open(mzml_file, 'rb')
    # mzml_out_fp = open(mzml_out, 'w', newline='\n', encoding='utf-8')

    data_format = {}

    # Iterate mzML file and get all data into data_position[]
    current_encoding = '32f'
    current_compress = 'no compression'
    spec_no = 0
    for event, elem in iter(ET.iterparse(mzml_file, events=('start',))):
        if elem.tag.endswith('}cvParam'):
            if elem.get('accession').endswith(('MS:1000521', 'MS:1000522', 'MS:1000523', 'MS:1000519', 'MS:1000520')):  # retrieves the datatype based on MS accession
                current_encoding = accession_dict[elem.get('accession')]
            if elem.get('accession').endswith(('MS:1000576', 'MS:1000574')):  # retrieves the compression based on MS accession
                current_compress = accession_dict[elem.get('accession')]
            if elem.get('accession').endswith(('MS:1000515', 'MS:1000514')):  # retrives array_type
                current_type = accession_dict[elem.get('accession')]
                data_format[current_type] = {'data_encoding': current_encoding, 'data_compression': current_compress, 'data_type': current_type}

        elif elem.tag.endswith('}spectrumList'):
            spec_no = int(elem.get('count'))
            print("Total binary data:", spec_no * 2)

        elif len(data_format.keys()) == 2:
            print("Mass and intensity data format", data_format)
            break

    data_position = find_string(mzml_read_fp, '<binary>', '</binary>', data_format, spec_no)
    data_chunks = chunks(data_position, 2)

    mass_num_data_list = []
    int_num_data_list = []

    for each_data in data_chunks:
        mass_pos, int_pos = each_data

        mass_num_data, mass_data_to_write = decode_pos(mass_pos, mzml_read_fp, data_format['mass'], 'mass')
        mass_binary_out_fp.write(mass_data_to_write)
        mass_num_data_list.append(mass_num_data)

        int_num_data, int_data_to_write = decode_pos(int_pos, mzml_read_fp, data_format['intensity'], 'intensity')
        int_binary_out_fp.write(int_data_to_write)
        int_num_data_list.append(int_num_data)

    mass_num_data_list.append(len(mass_num_data_list))
    int_num_data_list.append(len(int_num_data_list))

    mass_binary_out_fp.write(np.array(mass_num_data_list, dtype=np.int32).tobytes())
    int_binary_out_fp.write(np.array(int_num_data_list, dtype=np.int32).tobytes())

    mzml_read_fp.close()
    int_binary_out_fp.close()
    mass_binary_out_fp.close()


def mzml_decoder(smzml_file, bmass_file, bint_file, mzml_file):
    # Create file pointers
    mass_in_fp = open(bmass_file, 'rb')
    int_in_fp = open(bint_file, 'rb')
    smzml_file = open(smzml_file, 'r', encoding='utf-8')

    int_in_fp.seek(-4, os.SEEK_END)
    total_spec_no = int.from_bytes(int_in_fp.read(4), byteorder='little')
    int_in_fp.seek(-4 * (total_spec_no + 1), os.SEEK_END)
    spec_no_array = np.frombuffer(int_in_fp.read(4 * total_spec_no), np.int32)
    int_in_fp.seek(0)

    i = 0
    # Restore mzML file from smzml
    with open(mzml_file, 'w', newline='\n', encoding='utf-8') as f_out:
        for line in smzml_file:
            if '<binary>' in line and i < total_spec_no:
                f_out.write(line.split('<binary>')[0])
                _, data_fmt, data_compression, data_type, _ = line.split('$')
                data_num = spec_no_array[i]
                if data_type == 'mass':
                    number_array = np.frombuffer(mass_in_fp.read(int(int(data_fmt[:2]) * int(data_num) / 8)), np_dtype_numtype[data_fmt[-1]])
                else:
                    read_number=int(int(data_fmt[:2]) * int(data_num) / 8)
                    number_array = np.frombuffer(int_in_fp.read(read_number), np_dtype_numtype[data_fmt[-1]])
                    i += 1
                f_out.write('<binary>%s</binary>\n' % base64_encoder(number_array, data_compression))

            else:
                f_out.write(line)

    # Close file pointers
    mass_in_fp.close()
    int_in_fp.close()
    smzml_file.close()




if __name__ == '__main__':
    # mzml_file = r'K:\test\64bit_zlib_profile\fusion_20200116_qjl_YLD_19.mzML'
    # smzml_file = r'J:\mass_cloud\projects_folder\1000111\test\fusion_20200116_qjl_YLD_19.smzml'
    # bmass_file = r'J:\mass_cloud\projects_folder\1000111\test\fusion_20200116_qjl_YLD_19.bmass'
    # bint_file = r'J:\mass_cloud\projects_folder\1000111\test\fusion_20200116_qjl_YLD_19.bint'
    # mzml_out = r'J:\mass_cloud\projects_folder\1000111\test\test2.mzML'

    __status__ = "Development"
    __version__ = "0.0.1"

    # Split to lossless files
    # mzml_splitter(mzml_file)

    # Split to lossy files
    # mzr.mzml_lossy_splitter(mzml_file)

    # Re-generate mzML file
    # mzr.mzml_decoder(smzml_file, bmass_file, bint_file, mzml_out)

    import argparse

    # Creating the parser object
    parser = argparse.ArgumentParser(
        description="MSCompress is used for highly efficient compression of mass spec raw data"
    )


    # Add required argument group
    required = parser.add_argument_group("required arguments")
    required.add_argument(
        "-i","--mzml", help="mzml file path",
        required=True
    )
    required.add_argument(
        "-o", "--mzs", help="output mzs file path",
        required=True
    )

    # add optional argument group
    optional = parser.add_argument_group("user-optional arguments")
    # Adding the version option
    optional.add_argument(
        "-v", "--version", action="version",
        version="%(prog)s {}".format(__version__)
    )
    optional.add_argument(
        "--spec_no", help="specify spec number or spec number range. if not specified, process all"
    )
    optional.add_argument(
        "--spec_type", help="MS1, MS2 or all", default="all", choices=('MS1','MS2','all')
    )
    optional.add_argument(
        "--loss_type", help="lossless or lossy compression", default="lossless", choices=("lossless",'lossy')
    )

    # parse arguments
    args = parser.parse_args()

    mzml_input = args.mzml
    mzs_output = args.mzs

    if args.spec_no:
        pass

    if args.spec_type:
        pass

    int_binary_file = mzml_input.replace('.mzML', '.bint')
    mass_binary_file = mzml_input.replace('.mzML', '.bmass')
    smzml = mzml_input.replace('.mzML', '.smzml')

    if args.loss_type == 'lossless':
        mzml_splitter(mzml_input)
    elif args.loss_type == 'lossy':
        mzml_lossy_splitter(mzml_input)
    else:
        raise ValueError("only lossless and lossy compression could be performed")

    mzml_decoder(smzml,mass_binary_file,int_binary_file,mzs_output)


    # import pandas as pd
    # compressed_time = 'D:/data/mscompress/testfile_decompresstime.xlsx'
    # df = pd.read_excel(compressed_time,index_col=0)
    # df_new = pd.DataFrame(index=df.index)
    # for each in df:
    #     if each == "convert type":
    #         df_new[each] = df[each]
    #     else:
    #         time_list = df[each].tolist()
    #         new_time_list = []
    #         for t in time_list:
    #             t = t.replace(',','').replace(' ','')
    #             if 'm' in t and 's' in t:
    #                 t = int(t.split('m')[0])*60+float(t.split('m')[1].split('s')[0])
    #             elif ':' in t:
    #                 t = int(t.split(':')[0])*60+float(t.split(':')[1])
    #             elif 's' in t and 'm' not in t:
    #                 t = float(t.split('s')[0])
    #             else:
    #                 t = float(t)
    #             new_time_list.append(t)
    #         df_new[each+'(s)'] = new_time_list
    # df_new.to_excel('D:/data/mscompress/testfile_decompresstime_cleaned.xlsx')
