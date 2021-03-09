
using XLStats


data_dir="./data"
xml_data="$data_dir/hitsDetail.xls"
xml_data="$data_dir/salbiii_hitsDetail.dat"
links = read_all(xml_file=xml_data,
                 topolink_log="$data_dir/salbiii_topolink.log",
                 topolink_input="$data_dir/topolink.inp",
                 xic_file_name="$data_dir/salbiii_xic.dat",
                 domain=2:134)

 




