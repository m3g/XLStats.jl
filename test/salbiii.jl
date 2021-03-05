
using XLStats

data_dir="./data"
links = read_all(xml_file="$data_dir/salbiii_hitsDetail.dat",
                 topolink_log="$data_dir/salbiii_topolink.log",
                 topolink_input="$data_dir/topolink.inp",
                 #xic_file_name="$data_dir/salbiii_xic.dat",
                 domain=2:134)

 




