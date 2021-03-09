using XLStats
using Test

using XLStats


data_dir="./data"
xml_data="$data_dir/hitsDetail.xls"
xml_data="$data_dir/salbiii_hitsDetail.dat"
links = read_all(xml_file=xml_data,
                 topolink_log="$data_dir/salbiii_topolink.log",
                 topolink_input="$data_dir/topolink.inp",
                 xic_file_name="$data_dir/salbiii_xic.dat",
                 domain=2:134)

@testset "XLStats.jl" begin

  @test name(links[1]) == "K17-K6"

  @test indexes(links[13]) == (133, 17)

  @test resnames(links[45]) == ("GLU","GLU")

  @test length(links) == 102

  @test sum(nscans(links)) == 1125

  @test count(consistency(links)) == 22

  @test sum(maxscore1(links)) ≈ 340.29699999999997

  @test sum(maxscore2(links)) ≈ 138.113

  @test sum(avgscore1(links)) ≈ 154.42675260891815

  @test sum(avgscore2(links)) ≈ 56.736571404726476

  @test sum(deuc(links)) ≈ 2017.317

  @test sum(dtop(links)) ≈ 2763.8179999999998

  @test sum(dmax(links)) ≈ 1510.4999999999998
  
  @test count(hasxic.(links)) == 32

  @test sum(maxxic(links)) ≈ 8.5522559799e8

  @test sum(avgxic(links)) ≈ 2.4584517959989226e8

  @test sum(nspecies.(links)) == 245

  @test sum(count_nspecies.(links,mtol=20.0)) == 231

  @test count(consistency(links,tol=2.0)) == 26

  @test ismatch(resnames(links[2]),("SER","LYS"))
  @test ismatch(resnames(links[2]),("LYS","SER"))

  @test ismatch(indexes(links[2]),(99,131))
  @test ismatch(indexes(links[2]),(131,99))

  @test ismatch("D143-E108","D143-E108")
  @test ismatch("D108-E143","D143-E108")
  @test ismatch("D108-D143","D143-E108") == false
  @test ismatch("D108-D143","D100-E108") == false

  @test indomain(links[2],1:131)
  @test indomain(links[2],1:120) == false

  @test point_biserial(consistency(links),deuc(links)) ≈ -0.5720743658488596

end
