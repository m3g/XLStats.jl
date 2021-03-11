using XLStats
using Test

data_dir="./data"
xml_data="$data_dir/hitsDetail.xls"
links = read_all(xml_file=xml_data,
                 topolink_log="$data_dir/salbiii_topolink.log",
                 topolink_input="$data_dir/topolink.inp",
                 xic_file_name="$data_dir/salbiii_xic.dat",
                 domain=2:134)

@testset "XLStats.jl" begin

  @test name(links[1]) == "D3-E53"

  @test indexes(links[13]) == (6, 17)

  @test resnames(links[45]) == ("GLU","GLU")

  @test length(links) == 113

  @test sum(nscans(links)) == 761

  @test count(consistency(links)) == 26

  @test sum(skipmissing(maxscore1(links))) ≈ 297.84999999999997

  @test sum(skipmissing(maxscore2(links))) ≈ 86.866

  @test sum(skipmissing(avgscore1(links))) ≈ 120.55580639396922

  @test sum(skipmissing(avgscore2(links))) ≈ 30.502086812287494

  @test sum(deuc(links)) ≈ 2208.976

  @test sum(dtop(links)) ≈ 3008.168 

  @test sum(dmax(links)) ≈ 1680.1999999999998
  
  @test count(hasxic.(links)) == 20

  @test sum(skipmissing(maxxic(links))) ≈ 6.2237511644e8

  @test sum(avgxic(links)) ≈ 2.686064260940843e8

  @test sum(nspecies.(links)) == 159

  @test sum(count_nspecies.(links,mtol=20.0)) == 154

  @test count(consistency(links,tol=2.0)) == 29

  @test ismatch(resnames(links[2]),("ASP","GLU"))
  @test ismatch(resnames(links[2]),("GLU","ASP"))

  @test ismatch(indexes(links[2]),(3,111))
  @test ismatch(indexes(links[2]),(111,3))

  @test ismatch("D143-E108","D143-E108")
  @test ismatch("D108-E143","D143-E108")
  @test ismatch("D108-D143","D143-E108") == false
  @test ismatch("D108-D143","D100-E108") == false

  @test indomain(links[2],1:131)
  @test indomain(links[112],1:120) == false

  @test point_biserial(consistency(links),deuc(links)) ≈ -0.5950247343012209

end
