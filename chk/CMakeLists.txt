FILE(GLOB tests test*.f90)

add_executable(test01 test01.f90 qc6a.f90)

add_executable(test02 test02.f90 sfnck.f90)

add_executable(test03 test03.f90 dfnck.f90)

add_executable(test04 test04.f90 cfnck.f90)

add_executable(test05 test05.f90 eg8ck.f90 bikck.f90 bjyck.f90 fein.f90)

add_executable(test06 test06.f90 deg8ck.f90 dbikck.f90 dbjyck.f90 dfein.f90)

add_executable(test07 test07.f90 qckin.f90 qcpsi.f90)

add_executable(test08 test08.f90 dqckin.f90 dqcpsi.f90)

add_executable(test09 test09.f90 cqcai.f90 cqcbh.f90 cqcbj.f90 cqcbk.f90 cqcby.f90 cqcbi.f90)

add_executable(test10 test10.f90 zqcai.f90 zqcbh.f90 zqcbj.f90 zqcbk.f90 zqcby.f90 zqcbi.f90)

add_executable(test11 test11.f90 fcnqx1.f90 xcsrt.f90)

add_executable(test12 test12.f90 fcnqx2.f90 dxcsrt.f90)

add_executable(test13 test13.f90 qcrc.f90 qcrd.f90 qcrf.f90 qcrj.f90)

add_executable(test14 test14.f90 qcdrc.f90 qcdrd.f90 qcdrf.f90 qcdrj.f90)

add_executable(test15 test15.f90 qc36j.f90)

add_executable(test16 test16.f90 dqc36j.f90)

add_executable(test17 test17.f90 blachk.f90 header.f90 check1.f90 check2.f90 check0.f90 itest.f90 stest.f90 dtest.f90)

add_executable(test18 test18.f90 sblat2.f90 sblat3.f90 smvch.f90 lse.f90 schke2.f90 schk12.f90 schk22.f90 schk32.f90 schk42.f90 schk52.f90 schk62.f90 smmch.f90 schke3.f90 schk13.f90 schk23.f90 schk33.f90 schk43.f90 schk53.f90 chkxer.f90 smake2.f90 lseres.f90 smake3.f90 sbeg.f90)

add_executable(test19 test19.f90 dblat2.f90 dblat3.f90 dmvch.f90 lde.f90 dchke2.f90 dchk12.f90 dchk22.f90 dchk32.f90 dchk42.f90 dchk52.f90 dchk62.f90 dmmch.f90 dchke3.f90 dchk13.f90 dchk23.f90 dchk33.f90 dchk43.f90 dchk53.f90 chkxer.f90 dmake2.f90 lderes.f90 dmake3.f90 dbeg.f90)

add_executable(test20 test20.f90 cblat2.f90 cblat3.f90 cmvch.f90 lce.f90 cchke2.f90 cchk12.f90 cchk22.f90 cchk32.f90 cchk42.f90 cchk52.f90 cchk62.f90 cmmch.f90 cchke3.f90 cchk13.f90 cchk23.f90 cchk33.f90 cchk43.f90 cchk53.f90 chkxer.f90 cmake2.f90 lceres.f90 cmake3.f90 cbeg.f90)

add_executable(test21 test21.f90 sgeqc.f90 dgeqc.f90 cgeqc.f90 delx.f90 )

add_executable(test22 test22.f90 sqck.f90 dqck.f90 cqck.f90)

add_executable(test23 test23.f90 cgeck.f90 cgbqc.f90 cpoqc.f90 cppqc.f90 cpbqc.f90 csiqc.f90 cspqc.f90 chiqc.f90 chpqc.f90 ctrqc.f90 cgtqc.f90 cptqc.f90 cchqc.f90 cqrqc.f90 csvqc.f90 delx.f90 )

add_executable(test24 test24.f90 eisqx1.f90 eisqx2.f90)

add_executable(test25 test25.f90 slapqc.f90 srmgen.f90 vfill.f90 outerr.f90 ismpl.f90)

add_executable(test26 test26.f90 dlapqc.f90 drmgen.f90 dfill.f90 duterr.f90 ismpl.f90)

add_executable(test27 test27.f90 lseiqx.f90 qcglss.f90)

add_executable(test28 test28.f90 dlseit.f90 dqcgls.f90)

add_executable(test29 test29.f90 pntchk.f90 dpntck.f90)

add_executable(test30 test30.f90 bspck.f90 fb.f90)

add_executable(test31 test31.f90 dbspck.f90 dfb.f90)

add_executable(test32 test32.f90 pchqk1.f90 pchqk2.f90 pchqk3.f90 pchqk4.f90 pchqk5.f90 comp.f90 evchck.f90 evpcck.f90 everck.f90 fdtrue.f90)

add_executable(test33 test33.f90 dpchq1.f90 dpchq2.f90 dpchq3.f90 dpchq4.f90 dpchq5.f90 comp.f90 devchk.f90 devpck.f90 deverk.f90 dfdtru.f90)

add_executable(test34 test34.f90 cprpqx.f90 fztest.f90 dfztst.f90 rqrtst.f90 cqrtst.f90 pass.f90)

add_executable(test35 test35.f90 snsqqk.f90 sosnqx.f90 pass.f90 sosfnc.f90 sqfcn2.f90 sqjac2.f90)

add_executable(test36 test36.f90 dnsqqk.f90 dsosqx.f90 pass.f90 dsosfn.f90 dqfcn2.f90 dqjac2.f90)

add_executable(test37 test37.f90 splpqx.f90 sbocqx.f90 pass.f90)

add_executable(test38 test38.f90 dplpqx.f90 dbocqx.f90 pass.f90)

add_executable(test39 test39.f90 cqng.f90 cqag.f90 cqags.f90 cqagi.f90 cqagp.f90 cqawo.f90 cqawf.f90 cqaws.f90 cqawc.f90 cprin.f90 f1n.f90 f2n.f90 f1g.f90 f2g.f90 f3g.f90 f0s.f90 f1s.f90 f2s.f90 f3s.f90 f4s.f90 f5s.f90 t0.f90 t1.f90 t2.f90 t3.f90 t4.f90 t5.f90 f1p.f90 f2p.f90 f3p.f90 f4p.f90 f0o.f90 f1o.f90 f2o.f90 f0f.f90 f1f.f90 f0ws.f90 f1ws.f90 f0c.f90 f1c.f90)

add_executable(test40 test40.f90 cdqng.f90 cdqag.f90 cdqags.f90 cdqagi.f90 cdqagp.f90 cdqawo.f90 cdqawf.f90 cdqaws.f90 cdqawc.f90 dprin.f90 df1n.f90 df2n.f90 df1g.f90 df2g.f90 df3g.f90 df0s.f90 df1s.f90 df2s.f90 df3s.f90 df4s.f90 df5s.f90 dt0.f90 dt1.f90 dt2.f90 dt3.f90 dt4.f90 dt5.f90 df1p.f90 df2p.f90 df3p.f90 df4p.f90 df0o.f90 df1o.f90 df2o.f90 df0f.f90 df1f.f90 df0ws.f90 df1ws.f90 df0c.f90 df1c.f90)

add_executable(test41 test41.f90 avntst.f90 qg8tst.f90 qn79qx.f90 fqd1.f90 fqd2.f90)

add_executable(test42 test42.f90 davnts.f90 dqg8ts.f90 dqn79q.f90 dfqd1.f90 dfqd2.f90)

add_executable(test43 test43.f90 qxabm.f90 qxbdf.f90 qxrkf.f90 qxbvsp.f90 pass.f90 jac.f90 fdeqc.f90)

add_executable(test44 test44.f90 qxdabm.f90 qxdbdf.f90 qxdrkf.f90 qxdbvs.f90 pass.f90 djac.f90 dfdeqc.f90)

add_executable(test45 test45.f90 sdqck.f90 sdf.f90)

add_executable(test46 test46.f90 ddqck.f90 ddf.f90)

add_executable(test47 test47.f90 cdqck.f90 cdf.f90)

add_executable(test48 test48.f90 sdasqc.f90 edit2.f90 sdjac1.f90 sdjac2.f90 sdres1.f90 sdres2.f90)

add_executable(test49 test49.f90 ddasqc.f90 dedit2.f90 ddjac1.f90 ddjac2.f90 ddres1.f90 ddres2.f90)

add_executable(test50 test50.f90 qxplr.f90 qxcrt.f90 qxcyl.f90 qxssp.f90 qxcsp.f90 qxgbun.f90 qxblkt.f90)

add_executable(test51 test51.f90 fftqx.f90)

add_executable(test52 test52.f90 fcqx.f90 snls1q.f90 pfitqx.f90 fcn1.f90 fcn2.f90 fcn3.f90 pass.f90 cmpare.f90)

add_executable(test53 test53.f90 dfcqx.f90 dnls1q.f90 dpfitt.f90 dfcn1.f90 dfcn2.f90 dfcn3.f90 pass.f90 dcmpar.f90)

add_executable(test54 test54.f90 isrtqc.f90 ssrtqc.f90 dsrtqc.f90 hsrtqc.f90)

foreach( tst ${tests} )
  STRING(REPLACE ".f90" "" exe1 ${tst})
  STRING(REPLACE "${CMAKE_CURRENT_SOURCE_DIR}/" "" exe ${exe1})
  target_link_libraries(${exe} slatec)
  if(UNIX)
    add_test(NAME ${exe} COMMAND sh -c "./${exe} < ${CMAKE_CURRENT_SOURCE_DIR}/kprint.txt" WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
  elseif(WIN32)
    add_test(NAME ${exe} COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_BINARY_DIR} $ENV{ComSpec} /c "${exe} < ${CMAKE_CURRENT_SOURCE_DIR}/kprint.txt")
  endif()
  set_tests_properties(${exe} PROPERTIES FAIL_REGULAR_EXPRESSION "[^a-z]Error;ERROR;Failed;FAILED")
endforeach()
