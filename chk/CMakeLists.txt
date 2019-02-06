FILE(GLOB tests test*.f)

add_executable(test01 test01.f qc6a.f)

add_executable(test02 test02.f sfnck.f)

add_executable(test03 test03.f dfnck.f)

add_executable(test04 test04.f cfnck.f)

add_executable(test05 test05.f eg8ck.f bikck.f bjyck.f fein.f)

add_executable(test06 test06.f deg8ck.f dbikck.f dbjyck.f dfein.f)

add_executable(test07 test07.f qckin.f qcpsi.f)

add_executable(test08 test08.f dqckin.f dqcpsi.f)

add_executable(test09 test09.f cqcai.f cqcbh.f cqcbj.f cqcbk.f cqcby.f cqcbi.f)

add_executable(test10 test10.f zqcai.f zqcbh.f zqcbj.f zqcbk.f zqcby.f zqcbi.f)

add_executable(test11 test11.f fcnqx1.f xcsrt.f)

add_executable(test12 test12.f fcnqx2.f dxcsrt.f)

add_executable(test13 test13.f qcrc.f qcrd.f qcrf.f qcrj.f)

add_executable(test14 test14.f qcdrc.f qcdrd.f qcdrf.f qcdrj.f)

add_executable(test15 test15.f qc36j.f)

add_executable(test16 test16.f dqc36j.f)

add_executable(test17 test17.f blachk.f header.f check1.f check2.f check0.f itest.f stest.f dtest.f)

add_executable(test18 test18.f sblat2.f sblat3.f smvch.f lse.f schke2.f schk12.f schk22.f schk32.f schk42.f schk52.f schk62.f smmch.f schke3.f schk13.f schk23.f schk33.f schk43.f schk53.f chkxer.f smake2.f lseres.f smake3.f sbeg.f)

add_executable(test19 test19.f dblat2.f dblat3.f dmvch.f lde.f dchke2.f dchk12.f dchk22.f dchk32.f dchk42.f dchk52.f dchk62.f dmmch.f dchke3.f dchk13.f dchk23.f dchk33.f dchk43.f dchk53.f chkxer.f dmake2.f lderes.f dmake3.f dbeg.f)

add_executable(test20 test20.f cblat2.f cblat3.f cmvch.f lce.f cchke2.f cchk12.f cchk22.f cchk32.f cchk42.f cchk52.f cchk62.f cmmch.f cchke3.f cchk13.f cchk23.f cchk33.f cchk43.f cchk53.f chkxer.f cmake2.f lceres.f cmake3.f cbeg.f)

add_executable(test21 test21.f sgeqc.f dgeqc.f cgeqc.f )

add_executable(test22 test22.f sqck.f dqck.f cqck.f)

add_executable(test23 test23.f cgeck.f cgbqc.f cpoqc.f cppqc.f cpbqc.f csiqc.f cspqc.f chiqc.f chpqc.f ctrqc.f cgtqc.f cptqc.f cchqc.f cqrqc.f csvqc.f )

add_executable(test24 test24.f eisqx1.f eisqx2.f)

add_executable(test25 test25.f slapqc.f srmgen.f vfill.f outerr.f ismpl.f)

add_executable(test26 test26.f dlapqc.f drmgen.f dfill.f duterr.f ismpl.f)

add_executable(test27 test27.f lseiqx.f qcglss.f)

add_executable(test28 test28.f dlseit.f dqcgls.f)

add_executable(test29 test29.f pntchk.f dpntck.f)

add_executable(test30 test30.f bspck.f fb.f)

add_executable(test31 test31.f dbspck.f dfb.f)

add_executable(test32 test32.f pchqk1.f pchqk2.f pchqk3.f pchqk4.f pchqk5.f comp.f evchck.f evpcck.f everck.f fdtrue.f)

add_executable(test33 test33.f dpchq1.f dpchq2.f dpchq3.f dpchq4.f dpchq5.f comp.f devchk.f devpck.f deverk.f dfdtru.f)

add_executable(test34 test34.f cprpqx.f fztest.f dfztst.f rqrtst.f cqrtst.f pass.f)

add_executable(test35 test35.f snsqqk.f sosnqx.f pass.f sosfnc.f sqfcn2.f sqjac2.f)

add_executable(test36 test36.f dnsqqk.f dsosqx.f pass.f dsosfn.f dqfcn2.f dqjac2.f)

add_executable(test37 test37.f splpqx.f sbocqx.f pass.f)

add_executable(test38 test38.f dplpqx.f dbocqx.f pass.f)

add_executable(test39 test39.f cqng.f cqag.f cqags.f cqagi.f cqagp.f cqawo.f cqawf.f cqaws.f cqawc.f cprin.f f1n.f f2n.f f1g.f f2g.f f3g.f f0s.f f1s.f f2s.f f3s.f f4s.f f5s.f t0.f t1.f t2.f t3.f t4.f t5.f f1p.f f2p.f f3p.f f4p.f f0o.f f1o.f f2o.f f0f.f f1f.f f0ws.f f1ws.f f0c.f f1c.f)

add_executable(test40 test40.f cdqng.f cdqag.f cdqags.f cdqagi.f cdqagp.f cdqawo.f cdqawf.f cdqaws.f cdqawc.f dprin.f df1n.f df2n.f df1g.f df2g.f df3g.f df0s.f df1s.f df2s.f df3s.f df4s.f df5s.f dt0.f dt1.f dt2.f dt3.f dt4.f dt5.f df1p.f df2p.f df3p.f df4p.f df0o.f df1o.f df2o.f df0f.f df1f.f df0ws.f df1ws.f df0c.f df1c.f)

add_executable(test41 test41.f avntst.f qg8tst.f qn79qx.f fqd1.f fqd2.f)

add_executable(test42 test42.f davnts.f dqg8ts.f dqn79q.f dfqd1.f dfqd2.f)

add_executable(test43 test43.f qxabm.f qxbdf.f qxrkf.f qxbvsp.f pass.f jac.f fdeqc.f)

add_executable(test44 test44.f qxdabm.f qxdbdf.f qxdrkf.f qxdbvs.f pass.f djac.f dfdeqc.f)

add_executable(test45 test45.f sdqck.f sdf.f)

add_executable(test46 test46.f ddqck.f ddf.f)

add_executable(test47 test47.f cdqck.f cdf.f)

add_executable(test48 test48.f sdasqc.f edit2.f sdjac1.f sdjac2.f sdres1.f sdres2.f)

add_executable(test49 test49.f ddasqc.f dedit2.f ddjac1.f ddjac2.f ddres1.f ddres2.f)

add_executable(test50 test50.f qxplr.f qxcrt.f qxcyl.f qxssp.f qxcsp.f qxgbun.f qxblkt.f)

add_executable(test51 test51.f fftqx.f)

add_executable(test52 test52.f fcqx.f snls1q.f pfitqx.f fcn1.f fcn2.f fcn3.f pass.f cmpare.f)

add_executable(test53 test53.f dfcqx.f dnls1q.f dpfitt.f dfcn1.f dfcn2.f dfcn3.f pass.f dcmpar.f)

add_executable(test54 test54.f isrtqc.f ssrtqc.f dsrtqc.f hsrtqc.f)

foreach( tst ${tests} )
  STRING(REPLACE ".f" "" exe1 ${tst})
  STRING(REPLACE "${CMAKE_CURRENT_SOURCE_DIR}/" "" exe ${exe1})
  target_link_libraries(${exe} slatec)
  if(UNIX)
    add_test(NAME ${exe} COMMAND sh -c "./${exe} < ${CMAKE_CURRENT_SOURCE_DIR}/kprint.txt" WORKING_DIRECTORY ${CMAKE_BINARY_DIR})
  elseif(WIN32)
    add_test(NAME ${exe} COMMAND ${CMAKE_COMMAND} -E chdir ${CMAKE_BINARY_DIR} $ENV{ComSpec} /c "${exe} < ${CMAKE_CURRENT_SOURCE_DIR}/kprint.txt")
  endif()
  set_tests_properties(${exe} PROPERTIES FAIL_REGULAR_EXPRESSION "[^a-z]Error;ERROR;Failed;FAILED")
endforeach()
