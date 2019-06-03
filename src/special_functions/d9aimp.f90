!** D9AIMP
SUBROUTINE D9AIMP(X,Ampl,Theta)
  !>
  !  Evaluate the Airy modulus and phase.
  !***
  ! **Library:**   SLATEC (FNLIB)
  !***
  ! **Category:**  C10D
  !***
  ! **Type:**      DOUBLE PRECISION (R9AIMP-S, D9AIMP-D)
  !***
  ! **Keywords:**  AIRY FUNCTION, FNLIB, MODULUS, PHASE, SPECIAL FUNCTIONS
  !***
  ! **Author:**  Fullerton, W., (LANL)
  !***
  ! **Description:**
  !
  ! Evaluate the Airy modulus and phase for X .LE. -1.0
  !
  ! Series for AM20       on the interval -1.56250E-02 to  0.
  !                                        with weighted error   3.12E-32
  !                                         log weighted error  31.51
  !                               significant figures required  29.24
  !                                    decimal places required  32.38
  !
  ! Series for ATH0       on the interval -1.56250E-02 to  0.
  !                                        with weighted error   2.75E-32
  !                                         log weighted error  31.56
  !                               significant figures required  30.17
  !                                    decimal places required  32.42
  !
  ! Series for AM21       on the interval -1.25000E-01 to -1.56250E-02
  !                                        with weighted error   3.40E-32
  !                                         log weighted error  31.47
  !                               significant figures required  29.02
  !                                    decimal places required  32.36
  !
  ! Series for ATH1       on the interval -1.25000E-01 to -1.56250E-02
  !                                        with weighted error   2.94E-32
  !                                         log weighted error  31.53
  !                               significant figures required  30.08
  !                                    decimal places required  32.41
  !
  ! Series for AM22       on the interval -1.00000E+00 to -1.25000E-01
  !                                        with weighted error   3.76E-32
  !                                         log weighted error  31.42
  !                               significant figures required  29.47
  !                                    decimal places required  32.36
  !
  ! Series for ATH2       on the interval -1.00000E+00 to -1.25000E-01
  !                                        with weighted error   4.97E-32
  !                                         log weighted error  31.30
  !                               significant figures required  29.79
  !                                    decimal places required  32.23
  !
  !***
  ! **References:**  (NONE)
  !***
  ! **Routines called:**  D1MACH, DCSEVL, INITDS, XERMSG

  !* REVISION HISTORY  (YYMMDD)
  !   770701  DATE WRITTEN
  !   890531  Changed all specific intrinsics to generic.  (WRB)
  !   890531  REVISION DATE from Version 3.2
  !   891214  Prologue converted to Version 4.0 format.  (BAB)
  !   900315  CALLs to XERROR changed to CALLs to XERMSG.  (THJ)
  !   900720  Routine changed from user-callable to subsidiary.  (WRB)
  USE service, ONLY : XERMSG, D1MACH
  REAL(8) :: X, Ampl, Theta
  REAL(8) :: sqrtx, z
  INTEGER, SAVE :: nam20, nath0, nam21, nath1, nam22, nath2
  REAL(8), PARAMETER :: eta = 0.1D0*D1MACH(3), xsml = -1.0D0/D1MACH(3)**0.3333D0
  REAL(8), PARAMETER :: am20cs(57) = [ +.108716749086561856615730588125D-1, &
    +.369489228982663555091728665146D-3, +.440680100484689563667507001327D-5, &
    +.143686762361911153929183952833D-6, +.824275552390078308670628855353D-8, &
    +.684426758893661606173927278180D-9, +.739566697282739287731004740213D-10, &
    +.974595633696825017638702600847D-11, +.150076885829405775650973119497D-11, &
    +.262147910221527634206252854802D-12, +.508354111376487180357278966914D-13, &
    +.107684753358811440492985997070D-13, +.246091286618433429335914062617D-14, &
    +.600786380358656418436110373550D-15, +.155449156102388071150651388384D-15, &
    +.423535125035576604426382780182D-16, +.120862166289299840154401109189D-16, &
    +.359609651214658240861499706423D-17, +.111134218386395638261774604677D-17, &
    +.355559532432366609893680289225D-18, +.117433021600139309998766947387D-18, &
    +.399397454661077561389162200966D-19, +.139576671528916310425606325640D-19, &
    +.500240055309236041393459280716D-20, +.183552760958132679184834866457D-20, &
    +.688490998179202743197790112404D-21, +.263631035611417012359996885105D-21, &
    +.102924890237338360287153563785D-21, +.409246966671594885489762960571D-22, &
    +.165558573406734651039727903828D-22, +.680797467063033356116599685727D-23, &
    +.284326559934079832419751134476D-23, +.120507398348965255097287818819D-23, &
    +.517961243287505217976613610424D-24, +.225622613427562816303268640887D-24, &
    +.995418801147745168832117078246D-25, +.444551696397342424308280582053D-25, &
    +.200865195461501101425916097338D-25, +.917786344151775165973885645402D-26, &
    +.423872958105589240661672197948D-26, +.197789272007846092370846251490D-26, &
    +.932116351284620665680435253373D-27, +.443482133249918099955611379722D-27, &
    +.212945672365573895594589552837D-27, +.103158569651075977552209344907D-27, &
    +.504023773022591199157904590029D-28, +.248301304570155945304046541005D-28, &
    +.123301783128562196054198238560D-28, +.617033449920521746121976730507D-29, &
    +.311092617415918897233869792213D-29, +.157983085201706173015269071503D-29, &
    +.807931987538283607678121339092D-30, +.415997394138667562722951360052D-30, &
    +.215610934097716900471935862504D-30, +.112468857265869178296752823613D-30, &
    +.590331560632838091123040811797D-31, +.311735667692928562046280505333D-31 ]
  REAL(8), PARAMETER :: ath0cs(53) = [ -.8172601764161634499840208700543D-1, &
    -.8004012824788273287596481113068D-3, -.3186525268782113203795553628242D-5, &
    -.6688388266477509330741698865033D-7, -.2931759284994564516506822463184D-8, &
    -.2011263760883621669049030307186D-9, -.1877522678055973426074008166652D-10, &
    -.2199637137704601251899002199848D-11, -.3071616682592272449025746605586D-12, &
    -.4936140553673418361025600985389D-13, -.8902833722583660416935236969866D-14, &
    -.1768987764615272613656814199467D-14, -.3817868689032277014678199609600D-15, &
    -.8851159014819947594156286509984D-16, -.2184818181414365953149677679568D-16, &
    -.5700849046986452380599442295119D-17, -.1563121122177875392516031795495D-17, &
    -.4481437996768995067906688776353D-18, -.1337794883736188022044566044098D-18, &
    -.4143340036874114453776852445442D-19, -.1327263385718805025080481164652D-19, &
    -.4385728589128440522215756835955D-20, -.1491360695952818067686201743956D-20, &
    -.5208104738630711377154238188773D-21, -.1864382222390498923872526604979D-21, &
    -.6830263751167969012975435381881D-22, -.2557117058029329629296207591347D-22, &
    -.9770158640254300218246907254046D-23, -.3805161433416679084068428254886D-23, &
    -.1509022750737054063493926482995D-23, -.6087551341242424929005568014525D-24, &
    -.2495879513809711495425982124058D-24, -.1039157654581920948909588084274D-24, &
    -.4390235913976846536974594969051D-25, -.1880790678447990211675826820582D-25, &
    -.8165070764199462948863022205753D-26, -.3589944503749750514266435585041D-26, &
    -.1597658126632132872981291608708D-26, -.7193250175703823969113802835305D-27, &
    -.3274943012727856506209351132721D-27, -.1507042445783690665816975047272D-27, &
    -.7006624198319904717843967949140D-28, -.3289907402983718226528815678356D-28, &
    -.1559518084365146526445322711496D-28, -.7460690508208254582833851119721D-29, &
    -.3600877034824662020563277249431D-29, -.1752851437473772257350402219197D-29, &
    -.8603275775188512909623778628724D-30, -.4256432603226946534668039480105D-30, &
    -.2122161865044262927723650698206D-30, -.1065996156704879052472060798561D-30, &
    -.5393568608816949116410688086892D-31, -.2748174851043954822278496517870D-31 ]
  REAL(8), PARAMETER :: am21cs(60) = [ +.592790266721309588375717482814D-2, &
    +.200569405393165186428695217690D-2, +.911081850262275893553072526291D-4, &
    +.849894306372047155633172107475D-5, +.113297908976913076637929215494D-5, &
    +.187517946100666496180950627804D-6, +.359306519018245832699035211192D-7, &
    +.765757714071683864039093517470D-8, +.176999967168039173925953460744D-8, &
    +.436259555654598932720546585535D-9, +.113291641337853230035520085219D-9, &
    +.307257690982419244137868398126D-10, +.864482416482201075541200465766D-11, &
    +.251015250060924402115104562212D-11, +.749102496764440371601802227751D-12, &
    +.228996928487994073089565214432D-12, +.715113658927987694949327491175D-13, &
    +.227607924959566841946395165061D-13, +.736942142760886513969953227782D-14, &
    +.242328675267827490463991742006D-14, +.808153774548239869283406558403D-15, &
    +.273008079804356086659174563386D-15, +.933236070891385318473519474326D-16, &
    +.322508099681084622213867546973D-16, +.112581932346444541217757573416D-16, &
    +.396699463986938821660259459530D-17, +.141006567944319504660865034527D-17, &
    +.505302086537851213375537393032D-18, +.182461523215945141197999102789D-18, &
    +.663584568262130466928029121642D-19, +.242963731631276179741747455826D-19, &
    +.895238915123687802013669922963D-20, +.331845289350050791260229250755D-20, &
    +.123706196188658315384437905922D-20, +.463636677012390840306767734243D-21, &
    +.174653135947764475469758765989D-21, +.661116810234991176307910643111D-22, &
    +.251409918994072486176125666459D-22, +.960274995571732568694034386998D-23, &
    +.368324952289296395686436898078D-23, +.141843138269159136145535939553D-23, &
    +.548342674276935830106345800990D-24, +.212761054623118806650372562616D-24, &
    +.828443700849418591487734760953D-25, +.323670563926127001421028600927D-25, &
    +.126868882963286057355055062493D-25, +.498843818992121626935068934362D-26, &
    +.196734584467649390967119381790D-26, +.778135971020326957713212064836D-27, &
    +.308633941498911152919192968451D-27, +.122744647045453119789338037234D-27, &
    +.489431279134292205885241216204D-28, +.195646879802909821175925099724D-28, &
    +.783988952922426171166311492266D-29, +.314896914002484223748298978099D-29, &
    +.126769763137250681307067842559D-29, +.511470691906900141641632107724D-30, &
    +.206801709795538770250900316706D-30, +.837891344768519001325996867583D-31, &
    +.340168991971489802052339079577D-31 ]
  REAL(8), PARAMETER :: ath1cs(58) = [ -.6972849916208883845888148415037D-1, &
    -.5108722790650044987073448077961D-2, -.8644335996989755094525334749512D-4, &
    -.5604720044235263542188698916125D-5, -.6045735125623897409156376640077D-6, &
    -.8639802632488334393219721138499D-7, -.1480809484309927157147782480780D-7, &
    -.2885809334577236039999449908712D-8, -.6191631975665699609309191231800D-9, &
    -.1431992808860957830931365259879D-9, -.3518141102137214721504616874321D-10, &
    -.9084761919955078290070339808051D-11, -.2446171672688598449343283664767D-11, &
    -.6826083203213446240828996710264D-12, -.1964579931194940171278546257802D-12, &
    -.5808933227139693164009191265856D-13, -.1759042249527441992795400959024D-13, &
    -.5440902932714896613632538945319D-14, -.1715247407486806802622358519451D-14, &
    -.5500929233576991546871101847161D-15, -.1791878287739317259495152638754D-15, &
    -.5920372520086694197778411062231D-16, -.1981713027876483962470972206590D-16, &
    -.6713232347016352262049984343790D-17, -.2299450243658281116122358619832D-17, &
    -.7957300928236376595304637145634D-18, -.2779994027291784157172290233739D-18, &
    -.9798924361326985224406795480814D-19, -.3482717006061574386702645565849D-19, &
    -.1247489122558599057173300058084D-19, -.4501210041478228113487751824452D-20, &
    -.1635346244013352135596114164667D-20, -.5980102897780336268098762265941D-21, &
    -.2200246286286123454028196295475D-21, -.8142463073515085897408205291519D-22, &
    -.3029924773660042537432330709674D-22, -.1133390098574623537722943969689D-22, &
    -.4260766024749295719283049889791D-23, -.1609363396278189718797500634453D-23, &
    -.6106377190825026293045330444287D-24, -.2326954318021694061836577887573D-24, &
    -.8903987877472252604474129558186D-25, -.3420558530005675024117914752341D-25, &
    -.1319026715257272659017212100607D-25, -.5104899493612043091316191177386D-26, &
    -.1982599478474547451242444663466D-26, -.7725702356880830535636111851519D-27, &
    -.3020234733664680100815776863573D-27, -.1184379739074169993712946380800D-27, &
    -.4658430227922308520573252840106D-28, -.1837554188100384647157502006613D-28, &
    -.7268566894427990953321876684800D-29, -.2882863120391468135527089875626D-29, &
    -.1146374629459906350417591664639D-29, -.4570031437748533058179991688533D-30, &
    -.1826276602045346104809934028799D-30, -.7315349993385250469111066350933D-31, &
    -.2936925599971429781637815773866D-31 ]
  REAL(8), PARAMETER :: am22cs(74) = [ -.156284448062534112753545828583D-1, &
    +.778336445239681307018943100334D-2, +.867057770477189528406072812110D-3, &
    +.156966273156113719469953482266D-3, +.356396257143286511324100666302D-4, &
    +.924598335425043154495080090994D-5, +.262110161850422389523194982066D-5, &
    +.791882216516012561489469982263D-6, +.251041527921011847803162690862D-6, &
    +.826522320665407734472997712940D-7, +.280571166281305264396384290014D-7, &
    +.976821090484680786674631273890D-8, +.347407923227710343287279035573D-8, &
    +.125828132169836914219092738164D-8, +.462988260641895264497330784625D-9, &
    +.172728258813604072468143128696D-9, +.652319200131154135148574124970D-10, &
    +.249047168520982056019881087112D-10, +.960156820553765948078189890126D-11, &
    +.373448002067726856974776596757D-11, +.146417565032053391722216189678D-11, &
    +.578265471168512825475827881553D-12, +.229915407244706118560254184494D-12, &
    +.919780711231997257150883662365D-13, +.370060068813090065807504045556D-13, &
    +.149675761698672987823326345205D-13, +.608361194938461148720451399443D-14, &
    +.248404087115121397635425326873D-14, +.101862476526769080727914465339D-14, &
    +.419383856352753989429640310957D-15, +.173318901762930756149702493501D-15, &
    +.718821902388508517820445406811D-16, +.299123633598403607712470896113D-16, &
    +.124868990433238627855713110880D-16, +.522829344609483661928651193632D-17, &
    +.219532961724713396595998454359D-17, +.924298325229777281154410024332D-18, &
    +.390157708236091407825543197309D-18, +.165093892693863707213759030367D-18, &
    +.700221815715994367565716554487D-19, +.297651833616786915573214963506D-19, &
    +.126796539086902072571134261229D-19, +.541243400697077628687581725061D-20, &
    +.231487350218155252296382133283D-20, +.991920288386566563462623851167D-21, &
    +.425803015323732357158897608174D-21, +.183101842973024501678402003088D-21, &
    +.788678712311075375564526811022D-22, +.340254607386229874956582997235D-22, &
    +.147020881405712530791860892535D-22, +.636211018324916957733348071767D-23, &
    +.275707050680980721919395987768D-23, +.119645858090104071356261780457D-23, &
    +.519912545729242147981768210567D-24, +.226217674847104475260575286850D-24, &
    +.985526113754431819448565068283D-25, +.429870630332508717223681286187D-25, &
    +.187723641661580639829657670189D-25, +.820721941772842137268801052115D-26, &
    +.359214665604615507812767944463D-26, +.157390594612773315611458940587D-26, &
    +.690329781039333834965319153586D-27, +.303092079078968534607859331415D-27, &
    +.133204934160481219185689121944D-27, +.585978836851523490117937981442D-28, &
    +.258016868489487806338425080457D-28, +.113712433637283667223632182863D-28, &
    +.501592557226068509236430548549D-29, +.221445829395509373322569708484D-29, &
    +.978470283886507289984691416411D-30, +.432695414934180170112000952983D-30, &
    +.191497288193994570612929860440D-30, +.848164622402392354171298331562D-31, &
    +.375947065173955919947455052934D-31 ]
  REAL(8), PARAMETER :: ath2cs(72) = [ +.4405273458718778997061127057775D-2, &
    -.3042919452318454608483844239873D-1, -.1385653283771793791602692842653D-2, &
    -.1804443908954952302670486910952D-3, -.3380847108327308671057465323618D-4, &
    -.7678183535229023055257676817765D-5, -.1967839443716035324690935417077D-5, &
    -.5483727115877700361586143659281D-6, -.1625461550532612452712696212258D-6, &
    -.5053049981268895015277637842078D-7, -.1631580701124066881183851715617D-7, &
    -.5434204112348517507963436694817D-8, -.1857398556409900325763850109630D-8, &
    -.6489512033326108816213513640676D-9, -.2310594885800944720482995987079D-9, &
    -.8363282183204411682819329546745D-10, -.3071196844890191462660661303891D-10, &
    -.1142367142432716819409514579892D-10, -.4298116066345803065822470108971D-11, &
    -.1633898699596715440601646086632D-11, -.6269328620016619432123443754076D-12, &
    -.2426052694816257357356159203991D-12, -.9461198321624039090742527765052D-13, &
    -.3716060313411504806847798281269D-13, -.1469155684097526763170138810309D-13, &
    -.5843694726140911944556401363094D-14, -.2337502595591951298832675034934D-14, &
    -.9399231371171435401160167358411D-15, -.3798014669372894500076335263715D-15, &
    -.1541731043984972524883443681775D-15, -.6285287079535307162925662365202D-16, &
    -.2572731812811455424755383992774D-16, -.1057098119354017809340974866555D-16, &
    -.4359080267402696966695992699964D-17, -.1803634315959978013953176945540D-17, &
    -.7486838064380536821719431676914D-18, -.3117261367347604656799597209985D-18, &
    -.1301687980927700734792871620696D-18, -.5450527587519522468973883909909D-19, &
    -.2288293490114231872268635931903D-19, -.9631059503829538655655060440088D-20, &
    -.4063281001524614089092195416434D-20, -.1718203980908026763900413858510D-20, &
    -.7281574619892536367415322473328D-21, -.3092352652680643127960680345790D-21, &
    -.1315917855965440490383417023254D-21, -.5610606786087055512664907412668D-22, &
    -.2396621894086355206020304337895D-22, -.1025574332390581200832954423924D-22, &
    -.4396264138143656476403607323663D-23, -.1887652998372577373342508719450D-23, &
    -.8118140359576807603579433230445D-24, -.3496734274366286856375952089214D-24, &
    -.1508402925156873215171751475867D-24, -.6516268284778671059787773834341D-25, &
    -.2818945797529207424505942114583D-25, -.1221127596512262744598094464505D-25, &
    -.5296674341169867168620011705073D-26, -.2300359270773673431358870971744D-26, &
    -.1000279482355367494781220348930D-26, -.4354760404180879394806893162179D-27, &
    -.1898056134741477522515482827030D-27, -.8282111868712974697554009309315D-28, &
    -.3617815493066569006586213484374D-28, -.1582018896178003654858941843636D-28, &
    -.6925068597802270011772820383247D-29, -.3034390239778629128908629727335D-29, &
    -.1330889568166725224761977446509D-29, -.5842848522173090120487606971706D-30, &
    -.2567488423238302631121274357678D-30, -.1129232322268882185791505819151D-30, &
    -.4970947029753336916550570105023D-31 ]
  REAL(8), PARAMETER ::  pi4 = 0.78539816339744830961566084581988D0
  LOGICAL, SAVE :: first = .TRUE.
  !* FIRST EXECUTABLE STATEMENT  D9AIMP
  IF ( first ) THEN
    nam20 = INITDS(am20cs,57,eta)
    nath0 = INITDS(ath0cs,53,eta)
    nam21 = INITDS(am21cs,60,eta)
    nath1 = INITDS(ath1cs,58,eta)
    nam22 = INITDS(am22cs,74,eta)
    nath2 = INITDS(ath2cs,72,eta)
    first = .FALSE.
  END IF
  !
  IF ( X<(-4.0D0) ) THEN
    z = 1.0D0
    IF ( X>xsml ) z = 128.D0/X**3 + 1.0D0
    Ampl = 0.3125D0 + DCSEVL(z,am20cs,nam20)
    Theta = -0.625D0 + DCSEVL(z,ath0cs,nath0)
    !
  ELSEIF ( X>=(-2.0D0) ) THEN
    !
    IF ( X>=(-1.0D0) ) CALL XERMSG('D9AIMP','X MUST BE LE -1.0',1,2)
    !
    z = (16.D0/X**3+9.0D0)/7.0D0
    Ampl = 0.3125D0 + DCSEVL(z,am22cs,nam22)
    Theta = -0.625D0 + DCSEVL(z,ath2cs,nath2)
  ELSE
    z = (128.D0/X**3+9.0D0)/7.0D0
    Ampl = 0.3125D0 + DCSEVL(z,am21cs,nam21)
    Theta = -0.625D0 + DCSEVL(z,ath1cs,nath1)
  END IF
  !
  sqrtx = SQRT(-X)
  Ampl = SQRT(Ampl/sqrtx)
  Theta = pi4 - X*sqrtx*Theta
  !
END SUBROUTINE D9AIMP
