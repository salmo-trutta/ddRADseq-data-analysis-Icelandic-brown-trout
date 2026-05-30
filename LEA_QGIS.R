##USING THE PACKAGE LEA##
library(LEA)

vcf2geno("iceland.recode.vcf", "iceland.geno")
pc = pca("iceland.geno", scale = TRUE)

tw = tracy.widom(pc)
tw$pvalues[1:5]

plot(tw$percentage)

#project<-load.snmfProject("iceland.snmfProject")
# it seems that K=17 given the knee at K=18 but we need more evidence
project = NULL
project = snmf("iceland.geno",
               K = 15:23,
               entropy = TRUE,
               repetitions = 2500,
               seed = 1236,
               CPU = 4,
               alpha = 100,
               project = "new")

par(family="Times")
plot(project, col = "blue", pch = 19, cex = 1.2)

## REPLACING NUMBER WITH SAMPLE ID

equivalencias <- c('1'	=	"Bak001"	,
                   '2'	=	"Bak003"	,
                   '3'	=	"Bak004"	,
                   '4'	=	"Bak006"	,
                   '5'	=	"Bak007"	,
                   '6'	=	"Bak008"	,
                   '7'	=	"Bak009"	,
                   '8'	=	"Bak010"	,
                   '9'	=	"Bak011"	,
                   '10'	=	"Duf002"	,
                   '11'	=	"Duf006"	,
                   '12'	=	"Duf007"	,
                   '13'	=	"Duf008"	,
                   '14'	=	"Duf012"	,
                   '15'	=	"Duf014"	,
                   '16'	=	"Duf015"	,
                   '17'	=	"Duf020"	,
                   '18'	=	"Efr001"	,
                   '19'	=	"Efr006"	,
                   '20'	=	"Efr009"	,
                   '21'	=	"Efr010"	,
                   '22'	=	"Efr011"	,
                   '23'	=	"For002"	,
                   '24'	=	"For007"	,
                   '25'	=	"For008"	,
                   '26'	=	"For010"	,
                   '27'	=	"For013"	,
                   '28'	=	"For014"	,
                   '29'	=	"For016"	,
                   '30'	=	"For017"	,
                   '31'	=	"For018"	,
                   '32'	=	"For020"	,
                   '33'	=	"For021"	,
                   '34'	=	"For023"	,
                   '35'	=	"For024"	,
                   '36'	=	"For025"	,
                   '37'	=	"For030"	,
                   '38'	=	"For031"	,
                   '39'	=	"For032"	,
                   '40'	=	"Fos001"	,
                   '41'	=	"Fos002"	,
                   '42'	=	"Fos003"	,
                   '43'	=	"Fos004"	,
                   '44'	=	"Fos005"	,
                   '45'	=	"Fos006"	,
                   '46'	=	"Fos007"	,
                   '47'	=	"Fos015"	,
                   '48'	=	"Fos009"	,
                   '49'	=	"Fos010"	,
                   '50'	=	"Fos011"	,
                   '51'	=	"Fos016"	,
                   '52'	=	"Fos017"	,
                   '53'	=	"Fos013"	,
                   '54'	=	"Fos014"	,
                   '55'	=	"Fre001"	,
                   '56'	=	"Fre002"	,
                   '57'	=	"Fre004"	,
                   '58'	=	"Fre005"	,
                   '59'	=	"Fre006"	,
                   '60'	=	"Fre007"	,
                   '61'	=	"Fre009"	,
                   '62'	=	"Fre010"	,
                   '63'	=	"Fre011"	,
                   '64'	=	"Fre012"	,
                   '65'	=	"Fre013"	,
                   '66'	=	"Fre014"	,
                   '67'	=	"Fre015"	,
                   '68'	=	"Fre016"	,
                   '69'	=	"Gra013"	,
                   '70'	=	"Gra014"	,
                   '71'	=	"Gra015"	,
                   '72'	=	"Gra017"	,
                   '73'	=	"Gra018"	,
                   '74'	=	"Gra019"	,
                   '75'	=	"Gra020"	,
                   '76'	=	"Gra022"	,
                   '77'	=	"Gra023"	,
                   '78'	=	"Gra025"	,
                   '79'	=	"Gra026"	,
                   '80'	=	"Gra027"	,
                   '81'	=	"Gra028"	,
                   '82'	=	"Gra029"	,
                   '83'	=	"Gra031"	,
                   '84'	=	"Gra032"	,
                   '85'	=	"Gre013"	,
                   '86'	=	"Gre014"	,
                   '87'	=	"Gre016"	,
                   '88'	=	"Gre017"	,
                   '89'	=	"Gre018"	,
                   '90'	=	"Gre071"	,
                   '91'	=	"Gre072"	,
                   '92'	=	"Gre073"	,
                   '93'	=	"Gre074"	,
                   '94'	=	"Gre089"	,
                   '95'	=	"Gre090"	,
                   '96'	=	"Gre091"	,
                   '97'	=	"Gre092"	,
                   '98'	=	"Gre075"	,
                   '99'	=	"Gre076"	,
                   '100'	=	"Gre077"	,
                   '101'	=	"Gre078"	,
                   '102'	=	"Gre079"	,
                   '103'	=	"Gre080"	,
                   '104'	=	"Gre081"	,
                   '105'	=	"Gre082"	,
                   '106'	=	"Gre083"	,
                   '107'	=	"Gre025"	,
                   '108'	=	"Gre026"	,
                   '109'	=	"Gre027"	,
                   '110'	=	"Gre028"	,
                   '111'	=	"Gre029"	,
                   '112'	=	"Gre030"	,
                   '113'	=	"Gre031"	,
                   '114'	=	"Gre032"	,
                   '115'	=	"Gre033"	,
                   '116'	=	"Gre034"	,
                   '117'	=	"Gre035"	,
                   '118'	=	"Gre084"	,
                   '119'	=	"Gre085"	,
                   '120'	=	"Gre086"	,
                   '121'	=	"Gre087"	,
                   '122'	=	"Gre088"	,
                   '123'	=	"Hen001"	,
                   '124'	=	"Hen003"	,
                   '125'	=	"Hen007"	,
                   '126'	=	"Hin001"	,
                   '127'	=	"Hin003"	,
                   '128'	=	"Hin005"	,
                   '129'	=	"Hes001"	,
                   '130'	=	"Hes005"	,
                   '131'	=	"Hes006"	,
                   '132'	=	"Hes007"	,
                   '133'	=	"Hes008"	,
                   '134'	=	"Hes009"	,
                   '135'	=	"Hes010"	,
                   '136'	=	"Hes011"	,
                   '137'	=	"Hes012"	,
                   '138'	=	"Hes013"	,
                   '139'	=	"Hes031"	,
                   '140'	=	"Hes032"	,
                   '141'	=	"Hes033"	,
                   '142'	=	"Hes034"	,
                   '143'	=	"Hes035"	,
                   '144'	=	"Hes036"	,
                   '145'	=	"Hit001"	,
                   '146'	=	"Hit002"	,
                   '147'	=	"Hit003"	,
                   '148'	=	"Hit004"	,
                   '149'	=	"Hit005"	,
                   '150'	=	"Hit006"	,
                   '151'	=	"Hit011"	,
                   '152'	=	"Hit013"	,
                   '153'	=	"Hit019"	,
                   '154'	=	"Hol001"	,
                   '155'	=	"Hol002"	,
                   '156'	=	"Hol004"	,
                   '157'	=	"Hol005"	,
                   '158'	=	"Hus001"	,
                   '159'	=	"Hus002"	,
                   '160'	=	"Hus003"	,
                   '161'	=	"Hvi001"	,
                   '162'	=	"Hvi002"	,
                   '163'	=	"Hvi004"	,
                   '164'	=	"Kal005"	,
                   '165'	=	"Kal006"	,
                   '166'	=	"Kal007"	,
                   '167'	=	"Kal008"	,
                   '168'	=	"Kal009"	,
                   '169'	=	"Kel001"	,
                   '170'	=	"Kel002"	,
                   '171'	=	"Kel003"	,
                   '172'	=	"Lei014"	,
                   '173'	=	"Lei021"	,
                   '174'	=	"Lei032"	,
                   '175'	=	"Lei045"	,
                   '176'	=	"Lei052"	,
                   '177'	=	"Lei063"	,
                   '178'	=	"Lei065"	,
                   '179'	=	"Lit001"	,
                   '180'	=	"Lis004"	,
                   '181'	=	"Lis013"	,
                   '182'	=	"Lis014"	,
                   '183'	=	"Lis015"	,
                   '184'	=	"Lis016"	,
                   '185'	=	"Lis017"	,
                   '186'	=	"Lis018"	,
                   '187'	=	"Lis020"	,
                   '188'	=	"Lis021"	,
                   '189'	=	"Lis022"	,
                   '190'	=	"Lis023"	,
                   '191'	=	"Lis024"	,
                   '192'	=	"Lis025"	,
                   '193'	=	"Lis026"	,
                   '194'	=	"Lis027"	,
                   '195'	=	"Lis028"	,
                   '196'	=	"Lis029"	,
                   '197'	=	"Lis030"	,
                   '198'	=	"Lis031"	,
                   '199'	=	"Lis032"	,
                   '200'	=	"Lis033"	,
                   '201'	=	"Lis034"	,
                   '202'	=	"Lis035"	,
                   '203'	=	"Lis036"	,
                   '204'	=	"Lis037"	,
                   '205'	=	"Lis038"	,
                   '206'	=	"Lis039"	,
                   '207'	=	"Lis041"	,
                   '208'	=	"Lis042"	,
                   '209'	=	"Lis043"	,
                   '210'	=	"Lis045"	,
                   '211'	=	"Lis046"	,
                   '212'	=	"Lis048"	,
                   '213'	=	"Lis049"	,
                   '214'	=	"Lis050"	,
                   '215'	=	"Lis051"	,
                   '216'	=	"Lis053"	,
                   '217'	=	"Lis054"	,
                   '218'	=	"Lis055"	,
                   '219'	=	"Lis056"	,
                   '220'	=	"Lis057"	,
                   '221'	=	"Lis058"	,
                   '222'	=	"Lis059"	,
                   '223'	=	"Lis060"	,
                   '224'	=	"Lis061"	,
                   '225'	=	"Lis062"	,
                   '226'	=	"Lis064"	,
                   '227'	=	"Lis065"	,
                   '228'	=	"Mid001"	,
                   '229'	=	"Mid003"	,
                   '230'	=	"Mid005"	,
                   '231'	=	"Mid009"	,
                   '232'	=	"Mid010"	,
                   '233'	=	"Mid012"	,
                   '234'	=	"Mid038"	,
                   '235'	=	"Mid039"	,
                   '236'	=	"Mid040"	,
                   '237'	=	"Mid055"	,
                   '238'	=	"Mid026"	,
                   '239'	=	"Mid028"	,
                   '240'	=	"Mid029"	,
                   '241'	=	"Mid031"	,
                   '242'	=	"Mid032"	,
                   '243'	=	"Mid033"	,
                   '244'	=	"Mid034"	,
                   '245'	=	"Mid036"	,
                   '246'	=	"Mid043"	,
                   '247'	=	"Mid044"	,
                   '248'	=	"Mid045"	,
                   '249'	=	"Mid046"	,
                   '250'	=	"Mid047"	,
                   '251'	=	"Mid048"	,
                   '252'	=	"Mid056"	,
                   '253'	=	"Mid060"	,
                   '254'	=	"Mid016"	,
                   '255'	=	"Mid017"	,
                   '256'	=	"Mid018"	,
                   '257'	=	"Mid022"	,
                   '258'	=	"Mid023"	,
                   '259'	=	"Mid049"	,
                   '260'	=	"Mid051"	,
                   '261'	=	"Mid063"	,
                   '262'	=	"Myv002"	,
                   '263'	=	"Myv003"	,
                   '264'	=	"Myv006"	,
                   '265'	=	"Myv007"	,
                   '266'	=	"Myv008"	,
                   '267'	=	"Myv009"	,
                   '268'	=	"Myv011"	,
                   '269'	=	"Myv012"	,
                   '270'	=	"Myv017"	,
                   '271'	=	"Myv020"	,
                   '272'	=	"Myv021"	,
                   '273'	=	"Nor002"	,
                   '274'	=	"Nor003"	,
                   '275'	=	"Nor007"	,
                   '276'	=	"Nor008"	,
                   '277'	=	"Nor009"	,
                   '278'	=	"Nor010"	,
                   '279'	=	"Nup002"	,
                   '280'	=	"Nup004"	,
                   '281'	=	"Nup007"	,
                   '282'	=	"Nup012"	,
                   '283'	=	"Nup013"	,
                   '284'	=	"Nup014"	,
                   '285'	=	"Nup016"	,
                   '286'	=	"Nup017"	,
                   '287'	=	"Olf015"	,
                   '288'	=	"Olf016"	,
                   '289'	=	"Olf017"	,
                   '290'	=	"Olf001"	,
                   '291'	=	"Olf002"	,
                   '292'	=	"Olf003"	,
                   '293'	=	"Olf010"	,
                   '294'	=	"Olf018"	,
                   '295'	=	"Olf004"	,
                   '296'	=	"Olf005"	,
                   '297'	=	"Olf006"	,
                   '298'	=	"Olf007"	,
                   '299'	=	"Olf008"	,
                   '300'	=	"Olf011"	,
                   '301'	=	"Olf012"	,
                   '302'	=	"Olf013"	,
                   '303'	=	"Fus004"	,
                   '304'	=	"Fus006"	,
                   '305'	=	"Fus008"	,
                   '306'	=	"Fus009"	,
                   '307'	=	"Fus014"	,
                   '308'	=	"Fus019"	,
                   '309'	=	"Fus020"	,
                   '310'	=	"Fus021"	,
                   '311'	=	"Fus024"	,
                   '312'	=	"Fus040"	,
                   '313'	=	"Fus042"	,
                   '314'	=	"Fus027"	,
                   '315'	=	"Fus028"	,
                   '316'	=	"Fus030"	,
                   '317'	=	"Fus034"	,
                   '318'	=	"Fus035"	,
                   '319'	=	"Fus036"	,
                   '320'	=	"Fus038"	,
                   '321'	=	"Fus043"	,
                   '322'	=	"Oxa001"	,
                   '323'	=	"Oxa002"	,
                   '324'	=	"Oxa003"	,
                   '325'	=	"Oxa004"	,
                   '326'	=	"Oxa005"	,
                   '327'	=	"Oxa006"	,
                   '328'	=	"Oxa007"	,
                   '329'	=	"Oxa008"	,
                   '330'	=	"Oxa009"	,
                   '331'	=	"Oxa010"	,
                   '332'	=	"Oxa011"	,
                   '333'	=	"Oxa012"	,
                   '334'	=	"Oxa013"	,
                   '335'	=	"Oxa014"	,
                   '336'	=	"Oxa025"	,
                   '337'	=	"Oxa026"	,
                   '338'	=	"Oxa027"	,
                   '339'	=	"Oxa028"	,
                   '340'	=	"Oxa031"	,
                   '341'	=	"Oxa015"	,
                   '342'	=	"Oxa016"	,
                   '343'	=	"Oxa017"	,
                   '344'	=	"Oxa018"	,
                   '345'	=	"Oxa019"	,
                   '346'	=	"Oxa020"	,
                   '347'	=	"Oxa021"	,
                   '348'	=	"Oxa022"	,
                   '349'	=	"Oxa023"	,
                   '350'	=	"Oxa024"	,
                   '351'	=	"Oxa029"	,
                   '352'	=	"Oxa030"	,
                   '353'	=	"San001"	,
                   '354'	=	"San002"	,
                   '355'	=	"San003"	,
                   '356'	=	"San004"	,
                   '357'	=	"Ska001"	,
                   '358'	=	"Ska002"	,
                   '359'	=	"Ska003"	,
                   '360'	=	"Ska004"	,
                   '361'	=	"Ska005"	,
                   '362'	=	"Ska006"	,
                   '363'	=	"Ska007"	,
                   '364'	=	"Ska008"	,
                   '365'	=	"Slp001"	,
                   '366'	=	"Slp002"	,
                   '367'	=	"Slp004"	,
                   '368'	=	"Slp006"	,
                   '369'	=	"Slp007"	,
                   '370'	=	"Slp008"	,
                   '371'	=	"Slp009"	,
                   '372'	=	"Slp010"	,
                   '373'	=	"Slp011"	,
                   '374'	=	"Slp012"	,
                   '375'	=	"Slp013"	,
                   '376'	=	"Slp014"	,
                   '377'	=	"Slp015"	,
                   '378'	=	"Slp016"	,
                   '379'	=	"Slp017"	,
                   '380'	=	"Slp018"	,
                   '381'	=	"Slp019"	,
                   '382'	=	"Slp020"	,
                   '383'	=	"Slp021"	,
                   '384'	=	"Slp022"	,
                   '385'	=	"Slp023"	,
                   '386'	=	"Slp024"	,
                   '387'	=	"Slp025"	,
                   '388'	=	"Slp026"	,
                   '389'	=	"Slp028"	,
                   '390'	=	"Slp029"	,
                   '391'	=	"Slp030"	,
                   '392'	=	"Slp031"	,
                   '393'	=	"Slp036"	,
                   '394'	=	"Slp038"	,
                   '395'	=	"Slp039"	,
                   '396'	=	"Slp043"	,
                   '397'	=	"Slp046"	,
                   '398'	=	"Slp047"	,
                   '399'	=	"Slp048"	,
                   '400'	=	"Slp052"	,
                   '401'	=	"Slp053"	,
                   '402'	=	"Slp054"	,
                   '403'	=	"Slp057"	,
                   '404'	=	"Slp058"	,
                   '405'	=	"Slp060"	,
                   '406'	=	"Slp064"	,
                   '407'	=	"Sog002"	,
                   '408'	=	"Sog003"	,
                   '409'	=	"Sog004"	,
                   '410'	=	"Sog005"	,
                   '411'	=	"Sog006"	,
                   '412'	=	"Sog008"	,
                   '413'	=	"Sog009"	,
                   '414'	=	"Sog010"	,
                   '415'	=	"Sog011"	,
                   '416'	=	"Sog013"	,
                   '417'	=	"Sog015"	,
                   '418'	=	"Sog016"	,
                   '419'	=	"Sog017"	,
                   '420'	=	"Sog018"	,
                   '421'	=	"Sog019"	,
                   '422'	=	"Sog020"	,
                   '423'	=	"Sog021"	,
                   '424'	=	"Sog022"	,
                   '425'	=	"Sog023"	,
                   '426'	=	"Sta001"	,
                   '427'	=	"Sta002"	,
                   '428'	=	"Sta003"	,
                   '429'	=	"Sta004"	,
                   '430'	=	"Sta005"	,
                   '431'	=	"Sta006"	,
                   '432'	=	"Sta007"	,
                   '433'	=	"Sta008"	,
                   '434'	=	"Sta009"	,
                   '435'	=	"Sta020"	,
                   '436'	=	"Sta021"	,
                   '437'	=	"Sta022"	,
                   '438'	=	"Sta023"	,
                   '439'	=	"Sta024"	,
                   '440'	=	"Sta025"	,
                   '441'	=	"Sta026"	,
                   '442'	=	"Sta029"	,
                   '443'	=	"Sta030"	,
                   '444'	=	"Stf013"	,
                   '445'	=	"Stf014"	,
                   '446'	=	"Stf015"	,
                   '447'	=	"Stf016"	,
                   '448'	=	"Stf017"	,
                   '449'	=	"Stf018"	,
                   '450'	=	"Stf019"	,
                   '451'	=	"Stf020"	,
                   '452'	=	"Stf022"	,
                   '453'	=	"Stf023"	,
                   '454'	=	"Stf024"	,
                   '455'	=	"Stf025"	,
                   '456'	=	"Stf026"	,
                   '457'	=	"Stf027"	,
                   '458'	=	"Stf035"	,
                   '459'	=	"Stf039"	,
                   '460'	=	"Stf040"	,
                   '461'	=	"Sth001"	,
                   '462'	=	"Sth002"	,
                   '463'	=	"Sth005"	,
                   '464'	=	"Sth006"	,
                   '465'	=	"Sth008"	,
                   '466'	=	"Sth009"	,
                   '467'	=	"Sth014"	,
                   '468'	=	"Sth015"	,
                   '469'	=	"Sth016"	,
                   '470'	=	"Sth017"	,
                   '471'	=	"Sth018"	,
                   '472'	=	"Sth019"	,
                   '473'	=	"Sth022"	,
                   '474'	=	"Sth023"	,
                   '475'	=	"Sth024"	,
                   '476'	=	"Sth026"	,
                   '477'	=	"Sth031"	,
                   '478'	=	"Sth032"	,
                   '479'	=	"Sth033"	,
                   '480'	=	"Sth034"	,
                   '481'	=	"Sth035"	,
                   '482'	=	"Sth036"	,
                   '483'	=	"Sun002"	,
                   '484'	=	"Sun005"	,
                   '485'	=	"Sun006"	,
                   '486'	=	"Sun011"	,
                   '487'	=	"Sun014"	,
                   '488'	=	"Sun015"	,
                   '489'	=	"Sun017"	,
                   '490'	=	"Sun018"	,
                   '491'	=	"Svi004"	,
                   '492'	=	"Svi006"	,
                   '493'	=	"Svi017"	,
                   '494'	=	"Svi018"	,
                   '495'	=	"Svi019"	,
                   '496'	=	"Svi020"	,
                   '497'	=	"Svi023"	,
                   '498'	=	"Svi024"	,
                   '499'	=	"Svi036"	,
                   '500'	=	"Thi001"	,
                   '501'	=	"Thi002"	,
                   '502'	=	"Thi003"	,
                   '503'	=	"Thi004"	,
                   '504'	=	"Thi006"	,
                   '505'	=	"Thi007"	,
                   '506'	=	"Thi008"	,
                   '507'	=	"Thi009"	,
                   '508'	=	"Thi010"	,
                   '509'	=	"Thi011"	,
                   '510'	=	"Thi012"	,
                   '511'	=	"Thi013"	,
                   '512'	=	"Thi014"	,
                   '513'	=	"Thi016"	,
                   '514'	=	"Thi020"	,
                   '515'	=	"Thi024"	,
                   '516'	=	"Thi025"	,
                   '517'	=	"Thi059"	,
                   '518'	=	"Thi060"	,
                   '519'	=	"Thi054"	,
                   '520'	=	"Thi055"	,
                   '521'	=	"Thi057"	,
                   '522'	=	"Tho001"	,
                   '523'	=	"Tho002"	,
                   '524'	=	"Tho003"	,
                   '525'	=	"Tho004"	,
                   '526'	=	"Tho005"	,
                   '527'	=	"Tho006"	,
                   '528'	=	"Tho007"	,
                   '529'	=	"Tho008"	,
                   '530'	=	"Tho009"	,
                   '531'	=	"Tho010"	,
                   '532'	=	"Tho011"	,
                   '533'	=	"Tho012"	,
                   '534'	=	"Tho013"	,
                   '535'	=	"Tho015"	,
                   '536'	=	"Tho017"	,
                   '537'	=	"Thl001"	,
                   '538'	=	"Thl002"	,
                   '539'	=	"Thl003"	,
                   '540'	=	"Thl004"	,
                   '541'	=	"Thl005"	,
                   '542'	=	"Thl006"	,
                   '543'	=	"Thl007"	,
                   '544'	=	"Thl008"	,
                   '545'	=	"Thl009"	,
                   '546'	=	"Thl010"	,
                   '547'	=	"Thl012"	,
                   '548'	=	"Thl013"	,
                   '549'	=	"Thl015"	,
                   '550'	=	"Thl016"	,
                   '551'	=	"Thl018"	,
                   '552'	=	"Thv003"	,
                   '553'	=	"Thv007"	,
                   '554'	=	"Thv009"	,
                   '555'	=	"Thv015"	,
                   '556'	=	"Thv018"	,
                   '557'	=	"Thv019"	,
                   '558'	=	"Thv021"	,
                   '559'	=	"Thv022"	,
                   '560'	=	"Thv023"	,
                   '561'	=	"Thv024"	,
                   '562'	=	"Thv026"	,
                   '563'	=	"Thv027"	,
                   '564'	=	"Thv029"	,
                   '565'	=	"Thv030"	,
                   '566'	=	"Ulf001"	,
                   '567'	=	"Ulf003"	,
                   '568'	=	"Ulf004"	,
                   '569'	=	"Ulf005"	,
                   '570'	=	"Ulf006"	,
                   '571'	=	"Ulf010"	,
                   '572'	=	"Ulf014"	,
                   '573'	=	"Ulf015"	,
                   '574'	=	"Ulf016"	,
                   '575'	=	"Ulf022"	,
                   '576'	=	"Ulf024"	,
                   '577'	=	"Ulf030"	,
                   '578'	=	"Ulf032"	,
                   '579'	=	"Var001"	,
                   '580'	=	"Var002"	,
                   '581'	=	"Var003"	,
                   '582'	=	"Var004"	,
                   '583'	=	"Var014"	,
                   '584'	=	"Var005"	,
                   '585'	=	"Var006"	,
                   '586'	=	"Var008"	,
                   '587'	=	"Var009"	,
                   '588'	=	"Var010"	,
                   '589'	=	"Var012"	,
                   '590'	=	"Var016"	,
                   '591'	=	"Var018"	,
                   '592'	=	"Vil002"	,
                   '593'	=	"Vil003"	,
                   '594'	=	"Vil004"	,
                   '595'	=	"Vil005"	,
                   '596'	=	"Vil006"	,
                   '597'	=	"Vil008"	,
                   '598'	=	"Vil010"	,
                   '599'	=	"Vil011"	)
library(dplyr)

### Exporting matrices

best = which.min(cross.entropy(project, K = 15))
matrix15 <- Q(project, K=15, run=best)
write.csv(matrix15, "matrix15")

best = which.min(cross.entropy(project, K = 16))
matrix16 <- Q(project, K=16, run=best)
write.csv(matrix16, "matrix16")

best = which.min(cross.entropy(project, K = 17))
matrix17 <- Q(project, K=17, run=best)
write.csv(matrix17, "matrix17")

best = which.min(cross.entropy(project, K = 18))
matrix18 <- Q(project, K=18, run=best)
write.csv(matrix18, "matrix18")

best = which.min(cross.entropy(project, K = 19))
matrix19 <- Q(project, K=19, run=best)
write.csv(matrix19, "matrix19")

best = which.min(cross.entropy(project, K = 20))
matrix20 <- Q(project, K=20, run=best)
write.csv(matrix20, "matrix20")

best = which.min(cross.entropy(project, K = 21))
matrix21 <- Q(project, K=21, run=best)
write.csv(matrix21, "matrix21")

best = which.min(cross.entropy(project, K = 22))
matrix22 <- Q(project, K=22, run=best)
write.csv(matrix22, "matrix22")

best = which.min(cross.entropy(project, K = 23))
matrix23 <- Q(project, K=23, run=best)
write.csv(matrix23, "matrix23")

project<-load.snmfProject("iceland.snmfProject")


################################ QGIS ################################
project<-load.snmfProject("iceland.snmfProject")

# I think I will use K = 17

best = which.min(cross.entropy(project, K = 17))

qmatrix = Q(project, K=17, run=best)
popmap <- read.table("iceland.QGIS.tsv", header=TRUE, sep="\t")
k<- data.frame(qmatrix)
k$ID <- popmap$POP

library(dplyr)

k %>%
        group_by(k$ID) %>%
        summarize(mean_V1 = mean(V1, na.rm=TRUE),
                  mean_V2 = mean(V2, na.rm=TRUE),
                  mean_V3 = mean(V3, na.rm=TRUE),
                  mean_V4 = mean(V4, na.rm=TRUE),
                  mean_V5 = mean(V5, na.rm=TRUE),
                  mean_V6 = mean(V6, na.rm=TRUE),
                  mean_V7 = mean(V7, na.rm=TRUE),
                  mean_V8 = mean(V8, na.rm=TRUE),
                  mean_V9 = mean(V9, na.rm=TRUE),
                  mean_V10 = mean(V10, na.rm=TRUE),
                  mean_V11 = mean(V11, na.rm=TRUE),
                  mean_V12 = mean(V12, na.rm=TRUE),
                  mean_V13 = mean(V13, na.rm=TRUE),
                  mean_V14 = mean(V14, na.rm=TRUE),
                  mean_V15 = mean(V15, na.rm=TRUE),
                  mean_V16 = mean(V16, na.rm=TRUE),
                  mean_V17 = mean(V17, na.rm=TRUE))

k_summary <- k %>%
        group_by(k$ID) %>%
        summarize(mean_V1 = mean(V1, na.rm=TRUE),
                  mean_V2 = mean(V2, na.rm=TRUE),
                  mean_V3 = mean(V3, na.rm=TRUE),
                  mean_V4 = mean(V4, na.rm=TRUE),
                  mean_V5 = mean(V5, na.rm=TRUE),
                  mean_V6 = mean(V6, na.rm=TRUE),
                  mean_V7 = mean(V7, na.rm=TRUE),
                  mean_V8 = mean(V8, na.rm=TRUE),
                  mean_V9 = mean(V9, na.rm=TRUE),
                  mean_V10 = mean(V10, na.rm=TRUE),
                  mean_V11 = mean(V11, na.rm=TRUE),
                  mean_V12 = mean(V12, na.rm=TRUE),
                  mean_V13 = mean(V13, na.rm=TRUE),
                  mean_V14 = mean(V14, na.rm=TRUE),
                  mean_V15 = mean(V15, na.rm=TRUE),
                  mean_V16 = mean(V16, na.rm=TRUE),
                  mean_V17 = mean(V17, na.rm=TRUE))


coord <- data.frame(popmap$POP, popmap$LAT, popmap$LONG)

coord_summary <- coord %>%
        group_by(popmap.POP) %>%
        summarize(mean_LAT = mean(popmap.LAT, na.rm=TRUE),
                  mean_LONG = mean(popmap.LONG, na.rm=TRUE))

k_number <- k %>%
        group_by(k$ID) %>%
        tally()
k_summary$number <- k_number$n

k_summary

data <- data.frame(k_summary)

data$LAT <- coord_summary$mean_LAT
data$LONG <- coord_summary$mean_LONG

data

##QGIS
par(mfrow=c(1,1))
# spatial
library(raster)
library(rgeos)
library(rgdal)
library(ggmap)
library(sp)

# colors
library(colorspace)

# Read in the state shapefile
IS1<-readOGR(dsn="./Maps/", layer="IS1")
IS2<-readOGR(dsn="./Maps/", layer="IS2")
IS3<-readOGR(dsn="./Maps/", layer="IS3")

plot(IS1)

#If you want to add the other layers then:
plot(IS2, add=TRUE, pch=1, cex=0.01, alpha=0, col="light gray",border="dark gray", yaxs = 'i', xaxs = 'i')
#plot(IS3, add=TRUE, col="light gray", yaxs = 'i', xaxs = 'i')

y<- data$LAT
x<- data$LONG

data$lon=x
data$lat=y
coordinates(data) <- c("lon", "lat")
proj4string(data) <- CRS("+init=epsg:4326") # WGS 84
CRS.new <- CRS("+proj=lcc +lat_1=64.25 +lat_2=65.75 +lat_0=65 +lon_0=-19 +x_0=1700000 +y_0=300000 +ellps=GRS80 +units=m +no_defs")
d <- spTransform(data, CRS.new)
d1 <- data.frame(d)

plot(d, add=TRUE, pch=1, cex=0.01, alpha=0, col="light gray",border="dark gray", yaxs = 'i', xaxs = 'i')

library(maps)
library(plotrix)
library(scales)
library(seqinr)

points(d1$lon, d1$lat, cex = d1$number/1800, col='white', pch=19)

my.colors <- c("#d800d8", "#009700", "#009a9a", "#cc0000", "#6533cb", "#580058", "#125699", "#666666", 
               "#004c4c", "#961919", "#a500a5", "#ff69b4", "#ff814c", "#ff4c00", "#00cccc", "#c1adea",
               "#b3b300")
for (i in 1:(2)) {my.colors[i]<-col2alpha(color=my.colors[i], alpha=1)}


#for proportional sizes
#for (x in 1:nrow(d)){floating.pie(d1$lon[x],d1$lat[x], c(d1$mean_V1[x],d1$mean_V2[x],d1$mean_V3[x],d1$mean_V4[x],d1$mean_V5[x],d1$mean_V6[x],d1$mean_V7[x], d1$mean_V8[x]),radius=d1$number[x]*36, col =my.colors)}

#for same sizes
for (x in 1:nrow(d)){floating.pie(d1$lon[x],d1$lat[x], c(d1$mean_V1[x],d1$mean_V2[x],d1$mean_V3[x],d1$mean_V4[x],d1$mean_V5[x], d1$mean_V6[x],d1$mean_V7[x],d1$mean_V8[x],
                                                         d1$mean_V9[x],d1$mean_V10[x],d1$mean_V11[x],d1$mean_V12[x],d1$mean_V13[x], d1$mean_V14[x],d1$mean_V15[x],d1$mean_V16[x],
                                                         d1$mean_V17[x]),
                                  radius=11000, col =my.colors)}

########################### THINGVALLAVATN INSET

bbox(IS1)

b <- bbox(IS1)
b

b[1,1] <- 1569500
b[1,2] <- 1650000
b[2,1] <- 170000
b[2,2] <- 230000
b <- bbox(t(b))
b

gClip <- function(shp, bb){
        if(class(bb) == "matrix") b_poly <- as(extent(as.vector(t(bb))), "SpatialPolygons")
        else b_poly <- as(extent(bb), "SpatialPolygons")
        gIntersection(shp, b_poly, byid = T)
}

zones_clipped_1 <- gClip(IS1, b)
zones_clipped_2 <- gClip(IS2, b)
zones_clipped_3 <- gClip(IS3, b)

par(mai=c(0.1,0.1,0.1,0.1))
plot(zones_clipped_1, yaxs = 'i', xaxs = 'i', lwd=2)
plot(zones_clipped_2, add=T, col="light gray",border="dark gray", yaxs = 'i', xaxs = 'i')
plot(zones_clipped_3, add=T, col="light gray", yaxs = 'i', xaxs = 'i')


plot(d, add=TRUE, pch=1, cex=0.01, alpha=0, col="light gray",border="dark gray", yaxs = 'i', xaxs = 'i')

library(maps)
library(plotrix)
library(scales)
library(seqinr)

points(d1$lon, d1$lat, cex = d1$number/1800, col='white', pch=19)

my.colors <- c("#d800d8", "#009700", "#009a9a", "#cc0000", "#6533cb", "#580058", "#125699", "#666666", 
               "#004c4c", "#961919", "#a500a5", "#ff69b4", "#ff814c", "#ff4c00", "#00cccc", "#c1adea",
               "#b3b300")
for (i in 1:(2)) {my.colors[i]<-col2alpha(color=my.colors[i], alpha=1)}


#for proportional sizes
#for (x in 1:nrow(d)){floating.pie(d1$lon[x],d1$lat[x], c(d1$mean_V1[x],d1$mean_V2[x],d1$mean_V3[x],d1$mean_V4[x],d1$mean_V5[x],d1$mean_V6[x],d1$mean_V7[x], d1$mean_V8[x]),radius=d1$number[x]*36, col =my.colors)}

#for same sizes
for (x in 1:nrow(d)){floating.pie(d1$lon[x],d1$lat[x], c(d1$mean_V1[x],d1$mean_V2[x],d1$mean_V3[x],d1$mean_V4[x],d1$mean_V5[x], d1$mean_V6[x],d1$mean_V7[x],d1$mean_V8[x],
                                                         d1$mean_V9[x],d1$mean_V10[x],d1$mean_V11[x],d1$mean_V12[x],d1$mean_V13[x], d1$mean_V14[x],d1$mean_V15[x],d1$mean_V16[x],
                                                         d1$mean_V17[x]),
                                  radius=1800, col =my.colors)}


########################### VEIDIVOTN INSET

bbox(IS1)

b <- bbox(IS1)
b

b[1,1] <- 1695000
b[1,2] <- 1723000
b[2,1] <- 194000
b[2,2] <- 226000
b <- bbox(t(b))
b

gClip <- function(shp, bb){
        if(class(bb) == "matrix") b_poly <- as(extent(as.vector(t(bb))), "SpatialPolygons")
        else b_poly <- as(extent(bb), "SpatialPolygons")
        gIntersection(shp, b_poly, byid = T)
}

#zones_clipped_1 <- gClip(IS1, b)
zones_clipped_2 <- gClip(IS2, b)
zones_clipped_3 <- gClip(IS3, b)

par(mai=c(0.1,0.1,0.1,0.1))
#plot(zones_clipped_1, yaxs = 'i', xaxs = 'i', lwd=2)
plot(zones_clipped_2, yaxs = 'i', xaxs = 'i', col="light gray",border="dark gray")
plot(zones_clipped_3, add=T, col="light gray", yaxs = 'i', xaxs = 'i')

plot(d, add=TRUE, pch=1, cex=0.01, alpha=0, col="light gray",border="dark gray", yaxs = 'i', xaxs = 'i')

library(maps)
library(plotrix)
library(scales)
library(seqinr)

points(d1$lon, d1$lat, cex = d1$number/1800, col='white', pch=19)

my.colors <- c("#d800d8", "#009700", "#009a9a", "#cc0000", "#6533cb", "#580058", "#125699", "#666666", 
               "#004c4c", "#961919", "#a500a5", "#ff69b4", "#ff814c", "#ff4c00", "#00cccc", "#c1adea",
               "#b3b300")
for (i in 1:(2)) {my.colors[i]<-col2alpha(color=my.colors[i], alpha=1)}


#for proportional sizes
#for (x in 1:nrow(d)){floating.pie(d1$lon[x],d1$lat[x], c(d1$mean_V1[x],d1$mean_V2[x],d1$mean_V3[x],d1$mean_V4[x],d1$mean_V5[x],d1$mean_V6[x],d1$mean_V7[x], d1$mean_V8[x]),radius=d1$number[x]*36, col =my.colors)}

#for same sizes
for (x in 1:nrow(d)){floating.pie(d1$lon[x],d1$lat[x], c(d1$mean_V1[x],d1$mean_V2[x],d1$mean_V3[x],d1$mean_V4[x],d1$mean_V5[x], d1$mean_V6[x],d1$mean_V7[x],d1$mean_V8[x],
                                                         d1$mean_V9[x],d1$mean_V10[x],d1$mean_V11[x],d1$mean_V12[x],d1$mean_V13[x], d1$mean_V14[x],d1$mean_V15[x],d1$mean_V16[x],
                                                         d1$mean_V17[x]),
                                  radius=1000, col =my.colors)}

