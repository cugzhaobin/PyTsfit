*
* Eq file to account to for change in antenna heights

# Bad clock in receiver (offset in East)

# Bad Snow effects at nett

# Added for some bad processing.


# KunLun Mag 8.1 earthquake event
  eq_def KL 35.92 90.53 1000 11 2001 11 14 01 26
  eq_rename KL

# Sumatra-Andaman Mw9.1 earthquake
  eq_def SA 3.316 95.854  4500 30 2004 12 26 00 58 
  eq_rename SA

  
# Wenchuan Mag 7.9 earthquake event
  eq_def WC 31.099 103.279 800 20 2008 5 12 6 28
  eq_rename WC
  eq_log WC 30 50

# Japan earthquake event
  eq_def JP 38.32  142.36 3500 40 2011 3 11  5 36
  eq_prename JP
# eq_log JP 30 0.05 0.05 0.05 0.05 0.05 0.05
  eq_log JP 30 60

# Yaan Mw6.6 earthquake event
  eq_def YA 30.28  102.95 100  12 2013 4 20  0  2
  eq_prename YA

# Min xian Mw6.6 earthquake event
  eq_def MX 34.52 104.23  100 20  2013 7 22  7  45
  eq_prename MX
# Yutian Mw7.0 earthquake event
  eq_def YT 36.10  82.50  100  12 2014 2 12  9 19 
  eq_prename YT
# Pu er M6.1 earthquake event 
  eq_def PR 23.383 100.470 80 8.5 2014 10 07 13 49
  eq_prename PR

# Nepal Mw7.9 earthquake event
  eq_def NP 28.147 84.708 1000 15 2015 4 25  0 00
  eq_prename NP

# Nepal Mw7.3 earthquake event
  eq_def NL 27.819 86.080 800  15 2015 5 12  7 05
  eq_prename NL
# TA jike Mw7.2
  eq_def TJ 38.39  72.91  600  15 2015 12 07 7 50
  eq_prename TJ
# AKetao M6.7 earthquake event
  eq_def KT 39.27  74.04  150  10 2016 11 25 22 24
  eq_prename KT
# eq_log KT 30 0.05 0.05 0.05 0.05 0.05 0.05
# Hutubi M6.2 earthquake event
  eq_def HB 43.83  86.35  100  6  2016 12 8  5 15
  eq_prename HB
# JingHe M6.7 earthquake event
  eq_def JH 44.302  82.832 100  10 2017 08 08 23 27
  eq_prename JH
#Jiuzhaigou Mw6.6 earthquake event
  eq_def JZ 33.20  103.82  50  15 2017 8  8  13 19
  eq_prename JZ

  eq_def MD 34.59 98.34 400 17 2021 05 21 18 04 19
  eq_prename MD

# Added for few data for CORS stations
  rename BJSH_GPS BJSH_XCL 2000  8  2 0 0  2000  8  2 24 0
  rename DLHA_GPS DLHA_XCL 2000  3 12 0 0  2000  3 12 24 0
  rename KMIN_GPS KMIN_XCL 2003 11 12 0 0  2003 11 12 24 0
  rename KUNM_GPS KUNM_XCL 2001  3 24 0 0  2001  3 24 24 0
  rename KUNM_GPS KUNM_XCL 2002 11 10 0 0  2002 11 10 24 0
  rename KUNM_GPS KUNM_XCL 2005 12 22 0 0  2005 12 22 24 0
  rename CHUN_GPS CHUN_XCL 1999  8  8 0 0  1999  8 15 24 0
  rename XIAA_GPS XIAA_XCL 2011  7  6 0 0  2011  7 10 24 0

# Added for some bad processing for CORS 1999 of CMONOC
  rename PETP_GPS PETP_XCL 1999  1 10 0 0  1999  1 10 24 0
  rename QION_GPS QION_XCL 1999  1 28 0 0  1999  1 28 24 0
  rename TAIN_GPS TAIN_XCL 1999  7  6 0 0  1999  7  6 24 0 
  rename TAIN_GPS TAIN_XCL 1999  8  9 0 0  1999  8  9 24 0 
  rename BJSH_GPS BJSH_XCL 2004  8  8 0 0  2004  8  8 24 0 
  rename LHAS_GPS LHAS_XCL 2004  9 20 0 0  2004  9 20 24 0 

# rename for IGS sites affected by great earthquakes, which are not used
  rename CONT_GPS CONT_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename CONZ_GPS CONZ_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename MIZU_GPS MIZU_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename LYTT_GPS LYTT_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename SAPA_GPS SAPA_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename MTKA_GPS MTKA_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename SANT_GPS SANT_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename USUD_GPS USUD_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename TSK2_GPS TSK2_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename TSKB_GPS TSKB_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename MQZG_GPS MQZG_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename STK2_GPS STK2_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename METS_GPS METS_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename DAEJ_GPS DAEJ_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename SUWN_GPS SUWN_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename DUND_GPS DUND_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename OUSD_GPS OUSD_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename CCJM_GPS CCJM_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename MONP_GPS MONP_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename FLIN_GPS FLIN_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename KAZT_GPS KAZT_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename AIRA_GPS AIRA_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename GUUG_GPS GUUG_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename KHAJ_GPS KHAJ_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename MCIL_GPS MCIL_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename PETS_GPS PETS_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename AREQ_GPS AREQ_XCL 1999  1  1 0 0  2100  1  1  0 0
  rename ASPA_GPS ASPA_XCL 1999  1  1 0 0  2100  1  1  0 0
  rename MANZ_GPS MANZ_XCL 1999  1  1 0 0  2100  1  1  0 0
  rename BOER_GPS BOER_XCL 1999  1  1 0 0  2100  1  1  0 0
# rename XIAN_GPS XIAN_XCL 2010  1  1 0 0  2100  1  1  0 0

# rename for IGS sites spanning short periods
  rename MFKU_GPS MFKU_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename PBR2_GPS PBR2_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename PBRI_GPS PBRI_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename SA40_GPS SA40_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename SBOK_GPS SBOK_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename TDOU_GPS TDOU_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename THU2_GPS THU2_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename TROM_GPS TROM_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename ULDI_GPS ULDI_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename SA42_GPS SA42_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename UMTA_GPS UMTA_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename VLCN_GPS VLCN_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename YIBL_GPS YIBL_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename MFKG_GPS MFKG_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename VALP_GPS VALP_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename ACHO_GPS ACHO_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename CN15_GPS CN15_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename CN32_GPS CN32_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename COCO_GPS COCO_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename DEAR_GPS DEAR_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename EIL2_GPS EIL2_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename FRDN_GPS FRDN_XCL 2009  1  1 0 0  2100  1  1  0 0
  rename GSMQ_GPS GSMQ_XCL 2014  11 1 0 0  2015  12 30 0 0

# Added antenna offset for CMONOC CORS sites
# rename DLHA_GPS DLHA_GPS 2010  8 16 0 0  2100  1  1  0 0 0.0 0.0 -0.01745
# rename DXIN_GPS DXIN_GPS 2010  8  4 0 0  2100  1  1  0 0 0.0 0.0 -0.05454
# rename JIXN_GPS JIXN_GPS 2010  7 10 0 0  2100  1  1  0 0 0.0 0.0 -0.01224
# rename YANC_GPS YANC_GPS 2010  6  4 0 0  2100  1  1  0 0 0.0 0.0 -0.04456
# rename ZHNZ_GPS ZHNZ_GPS 2010  6  4 0 0  2100  1  1  0 0 0.0 0.0 -1.41282
# rename LHAS_GPS LHAS_GPS 2011  3 24 0 0  2100  1  1  0 0 0.0 0.0  0.13555

# rename to exclude some data that are anormal
# rename JLCB_GPS JLCB_XCL 2012 11  1 0 0  2013  4  1  0 0


# Offsets coursed by changing antenna or others
  break BJSH 2000  2 23  0  0
# break CHUN 2000  6  4  0  0
  break DLHA 2000  2 16  0  0
  break DLHA 2010  8 16  0  0
  break DXIN 2000  9  6  0  0
  break DXIN 2010  8  4  0  0
  break GUAN 2010  1 18  0  0
  break GUAN 2014  9 24  0  0
  break HLAR 2010  6  2  0  0
  break JIXN 2010  7 10  0  0
  break LHAS 2007  3 12  0  0
  break LHAS 2010  2 28  0  0
  break LHAS 2011  3 25  0  0
  break LUZH 2010 11 13  0  0
  break QION 2010  6  9  0  0
  break TASH 2010 10  4  0  0
  break TASH 2011  9 29  0  0
  break TASH 2012  4 26  0  0
  break URUM 2008  9 16  0  0
  break URUM 2012  5  8  0  0
  break WUHN 2002  1 26  0  0
  break WUSH 2010  7 15  0  0
  break XIAA 2007  9 14  0  0
  break XIAA 2010  7  6  0  0
  break XIAA 2011  3  9  9  9   #2015-01-12
  break XIAA 2012  6 26  0  0
  break XIAA 2014 11 17  0  0   #2015-01-12
  break XNIN 2001  4  6  0  0
  break XNIN 2010  9  3  0  0
  break YANC 2010  6  4  0  0
  break YONG 2001 12 27  0  0
  break ZHNZ 2002  1  1  0  0
  break ZHNZ 2010  6  4  0  0
  break HLAR 2014  2 27  0  0   #2015-06-04

  break BJGB 2011  1 19  0  0
  break GXGL 2011  5 13  0  0
  break QHYS 2010  8 27  0  0
  break SCSM 2011  1 24  0  0
  break SCMB 2011 12 13  0  0
  break SCML 2010  6 29  0  0
# break SDLY 2010  9 28  0  0
# break SDLY 2010 12  9  0  0
  break SDQD 2010 11 28  0  0
  break SDQD 2010 12 28  0  0
  break SDZB 2010 12 28  0  0
  break YNCX 2010  9 10  0  0
  break XJQH 2012  4 16  0  0
# break XZCD 2011  1  1  0  0
# break HECC 2013  1  1  0  0
# break HLMH 2013  1  1  0  0
  break XZYD 2013 11  3  0  0
  break SXKL 2013 11 30  0  0
  break GSMX 2013  7 22  0  0
  break GDZH 2010 11 16  0  0     #2014-01-12
  break GSQS 2014  8 26  0  0     #2014-01-12
  break MMMZ 2013  9 11  0  0     #2014-01-12
  break NMWH 2013 12 10  0  0     #2014-01-12
  break SNYA 2013 12  6  0  0     #2014-01-12
  

  break HYDE 2007 11 30  0  0
  break GUAM 2002  4 27  0  0
  break NOVM 2012  4 13  0  0
  break ULAB 2011  9 24  0  0
  break GSQS 2016  1 14  0  0
  break YNDC 2016  3 25  0  0
  break GSMQ 2014  5 14  0  0
  break GSMQ 2016  1  1  0  0
  break XJBL 2015 12 07  0  0
  break HEYY 2012  6 11  0  0
  break HEYY 2015  8 16  0  0
  break QHME 2016  1 21  0  0
  break SCBZ 2017  5 12  0  0
  break SCGZ 2014  3 12  0  0
  break SCSN 2014 12 03  0  0
  break SXLQ 2016  3 18  0  0
  break SXXX 2016  5 19  0  0
# Added to estimate earthquake offsets
# break BJFS 2011  3 11  5 36
# break CHAN 2011  3 11  5 36
# break SUIY 2011  3 11  5 36
# break JIXN 2011  3 11  5 36
# break TAIN 2011  3 11  5 36
# break HLAR 2011  3 11  5 36
