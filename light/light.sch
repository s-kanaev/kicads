EESchema Schematic File Version 4
EELAYER 30 0
EELAYER END
$Descr A4 11693 8268
encoding utf-8
Sheet 1 1
Title ""
Date ""
Rev ""
Comp ""
Comment1 ""
Comment2 ""
Comment3 ""
Comment4 ""
$EndDescr
$Comp
L Timer:TLC555CPS U2
U 1 1 60173BDB
P 1500 2300
F 0 "U2" H 1800 2650 50  0000 C CNN
F 1 "TLC555CPS" H 1800 1950 50  0000 C CNN
F 2 "Package_DIP:DIP-8_W7.62mm_LongPads" H 1500 2300 50  0001 C CNN
F 3 "http://www.ti.com/lit/ds/symlink/tlc555.pdf" H 1500 2300 50  0001 C CNN
	1    1500 2300
	1    0    0    -1  
$EndComp
$Comp
L Diode:BAT42 D1
U 1 1 60174126
P 3000 2300
F 0 "D1" H 3000 2516 50  0000 C CNN
F 1 "BAT42" H 3000 2425 50  0000 C CNN
F 2 "Diode_THT:D_DO-35_SOD27_P7.62mm_Horizontal" H 3000 2125 50  0001 C CNN
F 3 "http://www.vishay.com/docs/85660/bat42.pdf" H 3000 2300 50  0001 C CNN
	1    3000 2300
	0    -1   -1   0   
$EndComp
$Comp
L Device:R R2
U 1 1 60259C58
P 3000 2650
F 0 "R2" H 3070 2696 50  0000 L CNN
F 1 "1k" H 3070 2605 50  0000 L CNN
F 2 "Resistor_THT:R_Axial_DIN0207_L6.3mm_D2.5mm_P10.16mm_Horizontal" V 2930 2650 50  0001 C CNN
F 3 "~" H 3000 2650 50  0001 C CNN
	1    3000 2650
	1    0    0    -1  
$EndComp
$Comp
L Device:R R3
U 1 1 60259E93
P 3000 3350
F 0 "R3" H 3070 3396 50  0000 L CNN
F 1 "1k" H 3070 3305 50  0000 L CNN
F 2 "Resistor_THT:R_Axial_DIN0207_L6.3mm_D2.5mm_P10.16mm_Horizontal" V 2930 3350 50  0001 C CNN
F 3 "~" H 3000 3350 50  0001 C CNN
	1    3000 3350
	1    0    0    -1  
$EndComp
$Comp
L Device:R_POT RV1
U 1 1 6025A732
P 3000 3000
F 0 "RV1" H 2930 3046 50  0000 R CNN
F 1 "100k" H 2930 2955 50  0000 R CNN
F 2 "Connector_Wire:SolderWirePad_1x03_P7.62mm_Drill2.5mm" H 3000 3000 50  0001 C CNN
F 3 "~" H 3000 3000 50  0001 C CNN
	1    3000 3000
	-1   0    0    -1  
$EndComp
$Comp
L Device:C C5
U 1 1 6025AB39
P 2200 2750
F 0 "C5" H 2315 2796 50  0000 L CNN
F 1 "33n" H 2315 2705 50  0000 L CNN
F 2 "Capacitor_THT:C_Disc_D10.0mm_W2.5mm_P5.00mm" H 2238 2600 50  0001 C CNN
F 3 "~" H 2200 2750 50  0001 C CNN
	1    2200 2750
	1    0    0    -1  
$EndComp
$Comp
L Device:C C4
U 1 1 6025AEF4
P 3650 1050
F 0 "C4" H 3765 1096 50  0000 L CNN
F 1 "0.1u" H 3765 1005 50  0000 L CNN
F 2 "Capacitor_THT:C_Disc_D10.0mm_W2.5mm_P5.00mm" H 3688 900 50  0001 C CNN
F 3 "~" H 3650 1050 50  0001 C CNN
	1    3650 1050
	1    0    0    -1  
$EndComp
$Comp
L Device:CP C3
U 1 1 6025B399
P 3300 1050
F 0 "C3" H 3418 1096 50  0000 L CNN
F 1 "1u" H 3418 1005 50  0000 L CNN
F 2 "Capacitor_THT:CP_Radial_D8.0mm_P2.50mm" H 3338 900 50  0001 C CNN
F 3 "~" H 3300 1050 50  0001 C CNN
	1    3300 1050
	1    0    0    -1  
$EndComp
$Comp
L Device:C C2
U 1 1 6025BCCE
P 2550 1100
F 0 "C2" H 2665 1146 50  0000 L CNN
F 1 "0.1u" H 2665 1055 50  0000 L CNN
F 2 "Capacitor_THT:C_Disc_D10.0mm_W2.5mm_P5.00mm" H 2588 950 50  0001 C CNN
F 3 "~" H 2550 1100 50  0001 C CNN
	1    2550 1100
	1    0    0    -1  
$EndComp
$Comp
L Device:R R1
U 1 1 6025BFF9
P 3350 2100
F 0 "R1" V 3143 2100 50  0000 C CNN
F 1 "100" V 3234 2100 50  0000 C CNN
F 2 "Resistor_THT:R_Axial_DIN0207_L6.3mm_D2.5mm_P10.16mm_Horizontal" V 3280 2100 50  0001 C CNN
F 3 "~" H 3350 2100 50  0001 C CNN
	1    3350 2100
	0    1    1    0   
$EndComp
$Comp
L Regulator_Linear:L7805 U1
U 1 1 6025C499
P 2000 900
F 0 "U1" H 2000 1142 50  0000 C CNN
F 1 "L7805" H 2000 1051 50  0000 C CNN
F 2 "MY_LIB:TO-220-3_Vertical-mod2" H 2025 750 50  0001 L CIN
F 3 "http://www.st.com/content/ccc/resource/technical/document/datasheet/41/4f/b3/b0/12/d4/47/88/CD00000444.pdf/files/CD00000444.pdf/jcr:content/translations/en.CD00000444.pdf" H 2000 850 50  0001 C CNN
	1    2000 900 
	1    0    0    -1  
$EndComp
$Comp
L power:+12V #PWR0101
U 1 1 60264496
P 700 800
F 0 "#PWR0101" H 700 650 50  0001 C CNN
F 1 "+12V" H 715 973 50  0000 C CNN
F 2 "" H 700 800 50  0001 C CNN
F 3 "" H 700 800 50  0001 C CNN
	1    700  800 
	1    0    0    -1  
$EndComp
$Comp
L power:+5V #PWR0102
U 1 1 60264B91
P 2650 800
F 0 "#PWR0102" H 2650 650 50  0001 C CNN
F 1 "+5V" H 2665 973 50  0000 C CNN
F 2 "" H 2650 800 50  0001 C CNN
F 3 "" H 2650 800 50  0001 C CNN
	1    2650 800 
	1    0    0    -1  
$EndComp
Wire Wire Line
	2650 800  2650 900 
Wire Wire Line
	2650 900  2550 900 
Wire Wire Line
	2550 950  2550 900 
Connection ~ 2550 900 
Wire Wire Line
	2550 900  2300 900 
Wire Wire Line
	2550 1250 2000 1250
Wire Wire Line
	2000 1200 2000 1250
Connection ~ 2000 1250
Wire Wire Line
	2000 1250 1450 1250
$Comp
L power:GND #PWR0103
U 1 1 60266307
P 2000 1350
F 0 "#PWR0103" H 2000 1100 50  0001 C CNN
F 1 "GND" H 2005 1177 50  0000 C CNN
F 2 "" H 2000 1350 50  0001 C CNN
F 3 "" H 2000 1350 50  0001 C CNN
	1    2000 1350
	1    0    0    -1  
$EndComp
Wire Wire Line
	2000 1350 2000 1250
$Comp
L power:+5V #PWR0104
U 1 1 6026ABB8
P 3500 800
F 0 "#PWR0104" H 3500 650 50  0001 C CNN
F 1 "+5V" H 3515 973 50  0000 C CNN
F 2 "" H 3500 800 50  0001 C CNN
F 3 "" H 3500 800 50  0001 C CNN
	1    3500 800 
	1    0    0    -1  
$EndComp
$Comp
L power:GND #PWR0105
U 1 1 6026F637
P 3500 1300
F 0 "#PWR0105" H 3500 1050 50  0001 C CNN
F 1 "GND" H 3505 1127 50  0000 C CNN
F 2 "" H 3500 1300 50  0001 C CNN
F 3 "" H 3500 1300 50  0001 C CNN
	1    3500 1300
	1    0    0    -1  
$EndComp
Wire Wire Line
	3500 1300 3500 1200
Wire Wire Line
	3500 1200 3650 1200
Wire Wire Line
	3300 1200 3500 1200
Connection ~ 3500 1200
Wire Wire Line
	3300 900  3500 900 
Wire Wire Line
	3500 800  3500 900 
Connection ~ 3500 900 
Wire Wire Line
	3500 900  3650 900 
$Comp
L Device:C C1
U 1 1 6025B9CB
P 1450 1100
F 0 "C1" H 1565 1146 50  0000 L CNN
F 1 "0.33u" H 1565 1055 50  0000 L CNN
F 2 "Capacitor_THT:C_Disc_D10.0mm_W2.5mm_P5.00mm" H 1488 950 50  0001 C CNN
F 3 "~" H 1450 1100 50  0001 C CNN
	1    1450 1100
	1    0    0    -1  
$EndComp
$Comp
L Switch:SW_SPST SW1
U 1 1 60274A2D
P 1100 900
F 0 "SW1" H 1100 1135 50  0000 C CNN
F 1 "SW_SPST" H 1100 1044 50  0000 C CNN
F 2 "Connector_Wire:SolderWirePad_1x02_P7.62mm_Drill2.5mm" H 1100 900 50  0001 C CNN
F 3 "~" H 1100 900 50  0001 C CNN
	1    1100 900 
	1    0    0    -1  
$EndComp
Wire Wire Line
	1300 900  1450 900 
Wire Wire Line
	1450 950  1450 900 
Connection ~ 1450 900 
Wire Wire Line
	1450 900  1700 900 
Wire Wire Line
	700  800  700  900 
Wire Wire Line
	700  900  900  900 
$Comp
L power:+5V #PWR0106
U 1 1 60288A10
P 1500 1800
F 0 "#PWR0106" H 1500 1650 50  0001 C CNN
F 1 "+5V" H 1515 1973 50  0000 C CNN
F 2 "" H 1500 1800 50  0001 C CNN
F 3 "" H 1500 1800 50  0001 C CNN
	1    1500 1800
	1    0    0    -1  
$EndComp
Wire Wire Line
	1500 1800 1500 1900
$Comp
L power:GND #PWR0107
U 1 1 6028B3BB
P 1500 2750
F 0 "#PWR0107" H 1500 2500 50  0001 C CNN
F 1 "GND" H 1505 2577 50  0000 C CNN
F 2 "" H 1500 2750 50  0001 C CNN
F 3 "" H 1500 2750 50  0001 C CNN
	1    1500 2750
	1    0    0    -1  
$EndComp
Wire Wire Line
	1500 2750 1500 2700
Text Label 2150 2300 0    50   ~ 0
PWM_disch
Wire Wire Line
	2150 2300 2000 2300
Wire Wire Line
	2000 2100 3000 2100
Wire Wire Line
	3000 2100 3000 2150
Wire Wire Line
	3000 3200 3000 3150
Wire Wire Line
	3000 2850 3000 2800
Wire Wire Line
	3000 2500 3000 2450
Wire Wire Line
	2850 3000 2650 3000
Wire Wire Line
	2650 3000 2650 2500
Wire Wire Line
	2650 2500 2200 2500
Wire Wire Line
	2200 2600 2200 2500
Connection ~ 2200 2500
Wire Wire Line
	2200 2500 2000 2500
Text Label 3000 3650 3    50   ~ 0
PWM_disch
Wire Wire Line
	3000 3650 3000 3500
$Comp
L power:GND #PWR0108
U 1 1 602AB1A8
P 2200 2950
F 0 "#PWR0108" H 2200 2700 50  0001 C CNN
F 1 "GND" H 2205 2777 50  0000 C CNN
F 2 "" H 2200 2950 50  0001 C CNN
F 3 "" H 2200 2950 50  0001 C CNN
	1    2200 2950
	1    0    0    -1  
$EndComp
Wire Wire Line
	2200 2950 2200 2900
$Comp
L power:+5V #PWR0109
U 1 1 602ABED1
P 800 2450
F 0 "#PWR0109" H 800 2300 50  0001 C CNN
F 1 "+5V" H 815 2623 50  0000 C CNN
F 2 "" H 800 2450 50  0001 C CNN
F 3 "" H 800 2450 50  0001 C CNN
	1    800  2450
	1    0    0    -1  
$EndComp
Wire Wire Line
	800  2450 800  2500
Wire Wire Line
	800  2500 1000 2500
Text Label 2350 2500 0    50   ~ 0
PWM_tr
Text Label 650  2100 0    50   ~ 0
PWM_tr
Wire Wire Line
	650  2100 1000 2100
NoConn ~ 1000 2300
$Comp
L Connector:Conn_01x02_Female J1
U 1 1 602AE5B7
P 4600 900
F 0 "J1" H 4628 876 50  0000 L CNN
F 1 "Conn_01x02_Female" H 4628 785 50  0000 L CNN
F 2 "Connector_Wire:SolderWirePad_1x02_P7.62mm_Drill2.5mm" H 4600 900 50  0001 C CNN
F 3 "~" H 4600 900 50  0001 C CNN
	1    4600 900 
	1    0    0    -1  
$EndComp
$Comp
L power:+12V #PWR0110
U 1 1 602AF057
P 4300 800
F 0 "#PWR0110" H 4300 650 50  0001 C CNN
F 1 "+12V" H 4315 973 50  0000 C CNN
F 2 "" H 4300 800 50  0001 C CNN
F 3 "" H 4300 800 50  0001 C CNN
	1    4300 800 
	1    0    0    -1  
$EndComp
Wire Wire Line
	4300 800  4300 900 
Wire Wire Line
	4300 900  4400 900 
$Comp
L power:GND #PWR0111
U 1 1 602AFFC5
P 4300 1050
F 0 "#PWR0111" H 4300 800 50  0001 C CNN
F 1 "GND" H 4305 877 50  0000 C CNN
F 2 "" H 4300 1050 50  0001 C CNN
F 3 "" H 4300 1050 50  0001 C CNN
	1    4300 1050
	1    0    0    -1  
$EndComp
Wire Wire Line
	4300 1050 4300 1000
Wire Wire Line
	4300 1000 4400 1000
Wire Wire Line
	3200 2100 3000 2100
Connection ~ 3000 2100
Text Label 3700 2100 0    50   ~ 0
pwm_out
Wire Wire Line
	3700 2100 3500 2100
$Comp
L LDB600:LDB600 U3
U 1 1 602CC6EB
P 9100 1250
F 0 "U3" H 9100 1815 50  0000 C CNN
F 1 "LDB600" H 9100 1724 50  0000 C CNN
F 2 "MY_LIB:MEANWELL_LDB-350L" H 8550 900 50  0001 C CNN
F 3 "" H 8550 900 50  0001 C CNN
	1    9100 1250
	1    0    0    -1  
$EndComp
$Comp
L Connector:Conn_01x02_Male J5
U 1 1 602D70BA
P 10300 1250
F 0 "J5" H 10272 1132 50  0000 R CNN
F 1 "Conn_01x02_Male" H 10272 1223 50  0000 R CNN
F 2 "Connector_Wire:SolderWirePad_1x02_P7.62mm_Drill2.5mm" H 10300 1250 50  0001 C CNN
F 3 "~" H 10300 1250 50  0001 C CNN
	1    10300 1250
	-1   0    0    1   
$EndComp
$Comp
L Connector:Conn_01x02_Female J3
U 1 1 602D80E9
P 7850 1100
F 0 "J3" H 7742 775 50  0000 C CNN
F 1 "Conn_01x02_Female" H 7742 866 50  0000 C CNN
F 2 "Connector_Wire:SolderWirePad_1x02_P7.62mm_Drill2.5mm" H 7850 1100 50  0001 C CNN
F 3 "~" H 7850 1100 50  0001 C CNN
	1    7850 1100
	-1   0    0    1   
$EndComp
Wire Wire Line
	8550 1000 8550 1100
Wire Wire Line
	8450 1100 8450 1250
Wire Wire Line
	8450 1250 8550 1250
Wire Wire Line
	8550 1250 8550 1350
Connection ~ 8550 1250
Wire Wire Line
	8050 1400 8350 1400
Wire Wire Line
	8350 1400 8350 1350
Wire Wire Line
	8350 1350 8550 1350
Wire Wire Line
	8050 1500 8550 1500
Wire Wire Line
	9650 1000 9650 1100
Wire Wire Line
	9650 1250 9650 1350
Wire Wire Line
	9650 1100 10100 1100
Wire Wire Line
	10100 1100 10100 1150
Wire Wire Line
	10100 1250 9650 1250
Connection ~ 9650 1250
Connection ~ 9650 1100
Connection ~ 8550 1350
Wire Wire Line
	8050 1000 8050 800 
Wire Wire Line
	8550 800  8550 1000
Connection ~ 8550 1000
$Comp
L Connector:Conn_01x02_Female J4
U 1 1 602D019D
P 7850 1400
F 0 "J4" H 7742 1585 50  0000 C CNN
F 1 "Conn_01x02_Female" H 7742 1494 50  0000 C CNN
F 2 "Connector_Wire:SolderWirePad_1x02_P7.62mm_Drill2.5mm" H 7850 1400 50  0001 C CNN
F 3 "~" H 7850 1400 50  0001 C CNN
	1    7850 1400
	-1   0    0    -1  
$EndComp
Connection ~ 8300 1100
Connection ~ 8300 800 
Wire Wire Line
	8300 800  8550 800 
Wire Wire Line
	8050 800  8300 800 
Wire Wire Line
	8300 1100 8450 1100
Wire Wire Line
	8050 1100 8300 1100
$Comp
L Device:C C6
U 1 1 6032745E
P 8300 950
F 0 "C6" H 8415 996 50  0000 L CNN
F 1 "1u" H 8415 905 50  0000 L CNN
F 2 "Capacitor_THT:C_Disc_D10.0mm_W2.5mm_P5.00mm" H 8338 800 50  0001 C CNN
F 3 "~" H 8300 950 50  0001 C CNN
	1    8300 950 
	1    0    0    -1  
$EndComp
$Comp
L LDB600:LDB600 U4
U 1 1 60298E38
P 9100 2250
F 0 "U4" H 9100 2815 50  0000 C CNN
F 1 "LDB600" H 9100 2724 50  0000 C CNN
F 2 "MY_LIB:MEANWELL_LDB-350L" H 8550 1900 50  0001 C CNN
F 3 "" H 8550 1900 50  0001 C CNN
	1    9100 2250
	1    0    0    -1  
$EndComp
$Comp
L Connector:Conn_01x02_Male J8
U 1 1 60298E3E
P 10300 2250
F 0 "J8" H 10272 2132 50  0000 R CNN
F 1 "Conn_01x02_Male" H 10272 2223 50  0000 R CNN
F 2 "Connector_Wire:SolderWirePad_1x02_P7.62mm_Drill2.5mm" H 10300 2250 50  0001 C CNN
F 3 "~" H 10300 2250 50  0001 C CNN
	1    10300 2250
	-1   0    0    1   
$EndComp
$Comp
L Connector:Conn_01x02_Female J6
U 1 1 60298E44
P 7850 2100
F 0 "J6" H 7742 1775 50  0000 C CNN
F 1 "Conn_01x02_Female" H 7742 1866 50  0000 C CNN
F 2 "Connector_Wire:SolderWirePad_1x02_P7.62mm_Drill2.5mm" H 7850 2100 50  0001 C CNN
F 3 "~" H 7850 2100 50  0001 C CNN
	1    7850 2100
	-1   0    0    1   
$EndComp
Wire Wire Line
	8550 2000 8550 2100
Wire Wire Line
	8450 2100 8450 2250
Wire Wire Line
	8450 2250 8550 2250
Wire Wire Line
	8550 2250 8550 2350
Connection ~ 8550 2250
Wire Wire Line
	8050 2400 8350 2400
Wire Wire Line
	8350 2400 8350 2350
Wire Wire Line
	8350 2350 8550 2350
Wire Wire Line
	8050 2500 8550 2500
Wire Wire Line
	9650 2000 9650 2100
Wire Wire Line
	9650 2250 9650 2350
Wire Wire Line
	9650 2100 10100 2100
Wire Wire Line
	10100 2100 10100 2150
Wire Wire Line
	10100 2250 9650 2250
Connection ~ 9650 2250
Connection ~ 9650 2100
Connection ~ 8550 2350
Wire Wire Line
	8050 2000 8050 1800
Wire Wire Line
	8550 1800 8550 2000
Connection ~ 8550 2000
$Comp
L Connector:Conn_01x02_Female J7
U 1 1 60298E5E
P 7850 2400
F 0 "J7" H 7742 2585 50  0000 C CNN
F 1 "Conn_01x02_Female" H 7742 2494 50  0000 C CNN
F 2 "Connector_Wire:SolderWirePad_1x02_P7.62mm_Drill2.5mm" H 7850 2400 50  0001 C CNN
F 3 "~" H 7850 2400 50  0001 C CNN
	1    7850 2400
	-1   0    0    -1  
$EndComp
Connection ~ 8300 2100
Connection ~ 8300 1800
Wire Wire Line
	8300 1800 8550 1800
Wire Wire Line
	8050 1800 8300 1800
Wire Wire Line
	8300 2100 8450 2100
Wire Wire Line
	8050 2100 8300 2100
$Comp
L Device:C C7
U 1 1 60298E6A
P 8300 1950
F 0 "C7" H 8415 1996 50  0000 L CNN
F 1 "1u" H 8415 1905 50  0000 L CNN
F 2 "Capacitor_THT:C_Disc_D10.0mm_W2.5mm_P5.00mm" H 8338 1800 50  0001 C CNN
F 3 "~" H 8300 1950 50  0001 C CNN
	1    8300 1950
	1    0    0    -1  
$EndComp
$Comp
L LDB600:LDB600 U5
U 1 1 6029CF02
P 9100 3250
F 0 "U5" H 9100 3815 50  0000 C CNN
F 1 "LDB600" H 9100 3724 50  0000 C CNN
F 2 "MY_LIB:MEANWELL_LDB-350L" H 8550 2900 50  0001 C CNN
F 3 "" H 8550 2900 50  0001 C CNN
	1    9100 3250
	1    0    0    -1  
$EndComp
$Comp
L Connector:Conn_01x02_Male J11
U 1 1 6029CF08
P 10300 3250
F 0 "J11" H 10272 3132 50  0000 R CNN
F 1 "Conn_01x02_Male" H 10272 3223 50  0000 R CNN
F 2 "Connector_Wire:SolderWirePad_1x02_P7.62mm_Drill2.5mm" H 10300 3250 50  0001 C CNN
F 3 "~" H 10300 3250 50  0001 C CNN
	1    10300 3250
	-1   0    0    1   
$EndComp
$Comp
L Connector:Conn_01x02_Female J9
U 1 1 6029CF0E
P 7850 3100
F 0 "J9" H 7742 2775 50  0000 C CNN
F 1 "Conn_01x02_Female" H 7742 2866 50  0000 C CNN
F 2 "Connector_Wire:SolderWirePad_1x02_P7.62mm_Drill2.5mm" H 7850 3100 50  0001 C CNN
F 3 "~" H 7850 3100 50  0001 C CNN
	1    7850 3100
	-1   0    0    1   
$EndComp
Wire Wire Line
	8550 3000 8550 3100
Wire Wire Line
	8450 3100 8450 3250
Wire Wire Line
	8450 3250 8550 3250
Wire Wire Line
	8550 3250 8550 3350
Connection ~ 8550 3250
Wire Wire Line
	8050 3400 8350 3400
Wire Wire Line
	8350 3400 8350 3350
Wire Wire Line
	8350 3350 8550 3350
Wire Wire Line
	8050 3500 8550 3500
Wire Wire Line
	9650 3000 9650 3100
Wire Wire Line
	9650 3250 9650 3350
Wire Wire Line
	9650 3100 10100 3100
Wire Wire Line
	10100 3100 10100 3150
Wire Wire Line
	10100 3250 9650 3250
Connection ~ 9650 3250
Connection ~ 9650 3100
Connection ~ 8550 3350
Wire Wire Line
	8050 3000 8050 2800
Wire Wire Line
	8550 2800 8550 3000
Connection ~ 8550 3000
$Comp
L Connector:Conn_01x02_Female J10
U 1 1 6029CF28
P 7850 3400
F 0 "J10" H 7742 3585 50  0000 C CNN
F 1 "Conn_01x02_Female" H 7742 3494 50  0000 C CNN
F 2 "Connector_Wire:SolderWirePad_1x02_P7.62mm_Drill2.5mm" H 7850 3400 50  0001 C CNN
F 3 "~" H 7850 3400 50  0001 C CNN
	1    7850 3400
	-1   0    0    -1  
$EndComp
Connection ~ 8300 3100
Connection ~ 8300 2800
Wire Wire Line
	8300 2800 8550 2800
Wire Wire Line
	8050 2800 8300 2800
Wire Wire Line
	8300 3100 8450 3100
Wire Wire Line
	8050 3100 8300 3100
$Comp
L Device:C C8
U 1 1 6029CF34
P 8300 2950
F 0 "C8" H 8415 2996 50  0000 L CNN
F 1 "1u" H 8415 2905 50  0000 L CNN
F 2 "Capacitor_THT:C_Disc_D10.0mm_W2.5mm_P5.00mm" H 8338 2800 50  0001 C CNN
F 3 "~" H 8300 2950 50  0001 C CNN
	1    8300 2950
	1    0    0    -1  
$EndComp
$Comp
L LDB600:LDB600 U6
U 1 1 602A0D59
P 9100 4250
F 0 "U6" H 9100 4815 50  0000 C CNN
F 1 "LDB600" H 9100 4724 50  0000 C CNN
F 2 "MY_LIB:MEANWELL_LDB-350L" H 8550 3900 50  0001 C CNN
F 3 "" H 8550 3900 50  0001 C CNN
	1    9100 4250
	1    0    0    -1  
$EndComp
$Comp
L Connector:Conn_01x02_Male J14
U 1 1 602A0D5F
P 10300 4250
F 0 "J14" H 10272 4132 50  0000 R CNN
F 1 "Conn_01x02_Male" H 10272 4223 50  0000 R CNN
F 2 "Connector_Wire:SolderWirePad_1x02_P7.62mm_Drill2.5mm" H 10300 4250 50  0001 C CNN
F 3 "~" H 10300 4250 50  0001 C CNN
	1    10300 4250
	-1   0    0    1   
$EndComp
$Comp
L Connector:Conn_01x02_Female J12
U 1 1 602A0D65
P 7850 4100
F 0 "J12" H 7742 3775 50  0000 C CNN
F 1 "Conn_01x02_Female" H 7742 3866 50  0000 C CNN
F 2 "Connector_Wire:SolderWirePad_1x02_P7.62mm_Drill2.5mm" H 7850 4100 50  0001 C CNN
F 3 "~" H 7850 4100 50  0001 C CNN
	1    7850 4100
	-1   0    0    1   
$EndComp
Wire Wire Line
	8550 4000 8550 4100
Wire Wire Line
	8450 4100 8450 4250
Wire Wire Line
	8450 4250 8550 4250
Wire Wire Line
	8550 4250 8550 4350
Connection ~ 8550 4250
Wire Wire Line
	8050 4400 8350 4400
Wire Wire Line
	8350 4400 8350 4350
Wire Wire Line
	8350 4350 8550 4350
Wire Wire Line
	8050 4500 8550 4500
Wire Wire Line
	9650 4000 9650 4100
Wire Wire Line
	9650 4250 9650 4350
Wire Wire Line
	9650 4100 10100 4100
Wire Wire Line
	10100 4100 10100 4150
Wire Wire Line
	10100 4250 9650 4250
Connection ~ 9650 4250
Connection ~ 9650 4100
Connection ~ 8550 4350
Wire Wire Line
	8050 4000 8050 3800
Wire Wire Line
	8550 3800 8550 4000
Connection ~ 8550 4000
$Comp
L Connector:Conn_01x02_Female J13
U 1 1 602A0D7F
P 7850 4400
F 0 "J13" H 7742 4585 50  0000 C CNN
F 1 "Conn_01x02_Female" H 7742 4494 50  0000 C CNN
F 2 "Connector_Wire:SolderWirePad_1x02_P7.62mm_Drill2.5mm" H 7850 4400 50  0001 C CNN
F 3 "~" H 7850 4400 50  0001 C CNN
	1    7850 4400
	-1   0    0    -1  
$EndComp
Connection ~ 8300 4100
Connection ~ 8300 3800
Wire Wire Line
	8300 3800 8550 3800
Wire Wire Line
	8050 3800 8300 3800
Wire Wire Line
	8300 4100 8450 4100
Wire Wire Line
	8050 4100 8300 4100
$Comp
L Device:C C9
U 1 1 602A0D8B
P 8300 3950
F 0 "C9" H 8415 3996 50  0000 L CNN
F 1 "1u" H 8415 3905 50  0000 L CNN
F 2 "Capacitor_THT:C_Disc_D10.0mm_W2.5mm_P5.00mm" H 8338 3800 50  0001 C CNN
F 3 "~" H 8300 3950 50  0001 C CNN
	1    8300 3950
	1    0    0    -1  
$EndComp
$Comp
L LDB600:LDB600 U7
U 1 1 602A5D2C
P 9100 5250
F 0 "U7" H 9100 5815 50  0000 C CNN
F 1 "LDB600" H 9100 5724 50  0000 C CNN
F 2 "MY_LIB:MEANWELL_LDB-350L" H 8550 4900 50  0001 C CNN
F 3 "" H 8550 4900 50  0001 C CNN
	1    9100 5250
	1    0    0    -1  
$EndComp
$Comp
L Connector:Conn_01x02_Male J17
U 1 1 602A5D32
P 10300 5250
F 0 "J17" H 10272 5132 50  0000 R CNN
F 1 "Conn_01x02_Male" H 10272 5223 50  0000 R CNN
F 2 "Connector_Wire:SolderWirePad_1x02_P7.62mm_Drill2.5mm" H 10300 5250 50  0001 C CNN
F 3 "~" H 10300 5250 50  0001 C CNN
	1    10300 5250
	-1   0    0    1   
$EndComp
$Comp
L Connector:Conn_01x02_Female J15
U 1 1 602A5D38
P 7850 5100
F 0 "J15" H 7742 4775 50  0000 C CNN
F 1 "Conn_01x02_Female" H 7742 4866 50  0000 C CNN
F 2 "Connector_Wire:SolderWirePad_1x02_P7.62mm_Drill2.5mm" H 7850 5100 50  0001 C CNN
F 3 "~" H 7850 5100 50  0001 C CNN
	1    7850 5100
	-1   0    0    1   
$EndComp
Wire Wire Line
	8550 5000 8550 5100
Wire Wire Line
	8450 5100 8450 5250
Wire Wire Line
	8450 5250 8550 5250
Wire Wire Line
	8550 5250 8550 5350
Connection ~ 8550 5250
Wire Wire Line
	8050 5400 8350 5400
Wire Wire Line
	8350 5400 8350 5350
Wire Wire Line
	8350 5350 8550 5350
Wire Wire Line
	8050 5500 8550 5500
Wire Wire Line
	9650 5000 9650 5100
Wire Wire Line
	9650 5250 9650 5350
Wire Wire Line
	9650 5100 10100 5100
Wire Wire Line
	10100 5100 10100 5150
Wire Wire Line
	10100 5250 9650 5250
Connection ~ 9650 5250
Connection ~ 9650 5100
Connection ~ 8550 5350
Wire Wire Line
	8050 5000 8050 4800
Wire Wire Line
	8550 4800 8550 5000
Connection ~ 8550 5000
$Comp
L Connector:Conn_01x02_Female J16
U 1 1 602A5D52
P 7850 5400
F 0 "J16" H 7742 5585 50  0000 C CNN
F 1 "Conn_01x02_Female" H 7742 5494 50  0000 C CNN
F 2 "Connector_Wire:SolderWirePad_1x02_P7.62mm_Drill2.5mm" H 7850 5400 50  0001 C CNN
F 3 "~" H 7850 5400 50  0001 C CNN
	1    7850 5400
	-1   0    0    -1  
$EndComp
Connection ~ 8300 5100
Connection ~ 8300 4800
Wire Wire Line
	8300 4800 8550 4800
Wire Wire Line
	8050 4800 8300 4800
Wire Wire Line
	8300 5100 8450 5100
Wire Wire Line
	8050 5100 8300 5100
$Comp
L Device:C C10
U 1 1 602A5D5E
P 8300 4950
F 0 "C10" H 8415 4996 50  0000 L CNN
F 1 "1u" H 8415 4905 50  0000 L CNN
F 2 "Capacitor_THT:C_Disc_D10.0mm_W2.5mm_P5.00mm" H 8338 4800 50  0001 C CNN
F 3 "~" H 8300 4950 50  0001 C CNN
	1    8300 4950
	1    0    0    -1  
$EndComp
$Comp
L LDB600:LDB600 U8
U 1 1 602C785D
P 2250 4800
F 0 "U8" H 2250 5365 50  0000 C CNN
F 1 "LDB600" H 2250 5274 50  0000 C CNN
F 2 "MY_LIB:MEANWELL_LDB-350L" H 1700 4450 50  0001 C CNN
F 3 "" H 1700 4450 50  0001 C CNN
	1    2250 4800
	1    0    0    -1  
$EndComp
$Comp
L Connector:Conn_01x02_Male J20
U 1 1 602C7863
P 3450 4800
F 0 "J20" H 3422 4682 50  0000 R CNN
F 1 "Conn_01x02_Male" H 3422 4773 50  0000 R CNN
F 2 "Connector_Wire:SolderWirePad_1x02_P7.62mm_Drill2.5mm" H 3450 4800 50  0001 C CNN
F 3 "~" H 3450 4800 50  0001 C CNN
	1    3450 4800
	-1   0    0    1   
$EndComp
$Comp
L Connector:Conn_01x02_Female J18
U 1 1 602C7869
P 1000 4650
F 0 "J18" H 892 4325 50  0000 C CNN
F 1 "Conn_01x02_Female" H 892 4416 50  0000 C CNN
F 2 "Connector_Wire:SolderWirePad_1x02_P7.62mm_Drill2.5mm" H 1000 4650 50  0001 C CNN
F 3 "~" H 1000 4650 50  0001 C CNN
	1    1000 4650
	-1   0    0    1   
$EndComp
Wire Wire Line
	1700 4550 1700 4650
Wire Wire Line
	1600 4650 1600 4800
Wire Wire Line
	1600 4800 1700 4800
Wire Wire Line
	1700 4800 1700 4900
Connection ~ 1700 4800
Wire Wire Line
	1200 4950 1500 4950
Wire Wire Line
	1500 4950 1500 4900
Wire Wire Line
	1500 4900 1700 4900
Wire Wire Line
	1200 5050 1700 5050
Wire Wire Line
	2800 4550 2800 4650
Wire Wire Line
	2800 4800 2800 4900
Wire Wire Line
	2800 4650 3250 4650
Wire Wire Line
	3250 4650 3250 4700
Wire Wire Line
	3250 4800 2800 4800
Connection ~ 2800 4800
Connection ~ 2800 4650
Connection ~ 1700 4900
Wire Wire Line
	1200 4550 1200 4350
Wire Wire Line
	1700 4350 1700 4550
Connection ~ 1700 4550
$Comp
L Connector:Conn_01x02_Female J19
U 1 1 602C7883
P 1000 4950
F 0 "J19" H 892 5135 50  0000 C CNN
F 1 "Conn_01x02_Female" H 892 5044 50  0000 C CNN
F 2 "Connector_Wire:SolderWirePad_1x02_P7.62mm_Drill2.5mm" H 1000 4950 50  0001 C CNN
F 3 "~" H 1000 4950 50  0001 C CNN
	1    1000 4950
	-1   0    0    -1  
$EndComp
Connection ~ 1450 4650
Connection ~ 1450 4350
Wire Wire Line
	1450 4350 1700 4350
Wire Wire Line
	1200 4350 1450 4350
Wire Wire Line
	1450 4650 1600 4650
Wire Wire Line
	1200 4650 1450 4650
$Comp
L Device:C C11
U 1 1 602C788F
P 1450 4500
F 0 "C11" H 1565 4546 50  0000 L CNN
F 1 "1u" H 1565 4455 50  0000 L CNN
F 2 "Capacitor_THT:C_Disc_D10.0mm_W2.5mm_P5.00mm" H 1488 4350 50  0001 C CNN
F 3 "~" H 1450 4500 50  0001 C CNN
	1    1450 4500
	1    0    0    -1  
$EndComp
$Comp
L LDB600:LDB600 U9
U 1 1 602D0BF9
P 2250 5800
F 0 "U9" H 2250 6365 50  0000 C CNN
F 1 "LDB600" H 2250 6274 50  0000 C CNN
F 2 "MY_LIB:MEANWELL_LDB-350L" H 1700 5450 50  0001 C CNN
F 3 "" H 1700 5450 50  0001 C CNN
	1    2250 5800
	1    0    0    -1  
$EndComp
$Comp
L Connector:Conn_01x02_Male J23
U 1 1 602D0BFF
P 3450 5800
F 0 "J23" H 3422 5682 50  0000 R CNN
F 1 "Conn_01x02_Male" H 3422 5773 50  0000 R CNN
F 2 "Connector_Wire:SolderWirePad_1x02_P7.62mm_Drill2.5mm" H 3450 5800 50  0001 C CNN
F 3 "~" H 3450 5800 50  0001 C CNN
	1    3450 5800
	-1   0    0    1   
$EndComp
$Comp
L Connector:Conn_01x02_Female J21
U 1 1 602D0C05
P 1000 5650
F 0 "J21" H 892 5325 50  0000 C CNN
F 1 "Conn_01x02_Female" H 892 5416 50  0000 C CNN
F 2 "Connector_Wire:SolderWirePad_1x02_P7.62mm_Drill2.5mm" H 1000 5650 50  0001 C CNN
F 3 "~" H 1000 5650 50  0001 C CNN
	1    1000 5650
	-1   0    0    1   
$EndComp
Wire Wire Line
	1700 5550 1700 5650
Wire Wire Line
	1600 5650 1600 5800
Wire Wire Line
	1600 5800 1700 5800
Wire Wire Line
	1700 5800 1700 5900
Connection ~ 1700 5800
Wire Wire Line
	1200 5950 1500 5950
Wire Wire Line
	1500 5950 1500 5900
Wire Wire Line
	1500 5900 1700 5900
Wire Wire Line
	1200 6050 1700 6050
Wire Wire Line
	2800 5550 2800 5650
Wire Wire Line
	2800 5800 2800 5900
Wire Wire Line
	2800 5650 3250 5650
Wire Wire Line
	3250 5650 3250 5700
Wire Wire Line
	3250 5800 2800 5800
Connection ~ 2800 5800
Connection ~ 2800 5650
Connection ~ 1700 5900
Wire Wire Line
	1200 5550 1200 5350
Wire Wire Line
	1700 5350 1700 5550
Connection ~ 1700 5550
$Comp
L Connector:Conn_01x02_Female J22
U 1 1 602D0C1F
P 1000 5950
F 0 "J22" H 892 6135 50  0000 C CNN
F 1 "Conn_01x02_Female" H 892 6044 50  0000 C CNN
F 2 "Connector_Wire:SolderWirePad_1x02_P7.62mm_Drill2.5mm" H 1000 5950 50  0001 C CNN
F 3 "~" H 1000 5950 50  0001 C CNN
	1    1000 5950
	-1   0    0    -1  
$EndComp
Connection ~ 1450 5650
Connection ~ 1450 5350
Wire Wire Line
	1450 5350 1700 5350
Wire Wire Line
	1200 5350 1450 5350
Wire Wire Line
	1450 5650 1600 5650
Wire Wire Line
	1200 5650 1450 5650
$Comp
L Device:C C12
U 1 1 602D0C2B
P 1450 5500
F 0 "C12" H 1565 5546 50  0000 L CNN
F 1 "1u" H 1565 5455 50  0000 L CNN
F 2 "Capacitor_THT:C_Disc_D10.0mm_W2.5mm_P5.00mm" H 1488 5350 50  0001 C CNN
F 3 "~" H 1450 5500 50  0001 C CNN
	1    1450 5500
	1    0    0    -1  
$EndComp
$Comp
L LDB600:LDB600 U10
U 1 1 602DB58F
P 2250 6800
F 0 "U10" H 2250 7365 50  0000 C CNN
F 1 "LDB600" H 2250 7274 50  0000 C CNN
F 2 "MY_LIB:MEANWELL_LDB-350L" H 1700 6450 50  0001 C CNN
F 3 "" H 1700 6450 50  0001 C CNN
	1    2250 6800
	1    0    0    -1  
$EndComp
$Comp
L Connector:Conn_01x02_Male J26
U 1 1 602DB595
P 3450 6800
F 0 "J26" H 3422 6682 50  0000 R CNN
F 1 "Conn_01x02_Male" H 3422 6773 50  0000 R CNN
F 2 "Connector_Wire:SolderWirePad_1x02_P7.62mm_Drill2.5mm" H 3450 6800 50  0001 C CNN
F 3 "~" H 3450 6800 50  0001 C CNN
	1    3450 6800
	-1   0    0    1   
$EndComp
$Comp
L Connector:Conn_01x02_Female J24
U 1 1 602DB59B
P 1000 6650
F 0 "J24" H 892 6325 50  0000 C CNN
F 1 "Conn_01x02_Female" H 892 6416 50  0000 C CNN
F 2 "Connector_Wire:SolderWirePad_1x02_P7.62mm_Drill2.5mm" H 1000 6650 50  0001 C CNN
F 3 "~" H 1000 6650 50  0001 C CNN
	1    1000 6650
	-1   0    0    1   
$EndComp
Wire Wire Line
	1700 6550 1700 6650
Wire Wire Line
	1600 6650 1600 6800
Wire Wire Line
	1600 6800 1700 6800
Wire Wire Line
	1700 6800 1700 6900
Connection ~ 1700 6800
Wire Wire Line
	1200 6950 1500 6950
Wire Wire Line
	1500 6950 1500 6900
Wire Wire Line
	1500 6900 1700 6900
Wire Wire Line
	1200 7050 1700 7050
Wire Wire Line
	2800 6550 2800 6650
Wire Wire Line
	2800 6800 2800 6900
Wire Wire Line
	2800 6650 3250 6650
Wire Wire Line
	3250 6650 3250 6700
Wire Wire Line
	3250 6800 2800 6800
Connection ~ 2800 6800
Connection ~ 2800 6650
Connection ~ 1700 6900
Wire Wire Line
	1200 6550 1200 6350
Wire Wire Line
	1700 6350 1700 6550
Connection ~ 1700 6550
Connection ~ 1450 6650
Connection ~ 1450 6350
Wire Wire Line
	1450 6350 1700 6350
Wire Wire Line
	1200 6350 1450 6350
Wire Wire Line
	1450 6650 1600 6650
Wire Wire Line
	1200 6650 1450 6650
$Comp
L Device:C C13
U 1 1 602DB5C1
P 1450 6500
F 0 "C13" H 1565 6546 50  0000 L CNN
F 1 "1u" H 1565 6455 50  0000 L CNN
F 2 "Capacitor_THT:C_Disc_D10.0mm_W2.5mm_P5.00mm" H 1488 6350 50  0001 C CNN
F 3 "~" H 1450 6500 50  0001 C CNN
	1    1450 6500
	1    0    0    -1  
$EndComp
$Comp
L Connector:Conn_01x02_Female J25
U 1 1 602DB5B5
P 1000 6950
F 0 "J25" H 892 7135 50  0000 C CNN
F 1 "Conn_01x02_Female" H 892 7044 50  0000 C CNN
F 2 "Connector_Wire:SolderWirePad_1x02_P7.62mm_Drill2.5mm" H 1000 6950 50  0001 C CNN
F 3 "~" H 1000 6950 50  0001 C CNN
	1    1000 6950
	-1   0    0    -1  
$EndComp
Wire Wire Line
	6100 950  5700 950 
Text Label 5700 950  0    50   ~ 0
pwm_out
$Comp
L power:GND #PWR0112
U 1 1 603F5F48
P 6600 1750
F 0 "#PWR0112" H 6600 1500 50  0001 C CNN
F 1 "GND" H 6605 1577 50  0000 C CNN
F 2 "" H 6600 1750 50  0001 C CNN
F 3 "" H 6600 1750 50  0001 C CNN
	1    6600 1750
	1    0    0    -1  
$EndComp
Wire Wire Line
	6600 950  6600 1050
$Comp
L Connector_Generic:Conn_02x08_Odd_Even J2
U 1 1 604184E7
P 6300 1250
F 0 "J2" H 6350 1767 50  0000 C CNN
F 1 "Conn_02x08_Odd_Even" H 6350 1676 50  0000 C CNN
F 2 "Connector_PinHeader_2.54mm:PinHeader_2x08_P2.54mm_Vertical" H 6300 1250 50  0001 C CNN
F 3 "~" H 6300 1250 50  0001 C CNN
	1    6300 1250
	1    0    0    -1  
$EndComp
Wire Wire Line
	6100 1550 6100 1650
Wire Wire Line
	6600 1650 6600 1750
Wire Wire Line
	6100 1550 6100 1450
Connection ~ 6100 1550
Connection ~ 6100 950 
Connection ~ 6100 1050
Wire Wire Line
	6100 1050 6100 950 
Connection ~ 6100 1150
Wire Wire Line
	6100 1150 6100 1050
Connection ~ 6100 1250
Wire Wire Line
	6100 1250 6100 1150
Connection ~ 6100 1350
Wire Wire Line
	6100 1350 6100 1250
Connection ~ 6100 1450
Wire Wire Line
	6100 1450 6100 1350
Wire Wire Line
	6600 1050 6600 1150
Connection ~ 6600 1050
Connection ~ 6600 1650
Connection ~ 6600 1150
Wire Wire Line
	6600 1150 6600 1250
Connection ~ 6600 1250
Wire Wire Line
	6600 1250 6600 1350
Connection ~ 6600 1350
Wire Wire Line
	6600 1350 6600 1450
Connection ~ 6600 1450
Wire Wire Line
	6600 1450 6600 1550
Connection ~ 6600 1550
Wire Wire Line
	6600 1550 6600 1650
$EndSCHEMATC
