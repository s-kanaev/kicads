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
Wire Wire Line
	5850 2000 5850 1900
$Comp
L smoke-machine-rescue:NE555-Timer U1
U 1 1 5E0CF761
P 2800 1700
F 0 "U1" H 3100 2050 50  0000 C CNN
F 1 "NE555" H 3050 1350 50  0000 C CNN
F 2 "Package_DIP:DIP-8_W7.62mm" H 2800 1700 50  0001 C CNN
F 3 "http://www.ti.com/lit/ds/symlink/ne555.pdf" H 2800 1700 50  0001 C CNN
	1    2800 1700
	1    0    0    -1  
$EndComp
$Comp
L Device:R_POT RV1
U 1 1 5E0D270E
P 4300 2650
F 0 "RV1" V 4093 2650 50  0000 C CNN
F 1 "100k" V 4184 2650 50  0000 C CNN
F 2 "Potentiometer_THT:Potentiometer_Bourns_3005_Horizontal" H 4300 2650 50  0001 C CNN
F 3 "~" H 4300 2650 50  0001 C CNN
	1    4300 2650
	1    0    0    -1  
$EndComp
$Comp
L Device:C C1
U 1 1 5E0D406B
P 3400 2550
F 0 "C1" H 3515 2596 50  0000 L CNN
F 1 "33n" H 3515 2505 50  0000 L CNN
F 2 "Capacitor_THT:C_Disc_D10.0mm_W2.5mm_P5.00mm" H 3438 2400 50  0001 C CNN
F 3 "~" H 3400 2550 50  0001 C CNN
	1    3400 2550
	1    0    0    -1  
$EndComp
Text Label 3400 1900 0    50   ~ 0
U1_TR
Wire Wire Line
	3400 1900 3300 1900
Text Label 2150 1500 2    50   ~ 0
U1_TR
Wire Wire Line
	2150 1500 2300 1500
Wire Wire Line
	3300 1500 4300 1500
Wire Wire Line
	4300 1500 4300 1550
Wire Wire Line
	4650 1500 4300 1500
Connection ~ 4300 1500
$Comp
L Device:CP C2
U 1 1 5E0FE709
P 950 1200
F 0 "C2" H 1068 1246 50  0000 L CNN
F 1 "1u" H 1068 1155 50  0000 L CNN
F 2 "Capacitor_THT:CP_Radial_D8.0mm_P2.50mm" H 988 1050 50  0001 C CNN
F 3 "~" H 950 1200 50  0001 C CNN
	1    950  1200
	1    0    0    -1  
$EndComp
Wire Wire Line
	950  1500 950  1450
$Comp
L Connector:TestPoint TPpwm1
U 1 1 5E10334A
P 4300 1350
F 0 "TPpwm1" H 4358 1468 50  0000 L CNN
F 1 "TestPoint" H 4358 1377 50  0000 L CNN
F 2 "TestPoint:TestPoint_Pad_4.0x4.0mm" H 4500 1350 50  0001 C CNN
F 3 "~" H 4500 1350 50  0001 C CNN
	1    4300 1350
	1    0    0    -1  
$EndComp
Wire Wire Line
	4300 1350 4300 1500
$Comp
L Connector:TestPoint TPpwmC1
U 1 1 5E10759A
P 4750 2650
F 0 "TPpwmC1" V 4704 2838 50  0000 L CNN
F 1 "TestPoint" V 4795 2838 50  0000 L CNN
F 2 "TestPoint:TestPoint_Pad_4.0x4.0mm" H 4950 2650 50  0001 C CNN
F 3 "~" H 4950 2650 50  0001 C CNN
	1    4750 2650
	0    1    1    0   
$EndComp
$Comp
L Device:R R1
U 1 1 5E12CDC7
P 4800 1500
F 0 "R1" V 4593 1500 50  0000 C CNN
F 1 "100..1k" V 4684 1500 50  0000 C CNN
F 2 "Resistor_THT:R_Axial_DIN0207_L6.3mm_D2.5mm_P10.16mm_Horizontal" V 4730 1500 50  0001 C CNN
F 3 "~" H 4800 1500 50  0001 C CNN
	1    4800 1500
	0    1    1    0   
$EndComp
Wire Wire Line
	3400 1900 3400 2400
Text Label 3450 1700 0    50   ~ 0
U1_DISCH
Wire Wire Line
	3450 1700 3300 1700
$Comp
L Device:R Rrb1
U 1 1 5E150ADC
P 4300 2200
F 0 "Rrb1" H 4370 2246 50  0000 L CNN
F 1 "1k" H 4370 2155 50  0000 L CNN
F 2 "Resistor_THT:R_Axial_DIN0207_L6.3mm_D2.5mm_P10.16mm_Horizontal" V 4230 2200 50  0001 C CNN
F 3 "~" H 4300 2200 50  0001 C CNN
	1    4300 2200
	1    0    0    -1  
$EndComp
$Comp
L Device:R Rlb1
U 1 1 5E1511ED
P 4300 3100
F 0 "Rlb1" H 4370 3146 50  0000 L CNN
F 1 "1k" H 4370 3055 50  0000 L CNN
F 2 "Resistor_THT:R_Axial_DIN0207_L6.3mm_D2.5mm_P10.16mm_Horizontal" V 4230 3100 50  0001 C CNN
F 3 "~" H 4300 3100 50  0001 C CNN
	1    4300 3100
	1    0    0    -1  
$EndComp
Wire Wire Line
	4300 1850 4300 2050
Wire Wire Line
	4300 2500 4300 2350
Wire Wire Line
	4300 2800 4300 2950
$Comp
L Diode:BAT42 D1
U 1 1 5E15581A
P 4300 1700
F 0 "D1" V 4346 1621 50  0000 R CNN
F 1 "BAT42" V 4255 1621 50  0000 R CNN
F 2 "Diode_THT:D_DO-35_SOD27_P7.62mm_Horizontal" H 4300 1525 50  0001 C CNN
F 3 "http://www.vishay.com/docs/85660/bat42.pdf" H 4300 1700 50  0001 C CNN
	1    4300 1700
	0    -1   -1   0   
$EndComp
Wire Wire Line
	4450 2650 4600 2650
Text Label 4650 2450 0    50   ~ 0
U1_TR
Wire Wire Line
	4650 2450 4600 2450
Wire Wire Line
	4600 2450 4600 2650
Connection ~ 4600 2650
Wire Wire Line
	4600 2650 4750 2650
Text Label 4300 3600 1    50   ~ 0
U1_DISCH
Wire Wire Line
	4300 3250 4300 3600
$Comp
L Connector:Conn_01x01_Female Jpwr_gnd2
U 1 1 5E198AF4
P 2050 3700
F 0 "Jpwr_gnd2" V 1988 3612 50  0000 R CNN
F 1 "Conn_01x01_Female" V 1897 3612 50  0000 R CNN
F 2 "Connector_Pin:Pin_D1.1mm_L8.5mm_W2.5mm_FlatFork" H 2050 3700 50  0001 C CNN
F 3 "~" H 2050 3700 50  0001 C CNN
	1    2050 3700
	0    -1   -1   0   
$EndComp
$Comp
L Connector:Conn_01x01_Female Jpwr_vcc_1
U 1 1 5E19E7DB
P 850 2500
F 0 "Jpwr_vcc_1" V 788 2412 50  0000 R CNN
F 1 "Conn_01x01_Female" V 697 2412 50  0000 R CNN
F 2 "Connector_Pin:Pin_D1.1mm_L8.5mm_W2.5mm_FlatFork" H 850 2500 50  0001 C CNN
F 3 "~" H 850 2500 50  0001 C CNN
	1    850  2500
	0    -1   -1   0   
$EndComp
Text Label 850  2950 3    50   ~ 0
PWR_VCC_1
Wire Wire Line
	850  2950 850  2700
$Comp
L Connector:Conn_01x01_Female Jpwr_gnd1
U 1 1 5E1A65EF
P 850 3700
F 0 "Jpwr_gnd1" V 788 3612 50  0000 R CNN
F 1 "Conn_01x01_Female" V 697 3612 50  0000 R CNN
F 2 "Connector_Pin:Pin_D1.1mm_L8.5mm_W2.5mm_FlatFork" H 850 3700 50  0001 C CNN
F 3 "~" H 850 3700 50  0001 C CNN
	1    850  3700
	0    -1   -1   0   
$EndComp
Text Label 850  4100 3    50   ~ 0
PWR_GND_1
Wire Wire Line
	850  3900 850  4100
Text Label 2050 4100 3    50   ~ 0
PWR_GND_2
Wire Wire Line
	2050 4100 2050 3900
Text Label 950  800  0    50   ~ 0
PWR_VCC_1
Text Label 2800 900  0    50   ~ 0
PWR_VCC_1
Text Label 5850 650  0    50   ~ 0
PWR_VCC_1
Wire Wire Line
	950  800  950  950 
Wire Wire Line
	2800 900  2800 1300
Text Label 2800 2450 3    50   ~ 0
PWR_GND_3
Wire Wire Line
	2800 2100 2800 2450
Text Label 3400 2900 3    50   ~ 0
PWR_GND_2
Wire Wire Line
	3400 2900 3400 2700
Text Label 5850 2000 3    50   ~ 0
PWR_GND_1
Text Label 950  1500 3    50   ~ 0
PWR_GND_3
Text Label 1800 1900 2    50   ~ 0
PWR_VCC_1
Wire Wire Line
	1800 1900 2300 1900
NoConn ~ 2300 1700
$Comp
L Diode:1.5KExxA D2
U 1 1 5E375CAA
P 6500 1050
F 0 "D2" V 6454 1129 50  0000 L CNN
F 1 "1.5KExxA" V 6545 1129 50  0000 L CNN
F 2 "Diode_THT:D_DO-201AE_P15.24mm_Horizontal" H 6500 850 50  0001 C CNN
F 3 "https://www.vishay.com/docs/88301/15ke.pdf" H 6450 1050 50  0001 C CNN
	1    6500 1050
	0    1    1    0   
$EndComp
Wire Wire Line
	5850 650  5850 800 
Wire Wire Line
	5850 1250 5850 1300
$Comp
L Diode:1.5KExxA D3
U 1 1 5E384673
P 6500 1650
F 0 "D3" V 6454 1729 50  0000 L CNN
F 1 "1.5KExxA" V 6545 1729 50  0000 L CNN
F 2 "Diode_THT:D_DO-201AE_P15.24mm_Horizontal" H 6500 1450 50  0001 C CNN
F 3 "https://www.vishay.com/docs/88301/15ke.pdf" H 6450 1650 50  0001 C CNN
	1    6500 1650
	0    1    1    0   
$EndComp
Wire Wire Line
	6500 900  6500 800 
Wire Wire Line
	5850 800  6500 800 
Wire Wire Line
	6500 1200 6500 1250
Wire Wire Line
	5850 1250 5950 1250
Connection ~ 6500 1250
Wire Wire Line
	6500 1250 6500 1500
Wire Wire Line
	6500 1800 6500 1900
Wire Wire Line
	6500 1900 5850 1900
Connection ~ 5850 1900
Wire Wire Line
	5850 1900 5850 1700
$Comp
L Device:C C3
U 1 1 5E38C05F
P 1450 1200
F 0 "C3" H 1565 1246 50  0000 L CNN
F 1 "0.1u" H 1565 1155 50  0000 L CNN
F 2 "Capacitor_THT:C_Disc_D10.0mm_W2.5mm_P5.00mm" H 1488 1050 50  0001 C CNN
F 3 "~" H 1450 1200 50  0001 C CNN
	1    1450 1200
	1    0    0    -1  
$EndComp
Wire Wire Line
	1450 1050 1450 950 
Wire Wire Line
	1450 950  950  950 
Connection ~ 950  950 
Wire Wire Line
	950  950  950  1050
Wire Wire Line
	1450 1350 1450 1450
Wire Wire Line
	1450 1450 950  1450
Connection ~ 950  1450
Wire Wire Line
	950  1450 950  1350
$Comp
L IRFB3004PbF:IRFB3004PbF Q1
U 1 1 5E374851
P 5850 1500
F 0 "Q1" H 6008 1546 50  0000 L CNN
F 1 "IRFB3004PbF" H 6008 1455 50  0000 L CNN
F 2 "smoke-machine:TO-220-3-4_Horizontal_TabDown" H 6000 1350 50  0001 C CNN
F 3 "" H 6000 1350 50  0001 C CNN
	1    5850 1500
	1    0    0    -1  
$EndComp
Wire Wire Line
	5950 1300 5950 1250
Connection ~ 5950 1250
Wire Wire Line
	5950 1250 6500 1250
$Comp
L Connector:Conn_01x01_Female Jpwr_gnd3
U 1 1 5E38430C
P 3250 3700
F 0 "Jpwr_gnd3" V 3188 3612 50  0000 R CNN
F 1 "Conn_01x01_Female" V 3097 3612 50  0000 R CNN
F 2 "Connector_Pin:Pin_D1.1mm_L8.5mm_W2.5mm_FlatFork" H 3250 3700 50  0001 C CNN
F 3 "~" H 3250 3700 50  0001 C CNN
	1    3250 3700
	0    -1   -1   0   
$EndComp
Text Label 3250 4100 3    50   ~ 0
PWR_GND_3
Wire Wire Line
	3250 4100 3250 3900
Connection ~ 5850 1250
Wire Wire Line
	5850 1050 5850 1250
Connection ~ 5850 800 
Wire Wire Line
	5850 800  5850 950 
$Comp
L Connector:Conn_01x02_Female J1
U 1 1 5E07B595
P 5650 1050
F 0 "J1" H 5678 1026 50  0000 L CNN
F 1 "Conn_01x02_Female" H 5678 935 50  0000 L CNN
F 2 "TerminalBlock:TerminalBlock_bornier-2_P5.08mm" H 5650 1050 50  0001 C CNN
F 3 "~" H 5650 1050 50  0001 C CNN
	1    5650 1050
	-1   0    0    1   
$EndComp
Wire Wire Line
	4950 1500 5550 1500
Text Label 5100 1500 0    50   ~ 0
PWM_OUT
$EndSCHEMATC
