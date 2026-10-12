PROGRAM: ACRYLAMIDE (AA) and GLYCIDAMIDE (GA) MODEL (AA_GA.CSL)
!'4-Compartment model for AA administration and GA production'
!'original by CR KIRMAN 03/01'
!'parameters initialized for 250 g rat' 
!revised by Lisa M. Sweeney October 21, 2008
!expand number of tissues, add urinary elimination of glycidamide
!revise partition coefficients
!aa-val and ga-val adducts

INITIAL 

CONSTANT   QPC = 14.  !$'Alveolar ventilation rate (L/hr-kg^0.74)'
CONSTANT   QCC = 14.  !$'Cardiac output (L/hr-kg^0.74)'
CONSTANT   QLC = 0.18  !$'Blood flow to liver'
CONSTANT   QBrC = 0.02  !Blood flow to brain
CONSTANT   QKC = 0.13 !Blood flow to kidney
CONSTANT   QSC = 0.34 !Blood flow to slowly perfused tissue (skin and muscle)
CONSTANT   QFC = 0.07 !Blood flow to adipose tissue (fat)
QRC = 1-(QLC+QBRC+QKC+QSC+QFC)
CONSTANT   PAS1 = 1.0  !Permeability-area cross product for AA in spt as frac of tissue blood flow

CONSTANT   BW = 0.25  !$'Body weight (kg)'
CONSTANT   VUC = 0.05 !fraction unperfused tissue
CONSTANT   VLC = 0.037 !$'Fraction liver tissue'
CONSTANT   VBRC = 0.006 !Fraction brain/nervous system tissues
CONSTANT   VKC = 0.0073 !Fraction kidney tissue
CONSTANT   VSC = 0.59   !Fraction slowly perfused tissue (skin and muscle)
CONSTANT   VSBC = 0.02  !Fraction of bw as blood in SPT
CONSTANT   VFC = 0.087  !Fraction adipose tissue (fat)
CONSTANT   VBC = 0.074 !$'Fraction of blood vol'
VRC = 1-(VUC+VLC+VBRC+VKC+VSC+VFC+VBC)  !fraction richly perfused tissue

CONSTANT  VDAAMAC = 1.0 !Volume of distribution of AAMA (L/kg)
CONSTANT  VDGAMAC = 1.0 !Volume of distribution of GAMA (L/kg)
CONSTANT  VDGAOHC = 1.0 !Volume of distribution of GA hydrolysis product (L/kg)

CONSTANT   VABC = 0.35  !$'Fraction arterial blood'
VVBC = 1 - VABC  !$'Fraction venous blood'
CONSTANT   FBC = 0.44  !$'Fraction blood cells'
FBS = 1-FBC  !$'Fraction serum'

CONSTANT   PL1 = 1.7    !$'liver/blood partition AA'
CONSTANT   Pbr1 = 0.75  !$'brain/blood partition AA'
CONSTANT   pk1 = 1.7     !Kidney/blood partition AA
CONSTANT   PR1 = 1.7     !rpt/blood partition AA
CONSTANT   PS1 = 0.69    !spt/blood partition AA
CONSTANT   PF1 =0.1      !fat/blood partition AA
CONSTANT   PB1 = 31000000  !$'blood/air partition AA'
CONSTANT   PB2 = 98000000  !$'blood/air partition GA'
CONSTANT   PL2 = 0.5   !$'liver/blood partition GA'
CONSTANT   Pbr2 = 0.5  !$'brain/blood partition GA'
CONSTANT   pk2 = 0.5     !Kidney/blood partition GA
CONSTANT   PR2 = 0.5     !rpt/blood partition GA
CONSTANT   PS2 = 1.0    !spt/blood partition GA
CONSTANT   PF2 =0.1      !fat/blood partition GA

CONSTANT   MW1 = 71.    
CONSTANT   MW2 = 87.    !'molecular weight (g/mol)'

CONSTANT   VMAXC1 = 9.  !metabolism of AA to GA mg/hr-kg^0.7
CONSTANT   KMC1 = 10.	!Km for metabolism of AA to GA mg/L    

CONSTANT   VMAXC2 = 0.3 !hyrolysis of GA (mg/hr-kg^0.7)'
CONSTANT   KMC2 = 100   !km for hydrolysis of GA mg/L

CONSTANT   VMAXgC1 = 1  !Vmax for GSH-AA conjugation mg/hr-kg^0.7
CONSTANT   KMGC1 = 100  !Km with respect to AA for GSH conjugation mg/L
CONSTANT   KMGG = 0.1   !KM with respect to GSH for AA or GA conjugation with GSH mM
CONSTANT   VMAXgC2 = 1  !Vmax for GSH-GA conjugation mg/hr-kg^0.7
CONSTANT   KMGC2 = 100  !Km with respect to GA for GSH conjugation mg/L

CONSTANT   KUC2 = 0.00062 !urinary elimination of GA--fraction of kidney blood flow

CONSTANT   KDSIN1=0.2    !transfer rate to deep tissue compartment of SPT 1/hr
CONSTANT   KDSOUT1=0.1  !transfer rate out of deept tissue compartment of SPT 1/hr

CONSTANT   KpbL1 = 0.1    !protein binding in liver
CONSTANT   kpbbr1 = 0.08  !protein binding in brain
CONSTANT   kpbk1 = 0.08   !protein binding in kidney
CONSTANT   kpbr1 = 0.08   !protein binding in rpt
CONSTANT   kpbs1 = 0.08   !protein binding in spt
CONSTANT   kpbf1 = 0.08   !protein binding in fat
CONSTANT   kpbrb1 = 0.01  !protein binding in rbc (hemoglobin binding)
CONSTANT   kpbpl1 = 0.01   !protein binding in plasma 
CONSTANT   reactratio = 2.0  !reactivity ratio with protein, AA/GA
   KpbL2 = kpbl1/reactratio    !protein binding in liver
   kpbbr2 = kpbbr1/reactratio  !protein binding in brain
   kpbk2 = kpbk1/reactratio   !protein binding in kidney
   kpbr2 = kpbr1/reactratio   !protein binding in rpt
   kpbs2 = kpbs1/reactratio   !protein binding in spt
   kpbf2 = kpbf1/reactratio   !protein binding in fat
   kpbrb2 = kpbrb1/reactratio  !protein binding in rbc (hemoglobin binding)
   kpbpl2 = kpbpl1/reactratio  !protein binding in plasma 

CONSTANT   KPTL = 0.015   !'protein turnover rate in liver'
CONSTANT   KPTBR = 0.013   !'protein turnover rate in brain'
CONSTANT   KPTK = 0.013   !'protein turnover rate in kidney'
CONSTANT   KPTR = 0.013   !'protein turnover rate in rpt'
CONSTANT   KPTS = 0.0039   !'protein turnover rate in spt'
CONSTANT   KPTF = 0.0039   !'protein turnover rate in fat'
CONSTANT   KPTRB = 0.0039   !'protein turnover rate in rbc'
CONSTANT   KPTPL = 0.0039   !'protein turnover rate in plasma'

CONSTANT   KFORMAAVAL = 7500  !fmol AA-val/mg globin per mM AA-hr
CONSTANT   KFORMGAVAL = 34000   !fmol GA-val/mg globin per mM GA-hr

CONSTANT   KREMAAVAL = 0.00231  !removal of AA-val adducts from RBC per hr
CONSTANT   KREMGAVAL = 0.00231  !removal of GA-val adducts from RBC per hr


CONSTANT   KUAAMA = 1.0 ! urinary elimination rate of AAMA (1/hr)
CONSTANT   KUGAMA = 1.0 ! Urinary elimination rate of AAGA (1/hr)
CONSTANT   KUGAOH = 1.0 ! Urinary elimination rate of GA hydrolysis product (1/hr)


!'Glutathione parameters, DSouza et al., 1988'

CONSTANT   KGSHL = 0.14       !$'Liver GSH turnover (/hr)'
CONSTANT   GSHL0 = 7.0  !$'Initial GSH (mmol/l)'
    AGSH0 = GSHL0*VLC*BW
    MINPGSHL=1.0
CONSTANT   DEPLETIONON = 0.0 !Liver GSH depletion turned on, zero for off and 1 for onn

CONSTANT   IVDOSE = 0.  !$'IV dose AA (mg/kg)' 
           IVAMT = IVDOSE*BW/MW1
           IVR = IVAMT/TINF
CONSTANT   IVDOSE2 = 0.  !$'IV dose GA (mg/kg)' 
           IVAMT2 = IVDOSE2*BW/MW2
           IVR2 = IVAMT2/TINF
CONSTANT   IPDOSE = 0. !$'IP AA dose (mg/kg)' 
           IPAMT = IPDOSE*BW/MW1
CONSTANT   IPDOSE2 = 0. !$'IP GA dose (mg/kg)' 
           IPAMT2 = IPDOSE2*BW/MW2
CONSTANT   CONC = 0.    !$'Inhaled concentration (ppm) AA'
      
!'Timing commands'
CONSTANT   TSTOP = 24      !$'Length of experiment (hrs)'
CONSTANT   TCHNG1 = 6      !$'Length of exposure (hrs)'
CONSTANT   TCHNG2 = 120.0  !$'allows for 5 day/week exposure'
CONSTANT   TINF = .003     !$'Length of IV infusion (hrs)'
CONSTANT   TING = 0.5      ! ingestion time for diet
   
!'Periodic Drinking Water Exposure Section:'
 ! 'Assume T=0 is 7 am in morning for reference (T+4 = 11 am etc.)'
     INTEGER I $ I = 1      !$'Counter for Drinking Arrays'
    ARRAY DRTIME(6)      !$'Store Drinking Times'
     CONSTANT DRTIME = 1.0,    5.0,    9.0,    13.0,    17.0,    21.0
     ARRAY DRPCT(6) !$'Store Drinking %'
     CONSTANT DRPCT = 0.02,   0.02,   0.02,   0.313,   0.313,   0.313
     CONSTANT DRCONC = 0.0  !$'Conc of Parent in water (mg/liter)'
     DRVOL = 0.102*BW**0.7  !$'Assume 70 kg man & 2 L/d'
     TAMT = DRVOL*DRCONC/MW1    !$'Total amount from water (mmol)'
     DRAMT = 0
     NEWDAY=0

     CONSTANT  ODOSE = 0.0     !$'Cont oral dose AA (mg/kg-day)'
     OAMT = ODOSE*BW/MW1   !$'Daily amount (mmol)'
     ORATE = OAMT/24 !$'Cont oral rate (mmol/hr)'

     CONSTANT  ODOSE2 = 0.0     !'Cont oral dose GA (mg/kg-day)'
     OAMT2 = ODOSE2*BW/MW2   !'Daily amount (mmol)'
     ORATE2 = OAMT2/24 !'Cont oral rate (mmol/hr)'

     CONSTANT DDOSE = 0.      !Dietary dose AA (mg/kg/day)
     DRATE0 = DDOSE*BW/MW1/TING   !Ingestion rate during feeding (mmol/hr)

     CONSTANT GDOSE=0.0	    !$'gavage dose AA (mg/kg)'
     GAMT = GDOSE*BW/MW1    !$'gavage amount AA (mmol)'

     CONSTANT GDOSE2=0.0	    !$'gavage dose AA (mg/kg)'
     GAMT2 = GDOSE2*BW/MW2    !$'gavage amount AA (mmol)'

     CONSTANT KA = 3.       !$'AA Oral or ip uptake rate (/hr)'
     CONSTANT KA2 = 3.      !$'GA Oral or ip uptake rate (/hr)'

!'Scaled parameters'
   
       QC = QCC*BW**0.74
       QP = QPC*BW**0.74
       QL = QLC*QC
       QBR = QBRC*QC
       QK = QKC*QC
       QR = QRC*QC
       QS = QSC*QC
       QF = QFC*QC

       VL = VLC*BW
       VBR = VBRC*BW
       VK = VKC*BW
       VR = VRC*BW
       VS = VSC*BW
       VSB = VSBC*BW
       VF = VFC*BW
       VB = VBC*BW
            VAB = VB*VABC
            VVB = VB*VVBC
                VABRBC = VAB*FBC
                VABSER = VAB*FBS
                  VVBRBC = VVB*FBC
                VVBSER = VVB*FBS

    VMAX1 = VMAXC1/MW1*BW**0.7     !$'Liver P450 AA to GA, mmol/hr'
    VMAX2 = VMAXC2/MW2*BW**0.7    !$'Liver GA Hydrolysis, mmol/hr'
      KM1 = KMC1/MW1
      KM2 = KMC2/MW2
    VMAXG1 = VMAXGC1/MW1*BW**0.7     !$'Liver AA-GSH rate'
    VMAXG2 = VMAXGC2/MW2*BW**0.7     !$'Liver GA-GSH rate'
     KMG1 = KMGC1/MW1  !Km with respect to AA for GSH conjugation mM
     KMG2 = KMGC2/MW2  !Km with respect to GA for GSH conjugation mM

   KU2 = KUC2*QK    !urinary elimination of GA from the kidney
   vdaama = vdaamac*bw    !volume of distribution for AA-gsh conjugate
   vdgama = vdgamac*bw
   vdgaoh = vdgaohc*bw



END     !$'End of initial' 
   
DYNAMIC 
 
ALGORITHM IALG = 2  !$'Gear method for stiff systems'
       
DERIVATIVE  

!'------------------------------ACRYLAMIDE--------------------------------'
!      'CI = Concentration in inhaled air (mmol/L)'
   CIZONE = PULSE(0.0,24.0,TCHNG1) * PULSE(0.0,168.0,TCHNG2)
       CI = CIZONE*CONC/24450.       
   
!      'AI = Amount inhaled (mmol)'
      RAI = QP*CI
       AI = INTEG(RAI,0.)

!     'CAL1 = Concentration in arterial lung blood (mmol/L)'
      CAL1 = (QC*CVB1+QP*CI)/(QC+(QP/PB1))

 !     'AX1 = Amount exhaled (mmol)' 
       CX1 = CAL1/PB1
      RAX1 = QP*CX1 
       AX1 = INTEG(RAX1,0.)
  
  !'CA1 = Conc of free aa in arterial blood (mmol/l)'
      RAB1 = (QC*CAL1)-(QC*CA1) - rpbarb1 - rpbapl1
       AB1 = INTEG(RAB1,0.)
       CA1 = AB1/VAB

!protein binding in arterial blood
  rpbarb1 = ca1*kpbrb1*vabrbc   !covalent binding in arterial RBC
  rpbapl1 =ca1*kpbpl1*vabser   !covalent binding in arterial plasma

!amount AA bound in arterial blood
  relimpbarb1 = apbarb1*kptrb !rate of elimination of AA bound to protein in arterial blood rbc
  rapbarb1 = rpbarb1 - relimpbarb1 + QC*(cpbvrb1-cpbarb1)
   apbarb1 = INTEG(RAPBARB1, 0.)  !Amount of AA bound to protein in arterial blood rbc
   cpbarb1 = apbarb1/vabrbc
   
  relimpbapl1 = apbapl1*kptpl !rate of elimination of AA bound to protein in arterial blood serum
  rapbapl1 = rpbapl1 - relimpbapl1  + QC*(cpbvpl1-cpbapl1)
   apbapl1 = INTEG(RAPBApl1, 0.)  !Amount of AA bound to protein in arterial blood serum
   cpbapl1 = apbapl1/vabser

 !     'ABr1 = Amount in brain (mmol)' 
      RABr1 = QBr*(CA1-CVBr1) - RpbBr1
       ABr1 = INTEG(RABr1,0.)
       CBr1 = ABr1/VBr
      CVBr1 = CBr1/PBr1
    AUCBr1 = INTEG(CBr1, 0.)
    RpbBr1 = KpbBr1*abr1
    RelimpbBr1 = apbBr1*KPTbr
     rApbBr1 = rpbbr1 - relimpbBr1
      apbBr1 = INTEG(RapbBr1, 0.)

  !    'Ar1 = Amount in richly perfused tissue (mmol)' 
      RAr1 = Qr*(CA1-CVr1) - Rpbr1
       Ar1 = INTEG(RAr1,0.)
       Cr1 = Ar1/Vr
      CVr1 = Cr1/Pr1
    AUCr1 = INTEG(Cr1, 0.)
    Rpbr1 = Kpbr1*Ar1
    RelimpbR1 = apbR1*KPTr
     rApbr1 = rpbr1 - relimpbr1
      apbr1 = INTEG(Rapbr1, 0.)

   !   'As1 = Amount in slowly perfused tissue (mmol)' 
   !   RASB1 = Qs*(CA1-CVs1) + QS*PAS1*(CS1/PS1-CVS1)
    !   ASB1 = INTEG(RASB1, 0.)
     !  CVS1 = ASB1/VSB
      
     ! RAS1 = QS*PAS1*(CVS1-CS1/PS1)- Rpbs1 
      RAS1 = Qs*(CA1-CVs1)- Rpbs1
!+sptdeepout -sptdeepin
        ! Deep tissue compartment in SPT
!      sptdeepout= KDSOUT1*ASD1
!      sptdeepin = KDSIN1*AS1
!      RASD1 = sptdeepin-sptdeepout
!       ASD1 = INTEG(RASD1, 0.)
       As1 = INTEG(RAs1,0.)
!      Cs1 = (AS1+ASD1)/VS
    CS1 = AS1/VS
    CVS1 = CS1/PS1
    AUCs1 = INTEG(Cs1, 0.)
    Rpbs1 = Kpbs1*as1
    Relimpbs1 = apbs1*KPTs
     rApbs1 = rpbs1 - relimpbs1
      apbs1 = INTEG(Rapbs1, 0.)

    !  'Af1 = Amount in adipose tissue (fat)(mmol)' 
      RAf1 = Qf*(CA1-CVf1) - Rpbf1
       Af1 = INTEG(RAf1,0.)
       Cf1 = Af1/Vf
      CVf1 = Cf1/Pf1
    AUCf1 = INTEG(Cf1, 0.)
    Rpbf1 = Kpbf1*af1
    Relimpbf1 = apbf1*KPTf
     rApbf1 = rpbf1 - relimpbf1
      apbf1 = INTEG(Rapbf1, 0.)

     ! 'Ak1 = Amount in kidney(mmol)' 
      RAk1 = Qk*(CA1-CVk1) - Rpbk1
       Ak1 = INTEG(RAk1,0.)
       Ck1 = Ak1/Vk
      CVk1 = Ck1/Pk1
    AUCk1 = INTEG(Ck1, 0.)
    Rpbk1 = Kpbk1*ak1
    Relimpbk1 = apbk1*KPTk
     rApbk1 = rpbk1 - relimpbk1
      apbk1 = INTEG(Rapbk1, 0.)

    !'Stomach Compartment for Drinking Water and diet Inputs'
 DRATE = DRATE0*(1-STEP(TING))
  DAMT = INTEG(DRATE, 0.)
     RSTOM = -KA*STOM + ORATE + DRATE
      STOM = INTEG(RSTOM, 0.0) +  GAMT

!'GSHL levels after normal turnover and AA conjugation'
    RAMGSHL = (AGSH0-AMGL)*KGSHL-(RGST1+RGST2)*DEPLETIONON
     AMGL = INTEG(RAMGSHL,AGSH0)
       GSHL = (AMGL/VL)
       PGSHL = GSHL/GSHL0
       IF (PGSHL.LT.MINPGSHL) MINPGSHL=PGSHL

 !   'AL1 = Amount in liver tissue (mmol)'
      RAL1 = QL*CA1 + KA*STOM - QL*CVL1 - RP450 - RGST1 - rpbl1
       AL1 = INTEG(RAL1,0.)
      CVL1 = AL1/(VL*PL1)
       CL1 = AL1/VL
       rpbl1 = kpbl1*al1
    Relimpbl1 = apbL1*KPTL
      RapbL1 = rpbL1 - relimpbl1
     apbL1 = INTEG(Rapbl1, 0.)
 
 !'AP450 = Amount metabolized by P450 (mmol)' 
      RP450 = (VMAX1*CVL1)/(KM1+CVL1)
     AP450 = INTEG(RP450,0.)

 !'AGST1 = Amount metabolized by GST (mmol)'
      RGST1 = VMAXG1*CVL1*GSHL/(CVL1 + KMG1)/(GSHL + KMGG)
     AGST1 = INTEG(RGST1,0.)

 !'IV1 = Intravenous infusion rate (mmol/hr)'
       IV1 = IVR*(1.-step(tinf))
      
 ! 'CV1 = Mixed venous blood concentration (mmol/L)'
       CV1 = (QL*CVL1+QS*CVS1+QBR*CVBR1+QK*CVK1+QR*CVR1+QF*CVF1+IV1)/QC 

  !'IP dosing'
      RPER1 = -KA*PER1
      PER1 = INTEG(RPER1, 0.) + IPAMT

 !'CVB1 = Mixed ven. conc after binding'
      RVB1 = (QC*CV1)-(QC*CVB1) - rpbvrb1 - rpbvpl1+ KA*PER1
       VB1 = INTEG(RVB1,0.)
      CVB1 = VB1/VVB
     AUCVB1 = INTEG(CVB1,0.)

!protein binding in venous blood
  rpbvrb1 = cv1*kpbrb1*vvbrbc   !covalent binding in venous RBC
  rpbvpl1 =cv1*kpbpl1*vvbser   !covalent binding in venous plasma


!amount AA bound in venous blood
  relimpbvrb1 = apbvrb1*kptrb !rate of elimination of AA bound to protein in venous blood rbc
  rapbvrb1 = rpbvrb1 - relimpbvrb1 + QC*(cpbarb1-cpbvrb1)
   apbvrb1 = INTEG(RAPBvRB1, 0.)  !Amount of AA bound to protein in venous blood rbc
   cpbvrb1 = apbvrb1/vvbrbc

!terminal valine hemoglobin adducts
  RAAVAL = KFORMAAVAL * CVB1 - KREMAAVAL*AAVAL
   AAVAL = INTEG(RAAVAL, 0.)  !fmol adducts per mg globin

   
  relimpbvpl1 = apbvpl1*kptpl !rate of elimination of AA bound to protein in arterial blood serum
  rapbvpl1 = rpbvpl1 - relimpbvpl1  + QC*(cpbapl1-cpbvpl1)
   apbvpl1 = INTEG(RAPBvpl1, 0.)  !Amount of AA bound to protein in arterial blood serum
   cpbvpl1 = apbvpl1/vvbser
 
!'--------------------------------GLYCIDAMIDE----------------------------'
         
 !  '  'CAL2 = Concentration GA in arterial lung blood (mmol/L)'
      CAL2 = (QC*CVB2)/(QC+(QP/PB2))
  
  !    'CA2 = Conc of free ga in arterial blood (mmol/l)'
      RAB2 = (QC*CAL2)-(QC*CA2) - rpbarb2 - rpbapl2
       AB2 = INTEG(RAB2,0.)
       CA2 = AB2/VAB

!protein binding in arterial blood
  rpbarb2 = ca2*kpbrb2*vabrbc   !covalent binding in arterial RBC
  rpbapl2 =ca2*kpbpl2*vabser   !covalent binding in arterial plasma

!amount GA bound in arterial blood
  relimpbarb2 = apbarb2*kptrb !rate of elimination of GA bound to protein in arterial blood rbc
  rapbarb2 = rpbarb2 - relimpbarb2 + QC*(cpbvrb2-cpbarb2)
   apbarb2 = INTEG(RAPBARB2, 0.)  !Amount of GA bound to protein in arterial blood rbc
   cpbarb2 = apbarb2/vabrbc
   
  relimpbapl2 = apbapl2*kptpl !rate of elimination of GA bound to protein in arterial blood serum
  rapbapl2 = rpbapl2 - relimpbapl2 + QC*(cpbvpl2-cpbapl2)
   apbapl2 = INTEG(RAPBApl2, 0.)  !Amount of AA bound to protein in arterial blood serum
   cpbapl2 = apbapl2/vabser


 !   'ABr2 = Amount in brain (mmol)' 
      RABr2 = QBr*(CA2-CVBr2) - RpbBr2
       ABr2 = INTEG(RABr2,0.)
       CBr2 = ABr2/VBr
      CVBr2 = CBr2/PBr2
    AUCBr2 = INTEG(CBr2, 0.)
    RpbBr2 = KpbBr2*abr2
    RelimpbBr2 = apbBr2*KPTbr
     rApbBr2 = rpbbr2 - relimpbBr2
      apbBr2 = INTEG(RapbBr2, 0.)

 !     'Ar2 = Amount in richly perfused tissue (mmol)' 
      RAr2 = Qr*(CA2-CVr2) - Rpbr2
       Ar2 = INTEG(RAr2,0.)
       Cr2 = Ar2/Vr
      CVr2 = Cr2/Pr2
    AUCr2 = INTEG(Cr2, 0.)
    Rpbr2 = Kpbr2*Ar2
    RelimpbR2 = apbR2*KPTr
     rApbr2 = rpbr2 - relimpbr2
      apbr2 = INTEG(Rapbr2, 0.)

!      'As2 = Amount in slowly perfused tissue (mmol)' 
      RAs2 = Qs*(CA2-CVs2) - Rpbs2
       As2 = INTEG(RAs2,0.)
       Cs2 = As2/Vs
      CVs2 = Cs2/Ps2
    AUCs2 = INTEG(Cs2, 0.)
    Rpbs2 = Kpbs2*as2
    Relimpbs2 = apbs2*KPTs
     rApbs2 = rpbs2 - relimpbs2
      apbs2 = INTEG(Rapbs2, 0.)

!      'Af2 = Amount in adipose tissue (fat)(mmol)' 
      RAf2 = Qf*(CA2-CVf2) - Rpbf2
       Af2 = INTEG(RAf2,0.)
       Cf2 = Af2/Vf
      CVf2 = Cf2/Pf2
    AUCf2 = INTEG(Cf2, 0.)
    Rpbf2 = Kpbf2*af2
    Relimpbf2 = apbf2*KPTf
     rApbf2 = rpbf2 - relimpbf2
      apbf2 = INTEG(Rapbf2, 0.)

 !     'Ak2 = Amount in kidney(mmol)' 
      RAk2 = Qk*(CA2-CVk2) - Rpbk2 - ru2
       Ak2 = INTEG(RAk2,0.)
       Ck2 = Ak2/Vk
      CVk2 = Ck2/Pk2
    AUCk2 = INTEG(Ck2, 0.)
    Rpbk2 = Kpbk2*ak2
    Relimpbk2 = apbk2*KPTk
     rApbk2 = rpbk2 - relimpbk2
      apbk2 = INTEG(Rapbk2, 0.)
      ru2 = ku2*ca2
      AU2 = INTEG(RU2, 0.)

     RSTOM2 = -KA2*STOM2 + ORATE2
      STOM2 = INTEG(RSTOM2, 0.0) + GAMT2

!    'AL2 = Amount in liver tissue (mmol)'
      RAL2 = QL*CA2 + KA2*STOM2 - QL*CVL2 + RP450 - RGST2 - REH - rpbl2
       AL2 = INTEG(RAL2,0.)
      CVL2 = AL2/(VL*PL2)
       CL2 = AL2/VL
       rpbl2 = kpbl2*al2
    Relimpbl2 = apbL2*KPTL
      RapbL2 = rpbL2 - relimpbl2
     apbL2 = INTEG(Rapbl2, 0.)
 
! ' 'AGST2 = Amount GA metabolized by GST (mmol)'
      RGST2 = VMAXG2*CVL2*GSHL/(CVL2 + KMG2)/(GSHL + KMGG)
     AGST2 = INTEG(RGST2,0.)

!      'AEH = Amount GA metabolized by EH (mmol)'
        REH = (VMAX2*CVL2)/(KM2+CVL2)
       AEH = INTEG(REH,0.)

!      'IP dosing'
      RPER2 = -KA2*PER2
      PER2 = INTEG(RPER2, 0.) + IPAMT2

! 'IV2 = Intravenous infusion rate (mmol/hr)'
       IV2 = IVR2*(1.-step(tinf))
      
!  'CV2 = Mixed venous blood concentration (mmol/L)'
       CV2 = (QL*CVL2+QS*CVS2+QBR*CVBR2+QK*CVK2+QR*CVR2+QF*CVF2+IV2)/QC 
   
 !'CVB2 = Mixed ven. conc after binding'
      RVB2 = (QC*CV2)-(QC*CVB2) - rpbvrb2 - rpbvpl2+ KA2*PER2
       VB2 = INTEG(RVB2,0.)
      CVB2 = VB2/VVB
      AUCVB2 = INTEG(CVB2,0.)

!protein binding in venous blood
  rpbvrb2 = cv2*kpbrb2*vvbrbc   !covalent binding in venous RBC
  rpbvpl2 =cv2*kpbpl2*vvbser   !covalent binding in venous plasma

!terminal valine hemoglobin adducts
  RGAVAL = KFORMGAVAL*CVB2 - KREMGAVAL*GAVAL
   GAVAL = INTEG(RGAVAL, 0.)  !fmol adducts per mg globin


!amount GA bound in venous blood
  relimpbvrb2 = apbvrb2*kptrb !rate of elimination of GA bound to protein in venous blood rbc
  rapbvrb2 = rpbvrb2 - relimpbvrb2 + QC*(cpbarb2-cpbvrb2)
   apbvrb2 = INTEG(RAPBvRB2, 0.)  !Amount of GA bound to protein in venous blood rbc
   cpbvrb2 = apbvrb2/vvbrbc
   
  relimpbvpl2 = apbvpl2*kptpl !rate of elimination of GA bound to protein in arterial blood serum
  rapbvpl2 = rpbvpl2 - relimpbvpl2  + QC*(cpbapl2-cpbvpl2)
   apbvpl2 = INTEG(RAPBvpl2, 0.)  !Amount of GA bound to protein in arterial blood serum
   cpbvpl2 = apbvpl2/vvbser

!aama

raama = rgst1 - raamau
aama = INTEG(raama, 0.)   !mmol
raamau = aama*kuaama  !rate of elimination into urine
aamau = INTEG(RAAMAU, 0.)     !mmol

rgama = rgst2 - rgamau
gama = INTEG(rgama, 0.)   !mmol
rgamau = gama*kugama  !rate of elimination into urine
gamau = INTEG(RgAMAU, 0.)     !mmol

rgaoh = reh-rgaohu
gaoh = INTEG(RGAOH, 0.)   !mmol
rgaohu = gaoh*kugaoh  !rate of elimination into urine
gaohu = INTEG(RGAOHU, 0.)

#!Combined Dose Metrics for AA + GA '
    TOTL = (AL1+APBL1+AL2+APBL2)/VL
    BNDL = (APBL1+APBL2)/VL
    TOTS = (AS1+APBS1+AS2+APBS2)/VS
    BNDS = (APBS1+APBS2)/VS
    TOTBR = (ABR1+APBBR1+ABR2+APBBR2)/VBR
    TOTF = (AF1+APBF1+AF2+APBF2)/VF
    TOTK = (AK1+APBK1+AK2+APBK2)/VK
    TOTVSER = cpbvpl1+ cpbvpl2 + ca1 + ca2 
    TOTVRBC = cpbvrb1+ cpbvrb2 + ca1 + ca2 
    TOTV = (VB1+VB2+APBVRB1+APBVRB2+APBVPL1+APBVPL2)/VVB
    OTHERMET = (AAMA+GAMA+GAOH)/(BW*(1.0-VUC))
    TOTCARC=(AS1+APBS1+AS2+APBS2+AR1+APBR1+AR2+APBR2+AF1+APBF1+AF2+APBF2)/(VS+VR+VF)

!amount eliminated (into urine?) from protein turnover
relimpb1 = relimpbarb1+relimpbapl1+RelimpbBr1+RelimpbL1+Relimpbk1+ Relimpbs1+ Relimpbr1+ Relimpbf1+ relimpbvrb1+relimpbvpl1
relimpb2 = relimpbarb2+relimpbapl2+RelimpbBr2+RelimpbL2+Relimpbk2+ Relimpbs2+ Relimpbr2+ Relimpbf2+ relimpbvrb2+relimpbvpl2
relimpb = relimpb1 + relimpb2
aelimpb = INTEG(RELIMPB, 0.)


!'Mass check '
molesin1 = OAMT+TAMT+AI+IVAMT+IPAMT+GAMT+DAMT
molesin2 = OAMT2 + IVAMT2 + IPAMT2 + GAMT2
molesin = molesin1+ molesin2

bodyaa = AL1+ABR1+AK1+AR1+AS1+AF1+AB1+VB1
bodyga = AL2+ABR2+AK2+AR2+AS2+AF2+AB2+VB2
apbtissue =apbbr1+apbs1+apbr1+apbl1+apbk1+apbf1+apbbr2+apbs2+apbr2+apbl2+apbk2+apbf2
apbblood = apbarb1+apbapl1+apbvrb1+apbvpl1+apbarb2+apbapl2+apbvrb2+apbvpl2
apb=apbtissue+apbblood
bodyothermet= aama+gama+gaoh
urine = aamau+gamau +gaohu+au2+aelimpb

!'Percent AA Through GST Pathway'

PGST1 = AGST1/(OAMT+TAMT+AI+IVAMT+IPAMT+GAMT+1.e-12)*100

PGST2 = AGST2/(OAMT+TAMT+AI+IVAMT+IPAMT+GAMT+1.e-12)*100
!'Percent AA Oxidized'
PP450 = AP450/(OAMT+TAMT+AI+IVAMT+IPAMT+GAMT+1.e-12)*100

!'Percent AA Oxidized by P450 and Hydrolyzed by EH'
PEH = AEH/(OAMT+TAMT+AI+IVAMT+IPAMT+GAMT+1.e-12)*100


TERMT(T.GE.TSTOP) 
      
END      ! $'End of derivative'
END      ! $'End of dynamic'

TERMINAL

END       !$'End of terminal'  
END       !$'End of program'
