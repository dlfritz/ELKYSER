      !
      !##################################################################
      !#								#
      !#		Elektrolyse Simulation				#
      !#		       09.07.2012				#
      !#								#
      !#	Written By: Dr. David L. Fritz				#
      !#	Contact: d.fritz@fz-juelich.de				#
      !#								#
      !#	Forschungszentrum Juelich				#
      !#	Instituet fuer Energie- und Klimaforschung (IEK)	#
      !#	IEK-3: Brennstoffzellen					#
      !#	Wilhelm-Johnen-Strasse					#
      !#	52428 Juelich						#
      !#								#
      !#================================================================#
      !#								#
      !# Definition of Parameters:					#
      !# ========================					#
      !#								#
      !# RGC	-> Universal Gas Constant [J/(mol-K)]			#
      !#								#
      !# F	-> Faraday's Constant [C/mol]				#
      !#								#
      !# Variable Definitions:						#
      !# =====================						#
      !# 								#
      !# TK	-> Cathode Temperature, currently set to constant will  #
      !#	   be changed when modified to incorporate transient    #
      !#	   analysis. [K]					#
      !# 								#
      !# TA	-> Anode Temperature, currently set to constant, will   #
      !#	   be changed when modified to incorpoarate transient   #
      !#	   analysis. [K]					#
      !# 								#
      !# TZ	-> Average cell temperature, mainly to be used as an 	#
      !#	   initial condition. [K]				#
      !#								#
      !# PK	-> Cathode pressure in the channel [MPa]		#
      !#								#
      !# PA	-> Anode pressure in the channel [MPa]			#
      !# 								#
      !# P0	-> Atmospheric pressure [MPa]				#
      !#								#
      !# DGS	-> Delta G*: Shift in Gibb's free energy for non-stand- #
      !#	   ard temperature and pressure. Calculated considering #
      !#	   the hydrogen is developed on the cathode and the ox- #
      !#	   ygen and water are present on the anode. [J/mol]	#
      !#								#
      !# DG	-> Delta G: Gibbs free energy at any temperature or pre-#
      !#	   ssure within a temperature range of 300-4000K [J/mol]#
      !#								#
      !# OCV	-> Open Circuit Voltage of PEM Electrolysis cell. [V]   #
      !#								#
      !# AJ(i)	-> Empirical coefficient for enthalpy and entropy calc- #
      !#	   ulation for (i) species. [non-dim]			#
      !#		i=1 -> Water					#
      !#		i=2 -> Hydrogen					#
      !#		i=3 -> Oxygen					#
      !#								#
      !# BJ(i)	-> Empirical coefficient for enthalpy and entropy calc- #
      !#	   ulation for (i) species. [non-dim]			#
      !#		i=1 -> Water					#
      !#		i=2 -> Hydrogen					#
      !#		i=3 -> Oxygen					#
      !# CJ(i)	-> Empirical coefficient for enthalpy and entropy calc- #
      !#	   ulation for (i) species. [non-dim]			#
      !#		i=1 -> Water					#
      !#		i=2 -> Hydrogen					#
      !#		i=3 -> Oxygen					#
      !#								#
      !# DJ(i)	-> Empirical coefficient for enthalpy and entropy calc- #
      !#	   ulation for (i) species. [non-dim]			#
      !#		i=1 -> Water					#
      !#		i=2 -> Hydrogen					#
      !#		i=3 -> Oxygen					#
      !#								#
      !# TG(i) 	-> Temperatures for calculation of open circuit voltage #
      !#	   to be used in the Gibbs free energy eqn. [K]		#
      !#								#
      !# H(i) 	-> Molar Enthalpy for species "i" [J/mol]		#
      !#		i=1 -> Water					#
      !#		i=2 -> Hydrogen					#
      !#		i=3 -> Oxygen					#
      !#								#
      !# S(i) 	-> Molar Entropy for species "i" [J/(mol-K)]		#
      !#		i=1 -> Water					#
      !#		i=2 -> Hydrogen					#
      !#		i=3 -> Oxygen					#
      !#								#
      !# C(i,j,k) -> Molar concentration of species "i" in layer "j" on #
      !#	   side "k." [mol/m^3]					#
      !#		i=1 -> Water					#
      !#		i=2 -> Hydordgen				#
      !#		i=3 -> Oxygen					#
      !#		---------------------------------------		#
      !#		j=1 -> Channel					#
      !#		j=2 -> Membrane (Place on Cathode side)		#
      !#		j=3 -> Electrode				#
      !#		j=4 -> Catalyst Layer				#
      !#		---------------------------------------		#
      !#		k=1 -> Anode Side				#
      !#		k=2 -> Cathode Side				#
      !#								#
      !# Y(i,j,k) -> Molar fraction of species "i" in layer "j" on side #
      !#	   "k." [mol/m^3]					#
      !#		i=1 -> Water					#
      !#		i=2 -> Hydordgen				#
      !#		i=3 -> Oxygen					#
      !#		---------------------------------------		#
      !#		j=1 -> Channel					#
      !#		j=2 -> Membrane (Place on Cathode side)		#
      !#		j=3 -> Electrode				#
      !#		j=4 -> Catalyst Layer				#
      !#		---------------------------------------		#
      !#		k=1 -> Anode Side				#
      !#		k=2 -> Cathode Side				#
      !#								#
      !# ETA(i,j) -> Parasitic losses from mass transport, electrode    #
      !#	   resistivity and activation energy [V]		#
      !#		i=1 -> Activation Voltage			#
      !#		i=2 -> Ohmic Losses				#
      !#		i=3 -> Mass Transport Losses			#
      !#		----------------------------------------------	#
      !#		j=1 -> Anode Side/Membrane for Ohmic Losses	#
      !#		j=2 -> Cathode Side/Electrode for Ohmic Losses	#
      !#								#
      !# AA	-> Effective transfer coefficient (alpha_anode)		#
      !#								#
      !# AK	-> Effective transfer coefficient (alpha_cathode)	#
      !#								#
      !# CD0(i)	-> Exchange current density [A/cm^2]	                #
      !#		i=1 -> Anode side       			#
      !#                i=2 -> Cathode side                             #
      !#								#
      !# LAM    -> Membrane hydration lambda in Mole H2O / Mole SO3     #
      !#                                                                #
      !# D(i,j) -> Cell dimmensions [m]                                 #
      !#                i=1 -> Channel width                            #
      !#                i=2 -> Channel height                           #
      !#                i=3 -> Land width                               #
      !#                i=4 -> Cell width                               #
      !#                i=5 -> Electrode height                         #
      !#                i=6 -> Plate height (BP thick - channel height) #
      !#                i=7 -> MEA Length (j=2, placed on "cathode"     #
      !#                i=8 -> Cell Area (j=2, placed on "cathode" side #
      !#                i=9 -> Membrane Thickness (j=2, on cathode side #
      !#                ----------------------------------------------- #
      !#                j=1 -> Anode side                               #
      !#                j=2 -> Cathode side                             #
      !#                                                                #
      !# R(i,j) -> Electrical resistance matrix [Ohm]                   #
      !#                i=1 -> Equivalent resistance                    #
      !#                i=2 -> Electrode resistance under land          #
      !#                i=3 -> Electrode resistance under channel       #
      !#                i=4 -> Electrode resistance in-plane under chan.#
      !#                i=5 -> Bipolar plate rib resistance             #
      !#                i=6 -> Bipolar plate backing resistance         #
      !#                ------------------------------------------------#
      !#                j=1 -> Anode side                               #
      !#                j=2 -> Cathode side                             #
      !#                ------------------------------------------------#
      !#                i=7, j=1 -> R1 Intermediate solution            #
      !#                i=7, j=2 -> RA Intermediate solution            #
      !#                i=8, j=1 -> RB Intermediate solution            #
      !#                i=8, j=2 -> RR Intermediate solution            #
      !#                i=9, j=1 -> RS Intermediate solution            #
      !#                i=9, j=2 -> RT Intermediate solution            #
      !#                i10, j=1 -> Req,ch anode                        #
      !#                i10, j=2 -> Req,ch cathode                      #
      !#                                                                #
      !# MR(i,j)-> Material resistivity [Ohm-M]                         #
      !#                i=1 -> Electrode                                #
      !#                i=2 -> Bipolar plate                            #
      !#                --------------------                            #
      !#                j=1 -> Anode                                    #
      !#                j=2 -> Cathode                                  #
      !#                                                                #
      !# N(i)   -> Number integer array                                 #
      !#                i=1 -> Number of channels anode                 #
      !#                i=2 -> Number of channels cathode               #
      !#                                                                #
      !# SIGM   -> Proton conductivity in the membrane                  #
      !#                                                                #
      !# CUR    -> Operating Current [A]                                #
      !#                                                                #
      !# CD     -> User defined current density [A/cm^2]                #
      !#                                                                #
      !# CSTP   -> Current step size determined by maximum current val- #
      !#           ue defined by the user in the input file.            #
      !#                                                                #
      !# RES    -> Number of steps in the current density to be solved  #
      !#           for through the simulation.                          #
      !#                                                                #
      !# ZCD    -> Current at a voltage step (K)                        #
      !#                                                                #
      !# ND     -> Electro-osmotic drag coefficient [molH2O/molH+]      #
      !#                                                                #
      !# MFR(i,j)> Molar flow rates [mol/s]                             #
      !#                i=1 -> Production/Consumption rate              #
      !#                i=2 -> Flow through membrane (currently calcul- #
      !#                       ated with constant densities not a func- #
      !#                       tion of temperature. Must add an empiri- #
      !#                       cal formulation for density as a func-   #
      !#                       tion of temperature.                     #
      !#                i=3 -> Diffusion through membrane               #
      !#                i=4 -> Pressure driven flow through membrane    #
      !#                i=5 -> Electro-osmotic drag through membrane    #
      !#                i=6 -> Flow out of anode                        #
      !#                i=7 -> Flow out of cathode                      #
      !#                --------------------------------------------    #
      !#                j=1 -> Water                                    #
      !#                j=2 -> Hydrogen                                 #
      !#                j=3 -> Oxygen                                   #
      !#                                                                #
      !# KD     -> Darcy's Law coefficient                              #
      !#                                                                #
      !# DW     -> Membrane diffusivity coefficient                     #
      !#                                                                #
      !# RHO(i) -> Density of water [kg/m^3]                            #
      !#                i=1 -> Anode side                               #
      !#                i=2 -> Cathode side                             #
      !#                                                                #
      !# MU     -> Viscosity of water                                   #
      !#                                                                #  
      !# EP(i)  -> Electrode Porosity                                   #
      !#                i=1 -> Anode side                               #
      !#                i=2 -> Cathode side                             #
      !#                                                                #
      !# EPP(i) -> Electrod Percolation threshold                       #
      !#                i=1 -> Anode side                               #
      !#                i=2 -> Cathode side                             #
      !#                                                                #
      !# M(i)   -> Molecular Weight [g/mol]                             #
      !#                i=1 -> Water                                    #
      !#                i=2 -> Hydrogen                                 #
      !#                i=3 -> Oxygen                                   #
      !#                                                                #
      !# DE(i,j)-> Diffusion coefficients [m^2/s]                       #
      !#                i=1 -> Effective diffusion coefficient          #
      !#                i=2 -> Non-porosity corrected diffusion coeff.  #
      !#                ----------------------------------------------  #
      !#                j=1 -> Anode side                               #
      !#                j=2 -> Cathode side                             #
      !#                                                                #
      !# TC(i)  -> Critical temperature of speicies:                    #
      !#                i=1 -> Water                                    #
      !#                i=2 -> Hydrogen                                 #
      !#                i=3 -> Oxygen                                   #
      !#                                                                #
      !# PC(i)  -> Critical pressure of species:                        #
      !#                i=1 -> Water                                    #
      !#                i=2 -> Hydrogen                                 #
      !#                i=3 -> Oxygen                                   #
      !#                                                                #
      !# FX(i,j)-> Molar Fluxes [mol/m^3]                               #
      !#                i=1 -> Water                                    #
      !#                i=2 -> Hydrogen                                 #
      !#                i=3 -> Oxygen                                   #
      !#                ---------------                                 #
      !#                j=1 -> Anode side                               #
      !#                j=2 -> Cathode side                             #
      !#                                                                #
      !#  DC(i) -> Binary diffusion constants                           #
      !#                i=1 -> Diffusion constant "a"                   #
      !#                i=2 -> Diffusion constant "b"                   #
      !#                                                                #   
      !##################################################################
      !
      !
      !
      PROGRAM ZERODELEKTROLYSE
      !
      IMPLICIT NONE
      !
      INTEGER :: J,K,RES,I
      INTEGER, ALLOCATABLE :: N(:)
      DOUBLE PRECISION :: TK,TA,TZ,PK,PA,P0,DGS,DG,OCV,AA,AK,CD0A, &
                          CD0K,LAM,SIGM,CD,CUR,CSTP,ZCD,ND,KD,DW, &
                          MU,E,HI,SI
      DOUBLE PRECISION, ALLOCATABLE :: AJ(:),BJ(:),CJ(:),DJ(:),TG(:), &
				       H(:),S(:),C(:,:,:),Y(:,:,:), &
				       ETA(:,:),D(:,:),R(:,:),MR(:,:), &
                                       CD0(:),MFR(:,:),M(:),DE(:,:), &
                                       DC(:),TC(:),PC(:),EP(:),EPP(:), &
                                       FX(:,:),RHO(:),P(:)
      DOUBLE PRECISION, PARAMETER :: PI=4.D0*DATAN(1.D0)
      DOUBLE PRECISION, PARAMETER :: RGC=8.3144621D0
      DOUBLE PRECISION, PARAMETER :: F=96485.D0
      CHARACTER (LEN=100) :: TITEL
      !
      ALLOCATE (Y(3,4,2))
      ALLOCATE (C(3,4,2))
      ALLOCATE (MFR(9,3))
      ALLOCATE (ETA(3,2))
      ALLOCATE (MR(2,2))
      ALLOCATE (R(10,2))
      ALLOCATE (DE(2,2))
      ALLOCATE (FX(3,2))
      ALLOCATE (EPP(2))
      ALLOCATE (RHO(2))
      ALLOCATE (D(9,2))
      ALLOCATE (CD0(2))
      ALLOCATE (EP(2))
      ALLOCATE (AJ(3))
      ALLOCATE (BJ(3))
      ALLOCATE (CJ(3))
      ALLOCATE (DJ(3))
      ALLOCATE (DC(2))
      ALLOCATE (TG(3))
      ALLOCATE (TC(3))
      ALLOCATE (PC(3))
      ALLOCATE (H(3))
      ALLOCATE (M(3))
      ALLOCATE (N(2))
      ALLOCATE (P(3))
      ALLOCATE (S(3))
      !
      CALL SYSTEM("rm -rf OUTPUT.dat")
      !
      !
      !
      OPEN(UNIT=10,FILE='PEM_einstellung.dat',ERR=100)
      !
      READ(10,*)
      READ(10,*)
      READ(10,*) TITEL
      !
      ! Solver Settings
      READ(10,*)
      READ(10,*)
      READ(10,*) RES
      !
      ! Initial Concentrations
      READ(10,*)
      READ(10,*)
      READ(10,*) Y(1,1,1),Y(2,1,2),Y(3,1,1)
      !
      ! Operating Conditions
      READ(10,*)
      READ(10,*)
      READ(10,*) CD
      !
      ! Initial Temperature Settings
      READ(10,*)
      READ(10,*)
      READ(10,*) TK,TA,TZ
      !
      ! Initial Pressure Settings
      READ(10,*)
      READ(10,*)
      READ(10,*) PK,PA,P0
      !
      ! OCV reaction coefficients and exchange current densities
      READ(10,*)
      READ(10,*)
      READ(10,*) AA,AK,CD0(1),CD0(2)
      !
      ! Membrane Parameters
      READ(10,*)
      READ(10,*)
      READ(10,*) LAM
      !
      ! Number of Channels (anode,cathode)
      READ(10,*)
      READ(10,*)
      READ(10,*) N(1),N(2)
      !
      !Cell/Stack Dimensions Anode Side
      READ(10,*)
      READ(10,*)
      READ(10,*) D(1,1),D(2,1),D(3,1),D(4,1),D(5,1),D(6,1)
      !
      ! Cell/Stack Dimension Cathode Side
      READ(10,*)
      READ(10,*)
      READ(10,*) D(1,2),D(2,2),D(3,2),D(4,2),D(5,2),D(6,2),D(7,2), &
                 D(8,2),D(9,2)
      !
      ! Electrical Resistivities
      READ(10,*)
      READ(10,*)
      READ(10,*) MR(1,1),MR(2,1),MR(1,2),MR(2,2)
      !
      ! Mass Transport Coefficients
      READ(10,*)
      READ(10,*)
      READ(10,*) ND,KD,DW,EP(1),EP(2),EPP(1),EPP(2)
      !
      READ(10,*)
      READ(10,*)
      READ(10,*) MU,M(1),M(2),M(3),DC(1),DC(2)
      !
      READ(10,*)
      READ(10,*)
      READ(10,*) TC(1),TC(2),TC(3),PC(1),PC(2),PC(3)
      CLOSE(10)
      !
      ! 
      !
      !#################################################################
      !
      ! Wasser Koeffizients
      !
      AJ(1)=180.D0
      BJ(1)=-85.4D0
      CJ(1)=15.6D0
      DJ(1)=-0.858D0
      TG(1)=TA
      P(1)=PA*Y(1,1,1)
      !
      ! Wasserstoff Koeffizients
      !
      AJ(2)=79.5D0
      BJ(2)=-26.3D0
      CJ(2)=4.23D0
      DJ(2)=-0.197D0
      TG(2)=TK
      P(2)=PK*Y(2,1,2)
      !
      ! Sauerstoff Koeffizients
      !
      AJ(3)=10.3D0
      BJ(3)=5.4D0
      CJ(3)=-0.18D0
      DJ(3)=0.D0
      TG(3)=TA
      P(3)=PA*Y(3,1,1)
      !
      !##################################################################
      !
      ! Operating Current
      CUR=CD*10.D1*D(8,2)
      !
      CSTP=CUR/RES
      PRINT*, CUR, CSTP
      !
      !##################################################################
      !#                                                                #
      !#                        Activation Losses                       #
      !#                                                                #
      !##################################################################
      !
      DO I=1,3
      ! ENTHALPY
	CALL ENTHALPY(TG(I),AJ(I),BJ(I),CJ(I),DJ(I),HI)
        H(I)=HI
      ! ENTROPY
	CALL ENTROPY(TG(I),P0,AJ(I),BJ(I),CJ(I),DJ(I),RGC,SI)
	S(I)=SI
      !
      PRINT*,H(I),S(I)
      END DO
      !
      ! Gibb's free energy
      DGS=(H(2)+0.5D0*H(3)-H(1))-TZ*(S(2)+0.5D0*S(3)-S(1))
      DG=DGS+RGC*TZ*DLOG((P(2)*P(3)**(0.5D0))/P(1))
      !
      ! Open Circuit Voltage
      OCV=-DG/(2*F)
      !
      !
      !##################################################################
      !#                                                                #
      !#                        Ohmic Losses                            #
      !#                                                                #
      !##################################################################
      !
      ! Electrical resistance of electrode under land anode
      R(2,1)=MR(1,1)*(D(5,1)/((D(3,1)/2.D0)*D(7,2)))
      !
      ! Electrical resistance of electrode under channel anode
      R(3,1)=MR(1,1)*(D(5,1)/(D(1,1)*D(7,2)))
      !
      ! Electrical resistance in-plane under channel anode
      R(4,1)=MR(1,1)*((D(1,1)/4.D0)/(D(5,1)*D(7,2)))
      !
      ! Electrical resistance Bipolar plate rib anode
      R(5,1)=MR(2,1)*(D(2,1)/((D(3,1)/2.D0)*D(7,2)))
      !
      ! Electrical resistance Bipolar plate backing anode
      R(6,1)=MR(2,1)*(D(6,1)/D(8,2))
      !
      ! Electrical resistance of electrode under land cathode
      R(2,2)=MR(1,2)*(D(5,2)/((D(3,2)/2.D0)*D(7,2)))
      !
      ! Electrical resistance of electrode under channel cathode
      R(3,2)=MR(1,2)*(D(5,2)/(D(1,2)*D(7,2)))
      !
      ! Electrical resistance in-plane under channel cathode
      R(4,2)=MR(1,2)*((D(1,2)/4.D0)/(D(5,2)*D(7,2)))
      !
      ! Electrical resistance Bipolar plate rib cathode
      R(5,2)=MR(2,2)*(D(2,2)/((D(3,2)/2.D0)*D(7,2)))
      !
      ! Electrical resistance Bipolar plate backing cathode
      R(6,2)=MR(2,2)*(D(6,2)/D(8,2))
      !
      !Anode Intermediate Steps
      R(7,1)=(2.D0*R(3,1)*R(4,1)+R(4,1)**2)/(R(3,1))
      R(7,2)=2.D0*R(3,1)+R(4,1)
      R(8,1)=(R(2,1)*R(7,2))/(R(2,1)+R(7,2))
      R(8,2)=(R(8,1)*R(5,1)+R(8,1)*R(7,1)+R(5,1)*R(7,1))/R(8,1)
      R(9,1)=(R(8,1)*R(5,1)+R(8,1)*R(7,1)+R(5,1)*R(7,1))/R(5,1)
      R(9,2)=(R(8,1)*R(5,1)+R(8,1)*R(7,1)+R(5,1)*R(7,1))/R(7,1)
      R(10,1)=(R(9,2)*(((R(8,1)*R(9,1))/(R(8,1)+R(9,1)))+((R(8,2)* &
              R(5,1))/(R(8,2)+R(5,1)))))/(R(9,2)+(((R(8,1)*R(9,1))/ &
              (R(8,1)+R(9,1)))+((R(8,2)*R(5,1))/(R(8,2)+R(5,1)))))
      !
      ! Equivalent Resistance Anode Side
      R(1,1)=(((R(2,1)+R(5,1))/2.D0)*((R(10,1)/N(1))+R(6,1)))/ &
             (((R(2,1)+R(5,1))/2.D0)+((R(10,1)/N(1))+R(6,1)))
      !
      !Anode Intermediate Steps
      R(7,1)=(2.D0*R(3,2)*R(4,2)+R(4,2)**2)/(R(3,2))
      R(7,2)=2.D0*R(3,2)+R(4,2)
      R(8,1)=(R(2,2)*R(7,2))/(R(2,2)+R(7,2))
      R(8,2)=(R(8,1)*R(5,2)+R(8,1)*R(7,1)+R(5,2)*R(7,1))/R(8,1)
      R(9,1)=(R(8,1)*R(5,2)+R(8,1)*R(7,1)+R(5,2)*R(7,1))/R(5,2)
      R(9,2)=(R(8,1)*R(5,2)+R(8,1)*R(7,1)+R(5,2)*R(7,1))/R(7,1)
      R(10,2)=(R(9,2)*(((R(8,1)*R(9,1))/(R(8,1)+R(9,1)))+((R(8,2)* &
              R(5,2))/(R(8,2)+R(5,2)))))/(R(9,2)+(((R(8,1)*R(9,1))/ &
              (R(8,1)+R(9,1)))+((R(8,2)*R(5,2))/(R(8,2)+R(5,2)))))
      !
      ! Equivalent Resistance Cathode Side
      R(1,2)=(((R(2,2)+R(5,2))/2.D0)*((R(10,2)/N(2))+R(6,2)))/ &
             (((R(2,2)+R(5,2))/2.D0)+((R(10,2)/N(2))+R(6,2)))
      !
      ! Membrane proton transport resistance
      SIGM=(0.005139D0*LAM-0.00326D0)*DEXP(1268.D0*((1.D0/303.D0)- &
           (1.D0/TZ)))
      !
      !##################################################################
      !#                                                                #
      !#                   Mass Transport Losses                        #
      !#                                                                #
      !##################################################################
      !
      ! Density of water anode side
      RHO(1)=1.D-05*TA**3-0.0139*TA**2+5.1319*TA+414.32

      !Density Subroutine (In development)
      CALL RHOW(TA,RHO(1))
      !
      ! Density of water cathode side
      RHO(2)=1.D-05*TK**3-0.0139*TK**2+5.1319*TK+414.32
      !
      ! Hydrostatic pressure flows
      MFR(4,1)=KD*((D(8,2)*RHO(1))/(MU*M(1)))
      !
      
      !##################################################################
      !#                                                                # 
      !#                 Diffusion Coefficients                         #
      !#                                                                #
      !##################################################################
      !
      ! Anode Diffusion Coefficient
      DE(2,1)=(DC(1)*(TA/(TC(3)*TC(1))**(0.5D0))**DC(2)*(PC(3)* &
              PC(1))**(1.D0/3.D0)*(TC(3)*TC(1))**(5.D0/12.D0)* &
              ((1.D0/M(3))+(1.D0/M(1)))**(0.5D0))/PA
      !
      !Effective Anode Diffusivity
      DE(1,1)=DE(2,1)*EP(1)*((EP(1)-EPP(1))/(1.D0-EPP(1)))**0.785D0
      !
      ! Cathode Diffusion Coefficient
      DE(2,2)=(DC(1)*(TK/(TC(2)*TC(1))**(0.5D0))**DC(2)*(PC(2)* &
              PC(1))**(1.D0/3.D0)*(TC(2)*TC(1))**(5.D0/12.D0)* &
              ((1.D0/M(2))+(1.D0/M(1)))**(0.5D0))/PK
      !
      !Effective Cathode Diffusivity
      DE(1,2)=DE(2,2)*EP(2)*((EP(2)-EPP(2))/(1.D0-EPP(2)))**0.785D0
      !
      !
      !
      !##################################################################
      !#                                                                #
      !#                        Current Loop                            #
      !#                                                                #
      !##################################################################
      !
      DO K=1,RES
      !
      ! Current at each time step
      ZCD=CSTP*DBLE(K)
      !
      ! Molar flow rates
      !
      ! Water consumption
      MFR(1,1)=ZCD/(2*F)
      !
      ! Hydrogen Production
      MFR(1,2)=ZCD/(2*F)
      !
      ! Oxygen Production
      MFR(1,3)=ZCD/(4*F)
      !
      ! Electro-osmotic drag
      MFR(5,1)=ND*(ZCD/F)
      ! 
      ! Water flow through the membrane (MUST EDIT WITH TEMP DEP. DENSITY)
      MFR(2,1)=(MFR(5,1)-MFR(4,1)+((D(8,2)*DW)/D(9,2))*(((RHO(2)- &
               RHO(1))/(M(1)))+((D(5,1)*MFR(1,1))/(DE(1,1)*D(8,2))))) &
               /(1.D0-(DW/D(9,2))*((D(5,2)/DE(1,2))+(D(5,2)/ &
                DE(1,1))))
      !
      !
      !
      !##################################################################
      !#                                                                #
      !#                           Flux Terms                           #
      !#                                                                #
      !##################################################################
      !
      ! Anodic Water Flux
      FX(1,1)=(MFR(2,1)+MFR(1,1))/D(8,2)
      !
      ! Cathodic Water Flux
      FX(1,2)=MFR(2,1)/D(8,2)
      !
      ! Oxygen Flux (Only Anode side)
      FX(3,1)=ZCD/(4.D0*F*D(8,2))
      !
      ! Hydrogen Flux (Only Cathode side)
      FX(2,2)=ZCD/(2.D0*F*D(8,2))
      !
      !
      !
      !##################################################################
      !#                                                                #
      !#                        Mole Fractions/                         #
      !#                        Concentrations                          #
      !#                                                                #
      !##################################################################
      !
      ! Mole fractions in the channels
      Y(2,1,2)=FX(2,2)/(FX(2,2)+FX(1,2))
      Y(1,1,2)=FX(1,2)/(FX(2,2)+FX(1,2))
      Y(3,1,1)=FX(3,1)/(FX(3,1)+FX(1,1))
      Y(1,1,1)=FX(1,1)/(FX(2,2)+FX(1,1))
      !
      ! Species concentrations in channel
      C(2,1,2)=(PK*Y(2,1,2))/(RGC*TK)
      C(1,1,2)=RHO(2)/M(1)
      C(3,1,1)=(PA*Y(3,1,1))/(RGC*TA)
      C(1,1,1)=RHO(1)/M(1)
      !
      ! Species concentrations at the membrane electrode interface 
      !(catalyst layer)
      C(2,4,2)=C(2,1,2)+((D(5,2)*FX(2,2))/DE(1,2))
      C(1,4,2)=C(1,1,2)+((D(5,2)*FX(1,2))/DE(1,2))
      C(3,4,1)=C(3,1,1)+((D(5,1)*FX(3,1))/DE(1,1))
      C(1,4,1)=C(1,1,1)+((D(5,1)*FX(1,1))/DE(1,1))
      !
      ! Mole fraction at the membrane electrode interfaces
      Y(2,4,2)=(RGC*TK*C(2,4,2))/PK
      Y(1,4,2)=(C(1,4,2)/(C(1,4,2)+C(2,4,2)))
      Y(3,4,1)=(RGC*TA*C(3,4,1))/PA
      Y(1,4,1)=(C(1,4,1)/(C(1,4,1)+C(3,4,1)))
      !
      !Partial Pressures (Updated)
      P(1)=Y(1,4,1)*PA
      P(2)=Y(2,4,2)*PK
      P(3)=Y(3,4,1)*PA
      !
      ! Mass transport losses
      ETA(3,1)=((RGC*TA)/(4.D0*F))*DLOG(C(3,4,1)/C(3,1,1))
      ETA(3,2)=((RGC*TK)/(2.D0*F))*DLOG(C(2,4,2)/C(2,1,2))
      !
      ! Activation Losses
      ETA(1,1)=((RGC*TA)/(AA*F))*DASINH((ZCD/D(8,2))/(2*CD0(1)))
      ETA(1,2)=((RGC*TK)/(AK*F))*DASINH((ZCD/D(8,2))/(2*CD0(2)))
      !
      ! Ohmic losses through the membrane
      ETA(2,1)=D(9,2)*ZCD/(D(8,2)*SIGM)
      !
      ! Ohmic losses through the electrode
      ETA(2,2)=(R(1,1)+R(1,2))*ZCD
      !
      !E=OCV+ETA(1,1)+ETA(1,2)+ETA(2,2)!+ETA(3,1)!+ETA(3,2)!+ETA(2,1)
      E=OCV+ETA(1,1)+ETA(1,2)+ETA(2,1)+ETA(2,2)+ETA(3,2)+ETA(3,1)
      OPEN(UNIT=20,FILE="OUTPUT.dat",ACCESS='APPEND')
       WRITE(20,*) ZCD, E
      CLOSE(20)
      END DO
      PRINT*, ZCD, E
      PRINT*, Y(1,1,1),Y(2,1,2),Y(3,1,1)
      !
      !
      !
      !FEHLER ANGABEN
      !
  100 PRINT*, 'Laden Fehler mit PEM_einstellung.dat' 
      ! 
      END PROGRAM ZERODELEKTROLYSE
      !
      !
      !
      SUBROUTINE RHOW(T,RHO)

      IMPLICIT NONE
      DOUBLE PRECISION :: T,RHO
      RHO=1.D-05*T**3-0.0139*T**2+5.1319*T+414.32
      RETURN
      END SUBROUTINE RHOW


      !##################################################################
      !#								#
      !#		SUBROUTINE ENTHALPY 				#
      !#								#
      !#================================================================#
      !# Empirical formulation taken from Stull DR, Porphet H. JANAF    #
      !# Thermochemical tables. NSRDSNBS June 1971;37			#
      !#								#
      !##################################################################
      !
      SUBROUTINE ENTHALPY(TG,AJ,BJ,CJ,DJ,HI)
      !
      IMPLICIT NONE
      !
      DOUBLE PRECISION :: TG,AJ,BJ,CJ,DJ,HI
      !
        HI=AJ*TG+0.8D0*(BJ)*TG**1.25D0+(2.D0/3.D0)*(CJ)*TG**1.5D0+ &
           (4.D0/7.D0)*(DJ)*TG**(7.D0/4.D0)
      !
      RETURN
      !
      END SUBROUTINE ENTHALPY
      !
      !
      !
      !##################################################################
      !#								#
      !#		SUBROUTINE ENTROPY 				#
      !#								#
      !#================================================================#
      !# Empirical formulation taken from Stull DR, Porphet H. JANAF    #
      !# Thermochemical tables. NSRDSNBS June 1971;37			#
      !#								#
      !##################################################################
      !
      SUBROUTINE ENTROPY(TG,P0,AJ,BJ,CJ,DJ,GASC,SI)
      !
      IMPLICIT NONE
      !
      DOUBLE PRECISION :: TG,P0,AJ,BJ,CJ,DJ,GASC,SI

        SI=AJ*DLOG(TG)+4.D0*BJ*TG**(0.25D0)+2.D0*CJ*TG**(0.5D0)+(4.D0/ &
           3.D0)*DJ*TG**(0.75D0)-GASC*DLOG(P0)
      !
      RETURN
      !
      END SUBROUTINE ENTROPY

