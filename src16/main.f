***********************************************************************
      PROGRAM EQ_CHEMISTRY
***********************************************************************
      use PARAMETERS,ONLY: model_dim,model_struc,model_eqcond,
     >                     useDatabase,auto_atmos,adapt_cond
      use EXCHANGE,ONLY: chemcall,chemiter,ieqcond,ieqconditer,
     >                   itransform,preEst,preUse,preIter,
     >                   DUALcorr,HCOcorr
      use DATABASE,ONLY: NLAST
      implicit none

      call READ_PARAMETER
      call INIT
      call INIT_CHEMISTRY
      call INIT_DUSTCHEM
      
      if (model_dim==0) then
        if (adapt_cond) then
          call ADAPT_CONDENSATES
        else
          call DEMO_CHEMISTRY
        endif  
      else if (model_dim==1) then  
        if (auto_atmos) then
          call AUTO_STRUCTURE
        else if (model_struc==0) then 
          call DEMO_SWEEP
        else  
          call DEMO_STRUCTURE
        endif 
      else if (model_dim==2) then  
        call DEMO_PHASEDIAGRAM
      else
c        print*,'*** model_dim=',model_dim,' ???'
        stop
      endif   
      
c      print'("            smchem calls = ",I8)',chemcall
c      print'("         iterations/call = ",0pF8.2)',
c     >                     REAL(chemiter)/REAL(chemcall)
c      print'("     pre-iterations/call = ",0pF9.3)',
c     >                      REAL(preIter)/REAL(chemcall)
c      print'("usage of saved estimates = ",0pF9.3," %")',
c     >                       REAL(preUse)/REAL(preEst)*100.0
c      print'("   dual corrections/call = ",0pF9.3)',
c     >                     REAL(DUALcorr)/REAL(chemcall)
c      print'("  H-C-O corrections/call = ",0pF9.3)',
c     >                      REAL(HCOcorr)/REAL(chemcall)
      if (model_eqcond) then
c        print'("   eq condensation calls = ",I8)',ieqcond
c        print'("      eq iterations/call = ",0pF8.2)',
c     >                   REAL(ieqconditer)/REAL(ieqcond)
c        print'("         transform calls = ",I8)',itransform
        NLAST=0         ! also save replaced database entries
        if (useDatabase) call SAVE_DBASE
      endif

      end


***********************************************************************
      SUBROUTINE DEMO_CHEMISTRY
***********************************************************************
      use PARAMETERS,ONLY: nHmax,Tmax,pmax,model_pconst,model_eqcond,
     >                     verbose,Nseq,Tseq 
      use CHEMISTRY,ONLY: NMOLE,NELM,m_kind,m_anz,elnum,cmol,el
      use DUST_DATA,ONLY: NELEM,NDUST,elnam,eps0,bk,bar,amu,
     >                    muH,mass,mel,
     >                    dust_nam,dust_mass,dust_Vol,
     >                    dust_nel,dust_el,dust_nu
      use EXCHANGE,ONLY: nel,nat,nion,nmol,mmol,H,C,N,O,Si
      implicit none
      integer,parameter  :: qp = selected_real_kind ( 33, 4931 )
      real(kind=qp) :: eps(NELEM),Sat(NDUST),eldust(NDUST)
      real(kind=qp) :: nmax,threshold,deps
      real*8,parameter :: pi=3.14159265358979D+0
      real*8  :: Tg,nHges,p,mu,muold,pgas,fold,ff,dfdmu,dmu,mugas,Vol
      real*8  :: rhog,rhoc,rhod,d2g,mcond,mgas,Vcon,Vcond,Vgas,ngas
      real*8  :: nkey,nkeyt,tchem,tchemtot,AoverV,mic,atyp,alpha,vth
      real*8  :: yr,Myr,molm,molmass,stoich,HH,OO,CC,NN
      integer :: i,imol,iraus,e,aIraus,aIIraus,j,verb,dk,it,stindex
      integer :: k,keyel,imax,dustst
      integer :: H2O,CO2,CH4,O2,H2,N2,NH3,CO,OCS,SO2,S2,H2S
      logical :: included,haeufig,raus(NMOLE)
      logical :: rausI(NELEM),rausII(NELEM)
      character(len=10) :: sp
      character(len=20) :: limcond
      integer :: iseq

c     variable declaration for writing of partial pressures
      real*8 :: n_atom_total, n_mole_total, n_total, pp_electron,
     >  pp_atoms, pp_molecules, pp_non, pp_sum
      
      nHges = nHmax
      Tseq(Nseq) = Tmax
      p     = pmax
      eps   = eps0
      mu    = muH
      eldust = 0.Q0
      verb = verbose
      if (model_eqcond) verb=0

      do iseq=1,Nseq         ! possibility to set a sequence of T-points hot->cold
        Tg = Tseq(iseq)      ! (not used when default Nseq=1)
        do it=1,999
          if (model_pconst) nHges = p*mu/(bk*Tg)/muH
          if (model_eqcond) then
            call EQUIL_COND(nHges,Tg,eps,Sat,eldust,verbose)
          endif
          call GGCHEM(nHges,Tg,eps,.false.,verb)
          ngas = nel
          do j=1,NELEM
            ngas = ngas + nat(j)
          enddo
          do j=1,NMOLE
            ngas = ngas + nmol(j)
          enddo
          pgas  = ngas*bk*Tg
          ff    = p-pgas
          if (it==1) then
            muold = mu
            mu = nHges/pgas*(bk*Tg)*muH
            dmu = mu-muold
            if (.not.model_pconst) exit
          else
            dfdmu = (ff-fold)/(mu-muold)
            dmu   = -ff/dfdmu
            !write(98,'(I3,99(1pE14.7))') it,muold,mu,fold,ff,dfdmu,dmu/mu
            muold = mu
            if ((dmu>0.0).or.ABS(dmu/mu)<0.7) then
              mu = muold+dmu
            else
              mu = nHges/pgas*(bk*Tg)*muH
            endif  
          endif
          fold = ff
c          print '("p-it=",i3,"  mu=",2(1pE20.12))',it,mu/amu,dmu/mu
          if (ABS(dmu/mu)<1.E-10) exit
        enddo  
      enddo

      print*
c      write(*,*) '----- total nuclei dens. and fractions in gas -----'
      do e=1,NELM
        if (e==el) cycle
        i = elnum(e) 
c        write(*,'(" n<",A2,">=",1pE10.4,2x,1pE10.4)')
c     >      elnam(i),nHges*eps(i),eps(i)/eps0(i)
      enddo  

      ngas = nel      ! cm-3
      rhog = nel*mel  ! g cm-3
      do j=1,NELEM
        ngas = ngas + nat(j)
        rhog = rhog + nat(j)*mass(j)
      enddo
      do j=1,NMOLE
        ngas = ngas + nmol(j)
        rhog = rhog + nmol(j)*mmol(j)
      enddo
      rhod = 0.d0     ! g cm-3  
      Vcon = 0.d0     ! cm3
      do i=1,NDUST
        if (eldust(i)<=0.Q0) cycle 
        rhod = rhod + nHges*eldust(i)*dust_mass(i)
        Vcon = Vcon + nHges*eldust(i)*dust_Vol(i)
      enddo  
      mugas = rhog/ngas             ! mean molecular weight [g]
      d2g   = rhod/rhog             ! dust/gas mass ratio [-]
      if (Vcon>0.d0) then
      rhoc  = rhod/Vcon             ! condensed matter density [g cm-3]
      mcond = 1.0                   ! mass of condensates [= 1 g]
      Vcond = mcond/rhoc            ! volume of condensates [cm3]
      mgas  = mcond/d2g             ! mass of gas [g]
      Vgas  = mgas/mugas*bk*Tg/pgas ! volume of gas [cm3]

c      print*
c      write(*,*) '----- bulk properties -----'
c      print'("for Tg[K]=",0pF8.2," and n<H>[cm-3]=",1pE10.3)',Tg,nHges
c      print'(16x,3(A12))',"condensates","gas","check"
c      print'("         mass[g]",2(1pE12.4))',mcond,mgas
c      print'("  density[g/cm3]",2(1pE12.4))',rhoc,rhog
c      print'("tot.dens.[g/cm3]",12x,2(1pE12.4))',nHges*muH,
c     >                                   (mcond+mgas)/Vgas
c      print'("     volume[cm3]",2(1pE12.4))',Vcond,Vgas
c      print'("   pressure[bar]",12x,1pE12.4)',pgas/bar
c      print'("  el.press.[bar]",12x,1pE12.4)',nel*bk*Tg/bar
c      print'(" mol.weight[amu]",12x,1pE12.4)', mugas/amu
c      print'("         mu[amu]",12x,2(1pE12.4))', mu/amu,
c     >                (mcond+mgas)/Vgas/pgas*(bk*Tg)/amu
c      print'("        muH[amu]",12x,2(1pE12.4))',muH/amu,
c     >                     (mcond+mgas)/(nHges*Vgas)/amu
      
c      print*
c      write(*,*) '----- condensates [cm3] [mfrac] [Vfrac] -----'
      raus = .false.
      do 
        iraus = 0
        nmax  = 0.Q0
        do i=1,NDUST
          if (raus(i).or.eldust(i)<=0.Q0) cycle 
          if (eldust(i)>nmax) then
            iraus = i
            nmax = eldust(i)
          endif
        enddo
        if (iraus==0) exit
        raus(iraus) = .true.
c        write(*,1020) ' n'//trim(dust_nam(iraus))//'=',
c     >                eldust(iraus)*nHges,
c     >                eldust(iraus)*dust_mass(iraus)*nHges/rhod,
c     >                eldust(iraus)*dust_Vol(iraus)*nHges/Vcon
      enddo
      endif
      
c      write(*,*) '----- atoms and ions [cm3] -----'
c      write(*,1000) ' nel=',nel
      do e=1,NELM
        if (e==el) cycle
        i = elnum(e) 
c        write(*,1010) ' n'//trim(elnam(i))//'I=',nat(i),
c     >               '  n'//trim(elnam(i))//'II=',nion(i)
      enddo
  
c      write(*,*) '----- most abundant species [cm3] [molfrac] -----'
      raus   = .false.
      rausI  = .false.
      rausII = .false.
      do
        iraus   = 0
        aIraus  = 0
        aIIraus = 0
        nmax  = 0.Q0
        do i=1,NMOLE
          if ((nmol(i).gt.nmax).and.(.not.raus(i))) then
            iraus = i
            nmax  = nmol(i)
          endif
        enddo
        do e=1,NELM
          if (e==el) cycle
          i = elnum(e) 
          if ((nat(i).gt.nmax).and.(.not.rausI(i))) then
            iraus = 0
            aIraus = i
            aIIraus = 0
            nmax  = nat(i)
          endif
          !if ((nion(i).gt.nmax).and.(.not.rausII(i))) then
          !  iraus = 0
          !  aIraus = 0
          !  aIIraus = i
          !  nmax  = nion(i)
          !endif
        enddo
        haeufig = (nmax.gt.ngas*1.Q-9)
        if (.not.haeufig) exit
        if (iraus>0) then
          raus(iraus) = .true.
c          write(*,4010) cmol(iraus),nmol(iraus),
c     >                  nmol(iraus)/ngas,nmol(iraus)/ngas
        else if (aIraus>0) then 
          rausI(aIraus) = .true.
c          write(*,4010) elnam(aIraus)//"I       ",nat(aIraus),
c     >                  nat(aIraus)/ngas,nat(aIraus)/ngas
        !else if (aIIraus>0) then 
        !  rausII(aIIraus) = .true.
        !  write(*,4010) elnam(aIIraus)//"II     ",nion(aIIraus),
     >  !                nion(aIIraus)/ngas,nion(aIIraus)/ngas
        endif  
      enddo
      iraus = stindex(cmol,NMOLE,'H2')
c      if (.not.raus(iraus))
c     >   write(*,4010) cmol(iraus),nmol(iraus),
c     >                 nmol(iraus)/ngas,nmol(iraus)/ngas
      iraus = stindex(cmol,NMOLE,'O2')
c      if (.not.raus(iraus))
c     >   write(*,4010) cmol(iraus),nmol(iraus),
c     >                 nmol(iraus)/ngas,nmol(iraus)/ngas
      iraus = stindex(cmol,NMOLE,'CH4')
c      if (.not.raus(iraus))
c     >   write(*,4010) cmol(iraus),nmol(iraus),
c     >                 nmol(iraus)/ngas,nmol(iraus)/ngas
      iraus = stindex(cmol,NMOLE,'CO2')
c      if (.not.raus(iraus))
c     >   write(*,4010) cmol(iraus),nmol(iraus),
c     >                 nmol(iraus)/ngas,nmol(iraus)/ngas
      iraus = stindex(cmol,NMOLE,'NH3')
c      if (.not.raus(iraus))
c     >   write(*,4010) cmol(iraus),nmol(iraus),
c     >                 nmol(iraus)/ngas,nmol(iraus)/ngas
          
c      write(*,*) '-----  where are the elements?  -----'
      do e=1,NELM
        i = elnum(e)
        if (e==el) then
c          write(*,'("    Element ",A2,1pE15.3)') 'el',0.Q0
c          write(*,'(1x,A18,1pE10.3)') "nel",nel
          threshold = 1.Q-3*nel
        else   
c          write(*,'("    Element ",A2,1pE15.3)') elnam(i),eps0(i)*nHges 
          threshold = eps(i)*nHges*1.D-5
          if (nat(i).gt.threshold) then
c           write(*,'(1x,A18,1pE10.3)') "n"//trim(elnam(i)), nat(i) 
          endif  
        endif  

        raus = .false.
        do 
          iraus = 0
          nmax  = 0.Q0
          do dk=1,NDUST
            if (eldust(dk)<=0.Q0) cycle 
            included = .false. 
            do j=1,dust_nel(dk)
              if (i==dust_el(dk,j)) then
                included = .true.
              endif
            enddo  
            if (included) then
              if ((eldust(dk).gt.nmax).and.(.not.raus(dk))) then
                iraus = dk
                nmax = eldust(dk)
              endif  
            endif
          enddo
          if (iraus==0) exit
          haeufig = (nmax.gt.eps0(i)*1.D-2)
          if (.not.haeufig) exit
c          write(*,'(1x,A18,1pE10.3)') 
c     >          "n"//trim(dust_nam(iraus)),eldust(iraus)*nHges 
          raus(iraus) = .true.
        enddo  

        raus = .false.
        do 
          iraus = 0
          nmax  = 0.Q0
          do imol=1,NMOLE
            sp = cmol(imol) 
            if ((nmol(imol).gt.nmax).and.(.not.raus(imol))) then
              included = .false. 
              do j=1,m_kind(0,imol)
                if (e==m_kind(j,imol)) included=.true.
              enddo  
              if (included) then
                iraus = imol
                nmax = nmol(imol)
              endif  
            endif
          enddo  
          haeufig = (nmax.gt.threshold)
          if (.not.haeufig) exit
c          write(*,'(1x,A18,1pE10.3)') "n"//trim(cmol(iraus)),nmol(iraus)
          raus(iraus) = .true.
        enddo
      enddo  

*     ------------------------------
      call SUPERSAT(Tg,nat,nmol,Sat)
*     ------------------------------
c      print*
c      write(*,*) '----- supersaturation ratios -----'
      do i=1,NDUST
        if (Sat(i)<1.Q-2) cycle 
c        write(*,5000) dust_nam(i),Sat(i) 
      enddo

*     ----------------------------------------------------
*     ***  atmosphere types and low-temp expectations  ***
*     ----------------------------------------------------
      if (.true.) then
c        write(*,'(" epsH=",2(1pE12.4))') eps0(H),eps(H)
c       write(*,'(" epsC=",2(1pE12.4),"   C/(H+O+C)=",1pE12.4)') 
c     >         eps0(C),eps(C),eps(C)/(eps(H)+eps(O)+eps(C))
c        write(*,'(" epsO=",2(1pE12.4)," (O-H)/(O+H)=",1pE12.4)') 
c     >     eps0(O),eps(O),(eps(O)-eps(H))/(eps(O)+eps(H))
c        write(*,'(" epsN=",2(1pE12.4))') eps0(N),eps(N)
c        print*
        HH = eps(H)
        OO = eps(O)
        CC = eps(C)
        NN = eps(N)
        H2O = stindex(cmol,NMOLE,'H2O')
        CO2 = stindex(cmol,NMOLE,'CO2')
        CH4 = stindex(cmol,NMOLE,'CH4')
        H2  = stindex(cmol,NMOLE,'H2')
        N2  = stindex(cmol,NMOLE,'N2')
        O2  = stindex(cmol,NMOLE,'O2')
        NH3 = stindex(cmol,NMOLE,'NH3')
        CO  = stindex(cmol,NMOLE,'CO')
        OCS = stindex(cmol,NMOLE,'COS')
        SO2 = stindex(cmol,NMOLE,'SO2')
        S2  = stindex(cmol,NMOLE,'S2')
        H2S = stindex(cmol,NMOLE,'H2S')
        if (HH>2*OO+4*CC) then
          if (3*NN<HH-2*OO-4*CC) then
c            print'("type A1")'
c            print'("H2O:",2(0pF8.4))',nmol(H2O)/ngas,
c     >                                     (2*OO)/(HH-NN-2*CC)
c            print'("CH4:",2(0pF8.4))',nmol(CH4)/ngas,
c     >                                     (2*CC)/(HH-NN-2*CC)
c            print'("NH3:",2(0pF8.4))',nmol(NH3)/ngas,
c     >                                     (2*NN)/(HH-NN-2*CC)
c            print'(" H2:",2(0pF8.4))',nmol(H2)/ngas,
c     >                        (HH-2*OO-4*CC-3*NN)/(HH-NN-2*CC)
          else  
c            print'("type A2")'
c            print'("H2O:",2(0pF8.4))',nmol(H2O)/ngas,
c     >                                     (6*OO)/(HH+2*CC+3*NN+4*OO)
c            print'("CH4:",2(0pF8.4))',nmol(CH4)/ngas,
c     >                                     (6*CC)/(HH+2*CC+3*NN+4*OO)
c            print'("NH3:",2(0pF8.4))',nmol(NH3)/ngas,
c     >                           (2*HH-8*CC-4*OO)/(HH+2*CC+3*NN+4*OO)
c            print'(" N2:",2(0pF8.4))',nmol(N2)/ngas,
c     >                        (3*NN+4*CC+2*OO-HH)/(HH+2*CC+3*NN+4*OO)
          endif
        else if (OO>0.5*HH+2*CC) then
c          print'("type B")'
c          print'("H2O:",2(0pF8.4))',nmol(H2O)/ngas,
c     >                                     (2*HH)/(HH+2*OO+2*NN)
c          print'("CO2:",2(0pF8.4))',nmol(CO2)/ngas,
c     >                                     (4*CC)/(HH+2*OO+2*NN)
c          print'(" N2:",2(0pF8.4))',nmol(N2)/ngas,
c     >                                     (2*NN)/(HH+2*OO+2*NN)
c          print'(" O2:",2(0pF8.4))',nmol(O2)/ngas,
c     >                             (2*OO-HH-4*CC)/(HH+2*OO+2*NN)
        else if (CC>0.25*HH+0.5*OO) then
c          print*,"forbidden by graphite condensation"
        else 
c          print'("type C")'
c          print'("H2O:",2(0pF8.4))',nmol(H2O)/ngas,
c     >                             (HH+2*OO-4*CC)/(HH+2*OO+2*NN)
c          print'("CO2:",2(0pF8.4))',nmol(CO2)/ngas,
c    >                           (OO+2*CC-0.5*HH)/(HH+2*OO+2*NN)
c          print'("CH4:",2(0pF8.4))',nmol(CH4)/ngas,
c     >                           (0.5*HH-OO+2*CC)/(HH+2*OO+2*NN)
c          print'(" N2:",2(0pF8.4))',nmol(N2)/ngas,
c     >                                     (2*NN)/(HH+2*OO+2*NN)
        endif
      endif

      !print'(99(A9))',"SO2[ppm]","H2O[ppm]","OCS[ppm]","CO[ppm]",
     >!                "H2S[ppb]","S2[ppb]","H2[ppb]","O2[ppb]"
      !print'(99(0pF9.2))',nmol(SO2)/ngas/1.E-6,
     >!                    nmol(H2O)/ngas/1.E-6,
     >!                    nmol(OCS)/ngas/1.E-6,
     >!                    nmol(CO) /ngas/1.E-6,
     >!                    nmol(H2S)/ngas/1.E-9,
     >!                    nmol(S2) /ngas/1.E-9,
     >!                    nmol(H2) /ngas/1.E-9,
     >!                    nmol(O2) /ngas/1.E-9      
      !print'(99(A9))',"H2O[%]","CO2[%]","N2[%]","O2[%]"
      !print'(99(0pF9.4))',nmol(H2O)/ngas/1.E-2,
     >!                    nmol(CO2)/ngas/1.E-2,
     >!                    nmol(N2) /ngas/1.E-2,
     >!                    nmol(O2) /ngas/1.E-2
      
*     -----------------------------------------------------
*     ***  Calculation of the condenstation timescales  ***
*     -----------------------------------------------------
      if (.false.) then
c      print'("----- condensation timescales -----")'
      yr  = 365.25*24.0*3600.0
      Myr = 1.E+6*yr
      mic = 1.E-4                ! one micrometer [cm]
      atyp = 1.0*mic             ! typical particle size
      AoverV = 3.0*atyp          ! A/V = 4pi a^2/(4pi/3 a^3)
      alpha = 1.0                ! sticking probability
      tchemtot = 1.E-99          ! chem.timescale of slowest condensate
      Vol =0.0                   ! total dust volume per cm3
      do i=1,NDUST
        if (eldust(i)<=0.0) cycle
        Vol = Vol + eldust(i)*nHges*dust_Vol(i)      ! [cm3/cm3]
      enddo  
      do i=1,NDUST
        !--- loop over present condensates ---
        if (eldust(i)<=0.0) cycle
c        write(*,*) 'n_'//trim(dust_nam(i))//' = ',eldust(i)*nHges
        !--- find least abundant element in the gas phase   ---
        !--- print stoichiometry and gas element abundances ---
        !--- nkey = min (eps(e )*nHges/s_e )  ---
        nkey = 1.E+99
        do k=1,dust_nel(i)
          e = dust_el(i,k)
          !write(*,'("    Element ",A2,1pE11.4,I2)')
     >    !        elnam(e), eps(e)*nHges, dust_nu(i,k)
          nkeyt = eps(e)*nHges/dust_nu(i,k)
          if (nkeyt<nkey) then
            nkey = nkeyt
            keyel = e
            dustst = dust_nu(i,k)
          endif
        enddo
        !--- find atom/molecule which contains most ---
        !--- of the key element                     ---
        sp = elnam(keyel)
        nmax = nat(keyel)/dustst
        molmass = mass(keyel)
        do imol=1,NMOLE
          included = .false. 
          molm = 0.0
          do j=1,m_kind(0,imol)
            if (m_kind(j,imol)==el) cycle      ! avoid free electrons
            e = elnum(m_kind(j,imol))
            if (e==keyel) then
              included = .true.
              stoich = m_anz(j,imol)
            endif  
            molm = molm + m_anz(j,imol)*mass(e)
          enddo  
          if (included.and.nmol(imol)*stoich/dustst>nmax) then
            sp = cmol(imol)
            nmax = nmol(imol)*stoich/dustst
            molmass = molm
          endif  
        enddo
        vth = SQRT(8.0*bk*Tg/(pi*molmass))     ! [cm/s]
        tchem = 1.d0/(vth*alpha*AoverV*Vol)    ! [s]
c        write(*,'(" limiting element = ",A2)') elnam(keyel)
c        write(*,'(" mostly present as ",A10)') sp
c        write(*,'("        nmax,nkey = ",2(1pE11.4)," cm-3")')
c     >        nmax,nkey
c        write(*,'("             mass = ",2(1pE11.4)," amu")')
c     >        molmass/amu,mass(keyel)/amu
c        write(*,'("        timescale = ",1pE11.4," yr")')
c     >        tchem/yr
        if (tchem>tchemtot) then
          tchemtot = tchem
          limcond  = dust_nam(i)
        endif
      enddo  
c      write(*,'("Limiting condensate ",A22,"  timescale/yr = ",
c     >          1pE11.3)') limcond, tchemtot/yr
      endif
      
c     Modifications added to write a file with all partial
c     pressures for each
c     element and molecule. (28/10/20)

c     file 707 is the one to be read by MARCS
c     it contains partial pressures for a single atmosphere layer
c     i.e. the one which marcs is iterating over atm
      open(unit=707, file='pp.dat', status='replace')
c     first calculate the total number density
c     calculate the contribution of elements:
      n_atom_total = 0.0
      do i=1, el-1
        n_atom_total = n_atom_total + nat(elnum(i))
      enddo
c     calculate the contribution of molecules and ions:
      n_mole_total = 0.0
      do i=1, NMOLE
        n_mole_total = n_mole_total + nmol(i)
      enddo
c     calculate total number density (include number density of
c     electrons! - already defined)
      n_total = nel + n_atom_total + n_mole_total

c     to obtain a partial pressure one multiplies the ratio of the
c     element or molecule number density to the total number density
c     by the total gas pressure
c     write partial pressure of elements:
      write(707,'(1p8e12.3)')(nat(elnum(i)) * pgas/n_total,i=1,el-1)
      write(707,*)
c     write partial pressure of molecules:
      write(707,'(1p8e12.3)')(nmol(i) * pgas/n_total,i=1,NMOLE)
c     write element number and corresponding name:
      write(707,'(i4,a20)')( (elnum(i),elnam(elnum(i))), i=1,el-1)
      write(707,*)
c     write molecule name:
      write(707,714) (cmol(i),i=1,NMOLE)
714   format(10a8)
      close(707)

      pp_electron = (nel/n_total) * pgas
      pp_atoms = (n_atom_total/n_total) * pgas
      pp_molecules = (n_mole_total/n_total) * pgas
      pp_non = 0.0
      pp_sum = pp_electron+pp_atoms+pp_molecules
      open(unit=990,file='GGchem_ppel',status='replace')
      write(990,*) Tg,pgas,pp_electron,mu/amu,rhog,rhod,pp_molecules,
     >  pp_molecules,pp_non,pp_atoms,pp_non,pp_molecules,pp_sum 
      close(990)


      
 1000 format(a6,1pE9.3)
 1010 format(a6,1pE9.3,a8,1pE9.3)
 1020 format(a22,1pE15.9,2(0pF10.5))
 1030 format(a22,1pE12.4)
 4000 format(a7,1pE10.4,a5,1pE10.4)     
 4010 format(' n',a8,1pE12.4,0pF13.9,1pE12.4)
 4020 format(' n',a8,1pE12.4,1pE13.3)
 5000 format(1x,a20,' S=',1pE11.3E4)
      RETURN
      end      
