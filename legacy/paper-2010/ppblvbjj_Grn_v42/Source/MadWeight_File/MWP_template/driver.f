      PROGRAM DRIVER
C**************************************************************************
C     THIS IS THE DRIVER FOR THE WHOLE CALCULATION
C
C**************************************************************************
C**************************************************************************
C
C    Structure of the driver
C
C       1) initialisation
C       2) first loop
C           2.1) loop on permutation
C           2.2) loop on configuration
C       3) second loop
C           3.1) loop on permutation
C           3.2) loop on configuration
C           3.3) write result
C          
C**************************************************************************
      IMPLICIT NONE
C
      include 'genps.inc'
      include 'nexternal.inc'
      include 'data.inc'
      include 'madweight_param.inc'
      include 'driver.inc'
C
C     PARAMETERS
C
C
C     LOCAL (Variable of integration)
C
      DOUBLE PRECISION SD,CHI,CROSS
      INTEGER I,convergence_status
      integer pos,channel_pos
c
c     variable for permutation, channel
c
      integer perm_id(nexternal-2)     !permutation of 1,2,...,nexternal-2
      integer num_per                  !total number of permutations
      double precision weight,weight_error
c
c     variable for selecting perm in loop 1
c
      double precision minimal_value1,minimal_value2
      double precision best_precision
c      integer best_config
c
c     logical
c
      logical debug
C
C     GLOBAL     (BLOCK IN COMMON WITH VEGAS)
C
      DOUBLE PRECISION Xl(20),XU(20),ACC
      INTEGER NDIM,NCALL,ITMX,NPRN
      COMMON/BVEG1/XL,XU,ACC, NDIM,NCALL,ITMX,NPRN


C**************************************************************************
C    step 1: initialisation     
C**************************************************************************
      call global_init

      if(use_perm) then
         call get_num_per(num_per) !calculate the number of permutations.
      else
         num_per=1
      endif
C
C       VARIABLE INITIALISATION
C
 1    call init_perm()
      do i=1,nexternal-2
         perm_id(i)=i
      enddo
      if ((num_per*nb_sol_config).gt.max_channel) then
         write(*,*) 'ERROR: TOO MANY CHANNEL OF INTEGRATION'
         STOP
      endif
      do channel_pos=1,num_per*nb_sol_config
         cross_perm(1,channel_pos)=0d0
         cross_perm(2,channel_pos)=9d99
         l_cross_max(channel_pos)=0
      enddo
      if (histo)   call histo_init(num_per*nb_sol_config)

C**************************************************************************
C    step 2: first loop -ordering  channel-    
C**************************************************************************     
      channel_pos=0
      do perm_pos=1,num_per
         call give_permut(perm_id)     
         call assign_perm(perm_id)
         do config_pos=1,nb_sol_config
            channel_pos=channel_pos+1
            call integral_estimation(channel_pos)
            enddo               ! end config run  
         enddo                  ! end permutation run

      if (cross_perm(1,l_cross_max(1)).eq.0d0.or.
     &cross_perm(2,l_cross_max(1))/cross_perm(1,l_cross_max(1)).gt.
     &min_prec_cut1)then
         nevents=5*nevents
         if (nevents.le.8d6) then
            write(*,*) "restart with more points"
            goto 1              !restart program
         elseif (nevents.gt.8d6) then
            nevents = 5d5
         endif
      endif

C**************************************************************************
C    step 3: prepare the grid and compute the weight     
C**************************************************************************
      write(*,*) '\nstart third step\n'

      weight=0d0
      weight_error=0d0
      do pos=1,num_per*nb_sol_config
         !compute position object
         channel_pos=l_cross_max(pos)
         perm_pos=((channel_pos-1)/nb_sol_config)+1   ! integer division performed
         config_pos=(channel_pos-1)-(perm_pos-1)*nb_sol_config+1

         !check if we have to continue or not (precision reached)
         minimal_value1=cross_perm(2,l_cross_max(1))*value_cut1 
         if(cross_perm(1,channel_pos).lt.minimal_value1)then 
            write(*,*) "other channel are too law to be important "
            write(*,*) pos,perm_pos,cross_perm(1,perm_pos),
     &                  minimal_value1
            goto 6             ! write result
         endif

         !init this permutation
         call init_perm()
         do i=1,perm_pos
            call give_permut(perm_id)     
         enddo
         call assign_perm(perm_id)
         write(*,*) 'permutation ',perm_pos,' config ',
     &config_pos,' : ',cross_perm(1,channel_pos)  

         !make the integration
         call final_integration_launcher(channel_pos,
     &                       cross_perm(2,l_cross_max(1)))
      enddo

C**************************************************************************
C    step 6: write file     
C**************************************************************************

 6    weight=0d0
      weight_error=0d0
      convergence_status=1
      do channel_pos=1,num_per*nb_sol_config
         weight=weight+cross_perm(1,channel_pos)
         if (cross_perm(1,channel_pos).gt.0d0) weight_error=
     & weight_error+cross_perm(2,channel_pos)**2
         if (cross_perm(2, channel_pos).lt.0d0) convergence_status=-1
      enddo
      weight_error=dsqrt(weight_error)

      OPEN(UNIT=21,FILE='./weights.out',STATUS='UNKNOWN')
C      OPEN(UNIT=22,FILE='./likelihood.out',STATUS='UNKNOWN')
      OPEN(UNIT=23,FILE='./details.out',STATUS='UNKNOWN')

      write(21,*)weight/num_per,'  ',
     &convergence_status*weight_error/num_per
      write(*,*) "======================================"
      write(*,*)'Weight: ',weight/dble(num_per),'+-',
     &     weight_error/dble(num_per)
      write(*,*) "======================================"

      write(23,*) ' permutation channel   value      error'
      call init_perm()
      channel_pos=0
      do perm_pos=1,num_per
         call give_permut(perm_id)
         write(23,*) "======================================"
         write(23,*) perm_id
         do config_pos=1,nb_sol_config
            channel_pos=channel_pos+1
            write(23,*) perm_pos,' || ',config_pos,' || ',
     & cross_perm(1,channel_pos),
     &           ' ||',cross_perm(2,channel_pos)
         enddo
      enddo
      write(23,*) "======================================"
      write(23,*)'Weight: ',weight/dble(num_per),'+-',
     &     weight_error/dble(num_per)
      write(23,*) "======================================"
      close(21)
      OPEN(UNIT=21,FILE='./stop',STATUS='UNKNOWN')
      close(21)
C**************************************************************************
C    step 7: write histogram file (not by default) 
C**************************************************************************
      if (histo)   call histo_final(num_per*nb_sol_config)

      END

C**************************************************************************************
C                             related function
C**************************************************************************************

      subroutine integral_estimation(channel_pos)
      implicit none
      include 'genps.inc'
      include 'nexternal.inc'
      include 'data.inc'
      include 'madweight_param.inc'
      include 'driver.inc'
c
c     this subroutine call vegas for the phase space integration
c
c
c     local
c
      double precision cross,sd,chi
      integer i,vegas_flag
      integer channel_pos
C
C     GLOBAL     (BLOCK IN COMMON WITH VEGAS)
C
      DOUBLE PRECISION Xl(20),XU(20),ACC
      INTEGER NDIM,NCALL,ITMX,NPRN
      COMMON/BVEG1/XL,XU,ACC, NDIM,NCALL,ITMX,NPRN
      double precision calls,ti,tsi
      COMMON/BVEG4/CALLS,TI,TSI
C
C     Global  number of dimensions
C
      integer Ndimens
      common /to_dimension/ Ndimens
c
c      external
c
      double precision fct
      external  fct
c
c
c
      integer          iseed
      common /to_seed/ iseed
c-----------
c Begin code
c-----------
      write(*,*) 'driver enter in integral first estimation'
      write(*,*) 'Ndimens', Ndimens
      call initialize

      iseed=iseed+1 !avoid all integration to have the same seed
      NDIM=Ndimens
      ACC  = min_prec_cut1 ! by default 0.5
      NPRN = 1
      NCALL= nevents / 10
      ITMX = max_it_step1  ! by default 3
      write(*,*) 'driver before first vegas call',nevents

      cross=0d0
      do i=1,2              
         call check_nan(CROSS)
         if (CROSS.eq.0d0) then
             NCALL=10*NCALL
             CALL VEGAS(fct,CROSS,SD,CHI)
             if (histo) call histo_combine_iteration(channel_pos)
             write(32,*) '1',perm_pos,config_pos,CALLS,ti,tsi,cross,sd
         endif
      enddo

      CALL check_nan(cross)
      if (cross.eq.0d0) sd=9d99
            
C store information and ordering the crontribution
      call order_cross(channel_pos, cross, sd)

      end

C**************************************************************************************

      subroutine final_integration_launcher(channel_pos, err_crossmax)
      implicit none
      include 'genps.inc'
      include 'nexternal.inc'
      include 'data.inc'
      include 'madweight_param.inc'
      include 'driver.inc'
c
c     this subroutine call vegas for the final phase space integration
c
c
c     local
c
      double precision cross,sd,chi,err_crossmax
      integer pos,channel_pos,i,j
      double precision cut_on_val,cut_on_err,old_error,old_cross

C
C     GLOBAL     (BLOCK IN COMMON WITH VEGAS)
C
      DOUBLE PRECISION Xl(20),XU(20),ACC
      INTEGER NDIM,NCALL,ITMX,NPRN
      COMMON/BVEG1/XL,XU,ACC, NDIM,NCALL,ITMX,NPRN
      Double precision calls,ti,tsi
      COMMON/BVEG4/calls,ti,tsi
      Double precision ALPH
      integer NDMX,MDS
      COMMON/BVEG3/ALPH,NDMX,MDS
C
C     Global  number of dimensions
C
      integer Ndimens
      common /to_dimension/ Ndimens
c
c      external
c
      double precision fct
      external  fct
c
c
c
      integer          iseed
      common /to_seed/ iseed
c-----------
c Begin code
c-----------
c      write(*,*) 'driver enter in final integration launcher !!'
c      write(*,*) 'Ndimens',Ndimens
c
c     no invisible particle
c
      iseed=iseed+1
      NDIM=Ndimens
      call initialize
c
c     compute the end condition (relative precision - absolute value - absolute error)      
c
      ACC=final_prec
      if (pos.eq.1) then
         cut_on_val=1d-99
         cut_on_err=1d-99
      else
         cut_on_val=err_crossmax*value_cut1
         cut_on_err=err_crossmax*value_cut2
      endif
c
c     compute the a more precise estimation
c
      NPRN = 1
      NCALL = nevents
      ITMX = 1
c      ALPH=1.5d0
      write(*,*) 'driver before first vegas call', nevents
      CALL VEGAS(fct,CROSS,SD,CHI)
      if (histo) call histo_combine_iteration(perm_pos*nb_sol_config+(config_pos-1))
      write(32,*) '3',perm_pos,config_pos,CALLS,ti,tsi,cross,sd
      do i=1,2
         CALL check_nan(cross)
         if (cross.eq.0d0) then
            NCALL=10*nevents
            CALL VEGAS(fct,CROSS,SD,CHI)
            if (histo) call histo_combine_iteration(perm_pos*nb_sol_config+(config_pos-1))
            write(32,*) '3',perm_pos,config_pos,CALLS,ti,tsi,cross,sd
         endif
      enddo
      CALL check_nan(cross)
      if (cross.eq.0d0) then
         sd=9d99
      endif

      old_error=sd
      old_cross=cross
c      ITMX = max_it_step3
      do i=2, max_it_step3       
         ITMX=ITMX+1
         if (cross.le.cut_on_val
     &        .or.sd.le.cut_on_err
     &        .or.(sd/cross).le.final_prec) then
           goto 42 ! end of this integration
         endif
         CALL VEGAS2(fct,CROSS,SD,CHI)
         if (histo) call histo_combine_iteration(perm_pos*nb_sol_config+(config_pos-1))
         write(32,*) '4',perm_pos,config_pos,CALLS,ti,tsi,cross,sd
         if(sd.gt.(10*old_error).and.old_error/old_cross.le.3d-1)then
            write(*,*) 'detection of overconstrained grid 
     &       relaunch point with less iteration in MW_run 733'
            sd=-1d0*sd
            goto 42 ! end of this estimation
         else
            old_error=sd
            old_cross=cross
         endif
      enddo
      if (pos.eq.1) then

         NCALL=int(nevents*max_it_step3*(sd/cross)**2/(final_prec**2*
     &          (max_it_step5-max_it_step3))*1.2d0)
      else
         if (final_prec.ge.cut_on_err/cross)then
            NCALL=int(nevents*max_it_step3*(sd/cross)**2/(final_prec**2*
     &           (max_it_step5-max_it_step3))*1.2d0)
         else
            NCALL=int(1.2d0 * nevents*max_it_step3*(sd/cross)**2
     &           / (cut_on_err/cross)**2
     &           /(max_it_step5-max_it_step3))
         endif
      endif
      if (NCALL.lt.0d0.or.NCALL.gt.5d7) then
         write(*,*) 'Warning Huge number of point launched by ',
     &  'integration... (its better to increase the number of ',
     &   'iteration)'
         NCALL=20000000
      endif
c      ALPH=0d0
      write(*,*) 'boost NCALL',NCALL
c      write(*,*) 'info',max_it_step3,max_it_step5,final_prec

      do i=1,max_it_step5       
         ITMX = ITMX+1
         if (cross.le.cut_on_val
     &       .or.sd.le.cut_on_err
     &       .or.(sd/cross).le.final_prec) then
            goto 42 ! end of this integration
         endif
             CALL VEGAS2(fct,CROSS,SD,CHI)
         if (histo) call histo_combine_iteration(perm_pos*nb_sol_config+(config_pos-1))
         write(32,*) '5',perm_pos,config_pos,CALLS,ti,tsi,cross,sd
         if(sd.gt.(10*old_error)
     &      .and.old_error/old_cross.le.3d-1)then
            write(*,*) 'detection of overconstrained grid
     &       relaunch point with less iteration in MW_run 733'
            sd=-1d0*sd
            goto 42 !return
         else
            old_error=sd
            old_cross=cross
         endif
      enddo
      

 42   cross_perm(1,channel_pos)=cross
      cross_perm(2,channel_pos)=sd




   
         end

C**************************************************************************************
      subroutine order_cross(pos,cross,sd)
      implicit none
      double precision cross,sd        ! new_input 
      integer pos                      ! current position in cross_perm
      integer i,j
      include 'driver.inc'

c     assign in cross_perm
      cross_perm(1,pos)=cross
      cross_perm(2,pos)=sd

c     assign order

      do i=1,pos
         if( cross_perm(1,l_cross_max(i)).lt.cross)then
            do j=pos,i,-1
               l_cross_max(j+1)=l_cross_max(j)
            enddo
            l_cross_max(i)=pos
            goto 10 !break
         endif
      enddo
 10   return
      end


C**************************************************************************************
        double precision function fct(x,wgt)
        implicit none

        include 'genps.inc'
        include 'nexternal.inc'
        include 'run.inc'
        include 'driver.inc'
c
c       this is the function which is called by the integrator
c
c       arguments
c
        double precision x(20),wgt
c
c       local
c
c        integer i,j ! debug mode
c
c       global
c
        double precision              S,X1,X2,PSWGT,JAC
        common /PHASESPACE/ S,X1,X2,PSWGT,JAC
        double precision momenta(0:3,-max_branch:2*nexternal)  ! records the momenta of external/intermediate legs     (MG order)
        double precision mvir2(-max_branch:2*nexternal)                  ! records the sq invariant masses of intermediate particles (MG order)
        common /to_diagram_kin/ momenta, mvir2


        logical histo
        common /to_histo/histo
c
c       external
c
        double precision dsig
        external dsig
        include 'data.inc'

         call get_PS_point(x)
         if (jac.gt.0d0) then
         fct=jac
         xbk(1)=X1
         xbk(2)=X2
         fct=fct*dsig(momenta(0,1),wgt)
         else
         fct=0d0
         return
         endif
         
         if (histo)  then
            call FILL_plot(fct,wgt,perm_pos*nb_sol_config+config_pos)
         endif
         end


C**************************************************************************************
      subroutine global_init
      implicit none
C
C     initialize all global variable
C
      include 'genps.inc'
      include 'nexternal.inc'
c
c     LHCO input
c
      integer tag_init(3:nexternal),type(nexternal),run_number,trigger
      double precision eta_init(nexternal),phi_init(nexternal),
     &pt_init(nexternal),j_mass(nexternal),ntrk(nexternal),
     &btag(nexternal),had_em(nexternal),dummy1(nexternal),
     &dummy2(nexternal)
      common/LHCO_input/eta_init,phi_init,pt_init,
     &j_mass,ntrk,btag,had_em,dummy1,dummy2,tag_init,type,run_number,
     &trigger

      integer MG,j,k,i

      integer matching_type_part(3:nexternal) !modif/link between our order by type for permutation
      integer inv_matching_type_part(3:nexternal)
      common/madgraph_order_type/matching_type_part,
     & inv_matching_type_part

C     topol info
      integer num_inv,num_jet,num_bjet,num_e,num_ae,num_mu,num_amu,
     +        num_ta,num_ata   !number of jet,elec,muon, undetectable
      COMMON/num_part/num_inv,num_jet,num_bjet,num_e,num_ae,
     & num_mu,num_amu,num_ta,num_ata !particle in the final state

      double precision pexp_init(0:3,nexternal)  !impulsion in original configuration
      common/to_pexp_init/pexp_init

      character*60 param_name

      double precision pmass(nexternal)
      common/to_mass/pmass
      double precision momenta(0:3,-max_branch:2*nexternal)  ! records the momenta of external/intermediate legs     (MG order)
      double precision mvir2(-max_branch:2*nexternal)        ! records the sq invariant masses of intermediate particles (MG order)
      common /to_diagram_kin/ momenta, mvir2

      integer met_lhco,opt_lhco
      common/LHCO_met_tag/met_lhco,opt_lhco

      integer nevents
      common/number_integration_points/nevents
c
c     var to control histo
c
      logical histo
      common /to_histo/histo

      integer          isum_hel
      logical                   multi_channel
      common/to_matrix/isum_hel,multi_channel

      multi_channel=.false.

      open(unit=90,file="./start")
      write(90,*) 'start'
      close(90)
      open(UNIT=32,FILE='vegas_value.out',STATUS='unknown')
c
c     set parameters for the run
c
      CALL MW_SETRUN
      open(unit=89,file="./param.dat")
      read(89,*) param_name
c      read(89,*) nevents
      close(89)
      write(*,*) 'param name',param_name
      call setpara(param_name,.true.)
      write(*,*) 'number of points',nevents
c
c    set cuts
c
      CALL SETCUTS
c
c     deternime the final state content
c
      call info_final_part()                       !init variable concerning outut particle (num/type)+ give order
c
c     set parameters for the transfer functions
c
c      write(*,*) 'test 1'
      call setTF
c      write(*,*) 'test 2'
      call init_d_assignement()
      CALL PRINTOUT
      CALL RUN_PRINTOUT
c      call graph_init

      write(*,*) "read verif"
      OPEN(UNIT=24,file='./verif.lhco',status='old',err=48) ! input file

      read(24,*,err=48,end=48) k,run_number,trigger !line with tag 0
      do j=3,nexternal-num_inv
         MG=inv_matching_type_part(j)
         read(24,*,err=48,end=48) k,type(k),eta_init(k),
     &           phi_init(k),pt_init(k),j_mass(k),ntrk(k),btag(k),
     &           had_em(k) ,dummy1(k),dummy2(k)

         call four_momentum_set2(eta_init(k),phi_init(k),pt_init(k),
     &         j_mass(k),pexp_init(0,MG))
         tag_init(MG)=k
      enddo
      read(24,*,err=48,end=48) k,type(k),eta_init(k), !line with type 6
     &           phi_init(k),pt_init(k),j_mass(k),ntrk(k),btag(k),
     &           had_em(k) ,dummy1(k),dummy2(k)
      read(24,*,err=47,end=47) k,type(k),eta_init(k), !line with type 7 (optional)        
     &           phi_init(k),pt_init(k),j_mass(k),ntrk(k),btag(k),  
     &           had_em(k) ,dummy1(k),dummy2(k)
      opt_lhco=k
      close(24)
 47   write(*,*) "read experimental data:"
      write(*,*) "-----------------------"
      do j=3,nexternal-num_inv
         MG=inv_matching_type_part(j)
         write(*,*) (pexp_init(i,MG),i=0,3)
      enddo
c      call set_pmass ! done in setcuts
      do i=3,nexternal
         mvir2(i)=pmass(i)**2
      enddo
c      write(*,*) 'pmass ok'
      return
 48   stop
      end

C**************************************************************************************
C      subroutine save_grid(perm_pos,config)
C
C
C
C      integer perm_pos,config
C      character*11 buffer
C
C     GLOBAL     (BLOCK IN COMMON WITH VEGAS)
C
C      DOUBLE PRECISION Xl(20),XU(20),ACC
C      INTEGER NDIM,NCALL,ITMX,NPRN
C      COMMON/BVEG1/XL,XU,ACC, NDIM,NCALL,ITMX,NPRN
C      DOUBLE PRECISION              S,X1,X2,PSWGT,JAC
C      COMMON /PHASESPACE/ S,X1,X2,PSWGT,JAC
C      integer*4 it,ndo
C      double precision xi(50,20),si,si2,swgt,schi
C      common/bveg2/xi,si,si2,swgt,schi,ndo,it   

C      buffer='grid_00_000'
C      if (config.ge.10)then
C         write(buffer(6:7),'(I2)') config
C      else
C         write(buffer(7:7),'(I1)') config
C      endif
C      if (perm_pos.ge.100)then
C         write(buffer(9:11),'(i3)')perm_pos
C      elseif(perm_pos.ge.10)then
C         write(buffer(10:11),'(i2)')perm_pos
C      else
C         write(buffer(11:11),'(i1)')perm_pos
C      endif
C         write(*,*) buffer         

C      open(unit=88,file=buffer)
C      write(88,100)ndo
C      do i=1,ndo
C         write(88,101) (xi(i,j),j=1,20)
C      enddo
C      close(88)



C 100  format(5X,I5)
C 101  format(5x,E11.5,2x,E11.5,2x,E11.5,2x,E11.5,2x,E11.5,
C     &       2x,E11.5,2x,E11.5,2x,E11.5,2x,E11.5,2x,E11.5,
C     &       2x,E11.5,2x,E11.5,2x,E11.5,2x,E11.5,2x,E11.5,
C     &       2x,E11.5,2x,E11.5,2x,E11.5,2x,E11.5,2x,E11.5)
C 102  format(A5,I2,A1,I3)
c      end


C**************************************************************************************
C      subroutine read_grid(perm_pos,config)
C
C
C
C      integer perm_pos,config
C      character*20 buffer
C
C     GLOBAL     (BLOCK IN COMMON WITH VEGAS)
C
C      DOUBLE PRECISION Xl(20),XU(20),ACC
C      INTEGER NDIM,NCALL,ITMX,NPRN
C      COMMON/BVEG1/XL,XU,ACC, NDIM,NCALL,ITMX,NPRN
C      DOUBLE PRECISION              S,X1,X2,PSWGT,JAC
C      COMMON /PHASESPACE/ S,X1,X2,PSWGT,JAC
C      integer*4 it,ndo
C      double precision xi(50,20),si,si2,swgt,schi
C      common/bveg2/xi,si,si2,swgt,schi,ndo,it   

C      buffer='grid_00_000'
C      if (config.ge.10)then
C         write(buffer(6:7),'(I2)') config
C      else
C         write(buffer(7:7),'(I1)') config
C      endif
C      if (perm_pos.ge.100)then
C         write(buffer(9:11),'(i3)')perm_pos
C      elseif(perm_pos.ge.10)then
C         write(buffer(10:11),'(i2)')perm_pos
C      else
C         write(buffer(11:11),'(i1)')perm_pos
C      endif
C         write(*,*) buffer         
C      open(unit=88,file=buffer)
C      READ(88,100) ndo
C      do i=1,config
C         READ(88,101)(xi(i,j),j=1,20)
C      enddo
C      close(88) 


C 100  format(5X,I5)
C 101  format(5x,E11.5,2x,E11.5,2x,E11.5,2x,E11.5,2x,E11.5,
C     &       2x,E11.5,2x,E11.5,2x,E11.5,2x,E11.5,2x,E11.5,
C     &       2x,E11.5,2x,E11.5,2x,E11.5,2x,E11.5,2x,E11.5,
C     &       2x,E11.5,2x,E11.5,2x,E11.5,2x,E11.5,2x,E11.5)
C 102  format(A5,I2,A1,I3)
C
C      end
C**************************************************************************************
