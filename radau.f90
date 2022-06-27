! 
! This file contains the integrator RADAU by Everhart (1985) with adaptative step.
! The numerical scheme is of order 15.
! 
! A few tips have been taken from Rein & Spiegel (2015) and added to the original routine:
! - the iteration number of the predictor-corrector loop are not predefined but pursued until convergence
! - the choice of stepsize is made using non-dimensional quantities
! 
! This version: M. Saillenfest, 2017.
! Adjustments made in 2022 for compatilibity with c.
! 
! 
!*************************************************************************************************************************

  ! RADAU integrator by Everhart (1985) for equations of type dx/dt = F(x,t)
  ! 
  ! LL                 (integer): Precision required. Typically between 6 (less good) and 14 (best)
  !                               for double precision arithmetics.
  !                               The stepsize is adjusted to respect this precision index.
  !                               If LL <= 0, then dTini is used as a constant time step.
  ! 
  ! dimX               (integer): dimension of the state vector.
  ! 
  ! X   (double precision array): state vector, which is updated during the integration. It contains dimX elements.
  ! 
  ! t         (double precision): time variable, which is updated during the integration.
  ! 
  ! forceC    (function pointer): function that takes X and t as arguments of the differential equation
  !                               and returns F = dX/dt as a double precision array of dimension dimX.
  ! 
  ! dTini     (double precision): initial time step. Its sign gives the direction of the integration.
  !                               If LL <= 0, then dTini is the constant timestep for the whole integration.
  !                               If LL  > 0, then the integration timestep is initialised to dTini and
  !                                           reduced or enlarged to fulfill the precision constraint.
  ! 
  ! Tspan     (double precision): integration timespan. Its sign does not matter.
  ! 
  ! printFunC (function pointer): This subroutine is run after each time step. It can be used to save the
  !                               time evolution of the state vector and other useful quantities.
  ! 
  subroutine RA15_I(LL, dimX, X, t, forceC, dTini, Tspan, printFunC) bind(C, name="RA15_I")

     use,intrinsic:: iso_c_binding ! for interoperability with C
     implicit none

     !----------------------------------- ARGUMENTS DEFINITION ----------------------------------

     integer(c_int),value,intent(in):: LL
     integer(c_int),value,intent(in):: dimX
     real(c_double),intent(inOut):: X(dimX), t
     type(c_funptr),value:: forceC
     real(c_double),value,intent(in):: dTini
     real(c_double),value,intent(in):: Tspan
     type(c_funptr),value:: printFunC

     !-------------------------- PROTOTYPES FOR THE FUNCTIONS PASSED AS ARGUMENT -------------------------

     abstract interface

        subroutine forceFunc(dimXX, XX, tt, FF) bind(C)
           use,intrinsic:: iso_c_binding
           integer(c_int),intent(in):: dimXX
           real(c_double),intent(in):: XX(dimXX), tt
           real(c_double),intent(out):: FF(dimXX)
        end subroutine forceFunc

        subroutine saveFunc(dimXX, XX, tt) bind(C)
           use,intrinsic:: iso_c_binding
           integer(c_int),intent(in):: dimXX
           real(c_double),intent(in):: XX(dimXX), tt
        end subroutine saveFunc

     end interface

     !------------------------------ INTERNAL ADJUSTABLE CONSTANT PARAMETERS ------------------------------

     logical(c_bool),parameter:: forcedConv = .true. ! convergence test of the predictor-corrector loop. If false, the
                                                     ! classic scheme of Everhart (1985) is used: first 6 iterations, then 2.

     logical(c_bool),parameter:: globEst = .true. ! global or local proxy for the convergence test of the predictor-corrector loop.
                                                  ! See Rein & Spiegel (2015).

     real(c_double),parameter:: epsErrEst = 5.0D-10 ! convergence criterion for the predictor-corrector loop.
                                                    ! This sets the relative size of neglected terms in the Taylor series.

     integer(c_int),parameter:: ITMAX = 20 ! maximum number of iterations for the predictor-corrector loop

     !------------------------------ INTERNAL FIXED CONSTANT PARAMETERS ------------------------------

     ! coefficients for the position (Eq. 11)
     real(c_double),parameter:: W(7) = (/ 0.5D0, 1.0D0/3.0D0, 0.25D0, 0.2D0, &
                                            1.0D0/6.0D0, 1.0D0/7.0D0, 0.125D0 /)

     ! coefficients for the velocity (Eq. 12)
     real(c_double),parameter:: U(7) = (/ 0.5D0, 1.0D0/3.0D0, 0.25D0, 0.2D0, &
                                            1.0D0/6.0D0, 1.0D0/7.0D0, 0.125D0 /)

     ! Gauss-Radau spacings over [0,1] (for a numerical scheme of order 15)
     real(c_double),parameter:: H(7) = (/ 0.05626256053692215D0, 0.18024069173689236D0, 0.35262471711316964D0, &
                                            0.54715362633055538D0, 0.73421017721541053D0, 0.88532094683909577D0, &
                                            0.97752061356128750D0 /)
     ! 1D array of Rij = 1/(Hi-Hj)
     integer(c_int):: ii, jj
     real(c_double),parameter:: R(21) = (/  ( (1.0/(H(ii+1)-H(jj)),jj=1,ii), ii=1,6 )  /)

     ! reccurring indexes (for sorting the Cij and Dij into 1D arrays)
     integer(c_int),parameter:: Nw(7) = (/0, 1, 3, 6, 10, 15, 21/)

     ! conversion coefficients from G to B (sorted in a 1D array)
     real(c_double),parameter:: C(21) = (/ -5.62625605369221500D-02,  1.01408028300636303D-02, &
                                             -2.36503252273814510D-01, -3.57589772925161775D-03, &
                                              9.35376952594620664D-02, -5.89127969386984150D-01, &
                                              1.95656540994722115D-03, -5.47553868890686866D-02, &
                                              4.15881200082306861D-01, -1.13628159571753953D+00, &
                                             -1.43653023637089159D-03,  4.21585277212687078D-02, &
                                             -3.60099596502056811D-01,  1.25015071184069102D+00, &
                                             -1.87049177293295006D+00,  1.27179030902686780D-03, &
                                             -3.87603579159067706D-02,  3.60962243452845983D-01, &
                                             -1.46688420840042696D+00,  2.90613625930842930D+00, &
                                             -2.75581271977204583D+00 /)

     ! conversion coefficients from B to G (sorted in a 1D array)
     real(c_double),parameter:: D(21) = (/ 5.62625605369221500D-02, 3.16547571817082965D-03, &
                                             2.36503252273814510D-01, 1.78097769221743422D-04, &
                                             4.57929855060279179D-02, 5.89127969386984150D-01, &
                                             1.00202365223291297D-05, 8.43185715352570116D-03, &
                                             2.53534069054569267D-01, 1.13628159571753953D+00, &
                                             5.63764163931820938D-07, 1.52978400250046571D-03, &
                                             9.78342365324440060D-02, 8.75254664684091088D-01, &
                                             1.87049177293295006D+00, 3.17188154017613784D-08, &
                                             2.76293090982647632D-04, 3.60285539837364600D-02, &
                                             5.76733000277078727D-01, 2.24858876076915978D+00, &
                                             2.75581271977204583D+00 /)


     !-------------------------------- INTERNAL VARIABLES -----------------------------------

     real(c_double):: F(dimX), FJ(dimX), Y(dimX)
     real(c_double):: B(7,dimX), E(7,dimX), BD(7,dimX), G(7,dimX)
     real(c_double):: temp(dimX), B7Old(dimX)

     real(c_double):: TspanV, eps, errEst, t0, Dt, DtOld, q, hj
     logical(c_bool):: firstStep, lastStep
     integer(c_int):: Nit, Ntry1
     integer(c_int):: i, j, jm1, k

     procedure(forceFunc),pointer:: forceF
     procedure(saveFunc),pointer:: printFunF

     ! concerting the c function pointers into fortran procedure pointers
     call c_f_procpointer(forceC,forceF)
     call c_f_procpointer(printFunC,printFunF)

     !-------------------------------- INTEGRATION PARAMETERS -----------------------------------

     TspanV = sign(Tspan,dTini) ! TspanV and dTini must have the same sign
     eps = 10.0D0**(-LL)        ! precision required for the 15-order series expansion

     !-------------------------------------- INITIALISATIONS -----------------------------------------

     BD(:,:) = 0.0D0
     B(:,:) = 0.0D0
     B7Old(:) = B(7,:)

     firstStep = .true. ! true during the first time step
     lastStep = .false. ! true during the last time step

     Ntry1 = 0 ! number of tries to perform the first timestep (if dTini is too large and must be reduced)
     Nit = 6   ! default initial number of iterations for the predictor-corrector loop (6 for the first step, 2 for other steps)

     t0 = t                  ! saving the initial time
     call forceF(dimX,X,t,F) ! force value at initial time
     Dt = dTini              ! next time step

     ! default values if Dt is constant
     DtOld = Dt
     q = 1.0D0

     !----------------------------------- BEGINNING OF THE NUMERICAL INTEGRATION -------------------------------------

     DO ! MAIN LOOP (one cycle = one time step)

        ! computation of the G coefficients using the previous values of the B coefficients (Eq. 7)
        G(1,:) = B(1,:) + D( 1)*B(2,:) + D( 2)*B(3,:) + D( 4)*B(4,:) + D( 7)*B(5,:) + D(11)*B(6,:) + D(16)*B(7,:)
        G(2,:) = B(2,:) + D( 3)*B(3,:) + D( 5)*B(4,:) + D( 8)*B(5,:) + D(12)*B(6,:) + D(17)*B(7,:)
        G(3,:) = B(3,:) + D( 6)*B(4,:) + D( 9)*B(5,:) + D(13)*B(6,:) + D(18)*B(7,:)
        G(4,:) = B(4,:) + D(10)*B(5,:) + D(14)*B(6,:) + D(19)*B(7,:)
        G(5,:) = B(5,:) + D(15)*B(6,:) + D(20)*B(7,:)
        G(6,:) = B(6,:) + D(21)*B(7,:)
        G(7,:) = B(7,:)

        !----------------------------- MOVING FORWARD BY ONE TIMESTEP Dt --------------------------------

        do i = 1,ITMAX ! predictor-corrector loop

           !***************** 7 substeps using the Gauss-Radau spacing *****************
           do j = 1,7

              hj = H(j)
              jm1 = j - 1

              ! ---- estimates for the position and the force at this substep (Eq. 9) ----
              Y(:) = X(:) + Dt*hj*( F(:) + hj*( W(1)*B(1,:) + hj*( W(2)*B(2,:) &
                                         + hj*( W(3)*B(3,:) + hj*( W(4)*B(4,:) &
                                         + hj*( W(5)*B(5,:) + hj*( W(6)*B(6,:) &
                                         + hj*W(7)*B(7,:))))))))

              call forceF(dimX, Y, t+hj*Dt, FJ) ! updates FJ(:)

              ! ---- values of G for the force FJ of this substep (Eq. 4) ----

              temp = G(j,:) ! saving the previous value of G(j,:)
              G(j,:) = (FJ(:) - F(:))/hj ! initialisation

              do k = 1,jm1
                 G(j,:) = ( G(j,:) - G(k,:) )*R(Nw(jm1)+k)
              endDo
              temp = G(j,:) - temp ! improvement obtained for the value of G(j,:)

              ! ---- corresponding improvement of the values of B (Eq. 5) ----
              do k = 1,jm1
                 B(k,:) = B(k,:) + C(Nw(jm1)+k)*temp
              endDo
              B(j,:) = B(j,:) + temp ! last term

           endDo
           ! ***************************** end of the 7 substeps ********************************

           if(forcedConv) then ! convergence criterion

              !----------------- convergence test of the predictor-corrector ------------------
              if(globEst) then
                 errEst = maxval(abs(B(7,:)-B7Old(:))) / maxval(abs(F(:))) ! global estimate
              else
                 errEst = maxval(abs( (B(7,:)-B7Old(:))/F(:) )) ! local estimate
              endIf
              B7Old(:) = B(7,:) ! update of B7Old

              !write(*,*) "RA15_I pred-corr : ", i, errEst
              if(errEst < epsErrEst) exit ! convergence achieved

           else if(i >= Nit) then ! classic method of Everhart (1985): 6 iterations for warm up, then 2 iterations
              exit
           endIf

        endDo ! end of the predictor-corrector loop

        !if(i >= ITMAX) then
        !   write(0,*) "WARNING IN RA15_I: convergence criterion hard to reach for the predictor-corrector."
        !   write(0,*) "Precision reached/required: ", errEst, epsErrEst
        !endIf

        !----------------------------- CHOICE OF THE NEXT TIMESTEP --------------------------------

        if(LL > 0) then ! if a given precision is required

           DtOld = Dt

           Dt = sqrt9( eps * 72.0D0*Dt**7 / maxval(abs(B(7,:))) ) ! next timestep adapted to the required precision
           q = Dt/DtOld

           ! PARTICULAR CASE OF THE INITIAL STEP
           ! We consider that the required precision is not reached if the second step is smaller than the first one.
           if(firstStep) then
              if(q < 1.0D0)  then
                 if(Ntry1 > 20) then
                    !write(0,*) "FATAL ERROR IN RA15_I: the initial time step is too large for the required precision."
                    !stop
                    return
                 endIf
                 Ntry1 = Ntry1 + 1
                 Dt = 0.8D0*Dt ! the new step is 1.25 times smaller than the new estimated timestep
                 cycle ! we start all over again
              else
                 !write(*,'(A,I0)') "Number of rectifications of the first step: ", Ntry1
              endIf
           endIf

           ! restricting the size of the next step (we avoid overly large variations of Dt for safety)
           if(q > 1.4D0) then
              Dt = DtOld*1.4D0
              q = 1.4D0
           endIf

        endIf

        !------------------------------- UPDATE OF THE STATE VECTOR, TIME, FORCE -------------------------------------

        ! position and velocity (Eqs. 11 and 12)
        ! As the time step has been adjusted, the step previously made is DtOld.
        X(:) = X(:) + DtOld*( F(:) + B(1,:)*W(1) + B(2,:)*W(2) + B(3,:)*W(3) + B(4,:)*W(4) &
                                + B(5,:)*W(5) + B(6,:)*W(6) + B(7,:)*W(7) )

        t = t + DtOld ! time
        call forceF(dimX, X, t, F) ! force: updates F(:)

        !------------------------------------ STOP OF THE INTEGRATION -------------------------------------------

        if(lastStep) then ! if this was the last step
           call printFunF(dimX,X, t) ! last data save
           exit
        endIf

        !--------------------------- POSSIBLE RE-ADJUSTMENT OF THE NEXT STEP -------------------------------

        ! If the next step will be the last one
        if( abs(t+Dt-t0) >= abs(TspanV) ) then
           if(LL>0) then
              Dt = TspanV - (t-t0) ! readjusting the last step in order to end exactly after Tspan
              q = Dt/DtOld
           endIf
           lastStep = .true.
        endIf

        !-------------------------------- DATA SAVE --------------------------------------

        call printFunF(dimX,X, t)

        !----------------------- PREDICTION OF THE B COEFFICIENTS FOR THE NEXT TIMESTEP ---------------------------

        if(firstStep) then
           firstStep = .false.
           Nit = 2 ! we switch to a default value of 2 iterations for the predictor-corrector loop
        else
           BD(:,:) = B(:,:) - E(:,:) ! difference between B and its estimate E for this timestep
        endIf

        ! estimate of the next values that the B coefficients will have (Eq. 13)
        E(1,:) = (B(1,:)+2.0D0*B(2,:)+ 3.0D0*B(3,:)+ 4.0D0*B(4,:)+ 5.0D0*B(5,:)+ 6.0D0*B(6,:)+7.0D0*B(7,:))*q
        E(2,:) = (B(2,:)+3.0D0*B(3,:)+ 6.0D0*B(4,:)+10.0D0*B(5,:)+15.0D0*B(6,:)+21.0D0*B(7,:))*q**2
        E(3,:) = (B(3,:)+4.0D0*B(4,:)+10.0D0*B(5,:)+20.0D0*B(6,:)+35.0D0*B(7,:))*q**3
        E(4,:) = (B(4,:)+5.0D0*B(5,:)+15.0D0*B(6,:)+35.0D0*B(7,:))*q**4
        E(5,:) = (B(5,:)+6.0D0*B(6,:)+21.0D0*B(7,:))*q**5
        E(6,:) = (B(6,:)+7.0D0*B(7,:))*q**6
        E(7,:) =  B(7,:)*q**7

        ! a more precise estimate is obtained by adding the correction obtained for the current step
        B(:,:) = E(:,:) + BD(:,:)

        !---------------------------------- WE GO TO THE NEXT TIMESTEP --------------------------------------
     ENDdO
     !  END OF THE NUMERICAL INTEGRATION

  contains
  !-----------------------------------------------------------------
  ! power 1/9 accepting a negative argument
  !-----------------------------------------------------------------
     function sqrt9(z)
        real(c_double):: sqrt9
        real(c_double),intent(in):: z
        real(c_double),parameter:: pw = 1.0D0/9.0D0

        sqrt9 = sign( abs(z)**pw , z)

     end function sqrt9
  !-----------------------------------------------------------------
  end subroutine RA15_I

!*************************************************************************************************************************
