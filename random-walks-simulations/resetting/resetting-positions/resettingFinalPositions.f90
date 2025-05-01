include 'resettingPositionsModule.f90'
program resettingFinalPositions
   !!  ------------------------------------------------------------------------
   !!  |  SIMULATION PROGRAM FOR THE RESETTING RANDOM WALKER FINAL POSITIONS  |
   !!  ------------------------------------------------------------------------
   !!
   use omp_lib
   use resettingPositionsModule
   implicit none

   !! VARIABLES TO MEASURE EXECUTION TIME
   integer :: t0, t1, CountRate;

   !! FIRST CHECKPOINT
   call system_clock(t0);

   !! INITIALIZING LISTS AS ZERO
   s1 = 0.; s2 = 0.; s3 = 0.; s4 = 0.;
   skw = 0.; krt = 0.; var = 0.;

   !! POSITIONS LIST
   do i=1,2*maxDistance+1
      positions(i) = i-maxDistance-1
   end do

   !! INITIALIZING MOMENTS COUNTS
   allocate(momentCounts(momentsTime,2*maxDistance+1));

   !! SIMULATION OF THE RANDOM WALKS
   do i=1,cases
      probabilityTemp = probability(i);
      momentCounts = 0;

      !$OMP PARALLEL PRIVATE(randomWalk,tempRand,ii,iii,j) &
      !$OMP& REDUCTION(+:momentCounts, counts)

      !!  RANDOM SEED INITIALIZER
      call init_random_seed();
      j = steps;
      
      !$OMP DO
      do ii=1,realizations
         randomWalk = 0;

         if (ii .ge. int(realizations/2)) then
            j = steps+1;
         end if

         do iii=1,j

            if (iii<=momentsTime.and.abs(randomWalk)<=maxDistance) then
               momentCounts(iii,randomWalk+maxDistance+1) = momentCounts(iii,randomWalk+maxDistance+1) + 1;
            end if

            !! STEP DYNAMICS
            call random_number(tempRand);
            if (tempRand <= probabilityTemp) then
               randomWalk = 0;
            else
               call random_number(tempRand);
               randomWalk = randomWalk + 2*nint(tempRand)-1;
            end if

         end do

         !! COUNTING THE FINAL POSITIONS
         if (abs(randomWalk)<=maxDistance) then
            counts(i,randomWalk+maxDistance+1) = counts(i,randomWalk+maxDistance+1) + 1
         end if

      end do
      !$OMP END DO
      !$OMP END PARALLEL

      !! CALCULATING THE MOMENTS
      do ii=1,momentsTime

         do j=1,2*maxDistance+1

            s1(i,ii) = s1(i,ii) + 1.*momentCounts(ii,j)*(positions(j))/realizations
            s2(i,ii) = s2(i,ii) + 1.*momentCounts(ii,j)*(positions(j)**2.)/realizations
            s3(i,ii) = s3(i,ii) + 1.*momentCounts(ii,j)*(positions(j)**3.)/realizations
            s4(i,ii) = s4(i,ii) + 1.*momentCounts(ii,j)*(positions(j)**4.)/realizations

         end do

         var(i,ii) = ( s2(i,ii) - s1(i,ii)**2. )
         skw(i,ii) = ( s3(i,ii) - 3.*s1(i,ii)*( var(i,ii) )&
                  - s1(i,ii)**3. )/( var(i,ii)**(3./2.) )
         krt(i,ii) = ( s4(i,ii) - 4.*s1(i,ii)*s3(i,ii) &
                  + 6.*(s1(i,ii)**2.)*s2(i,ii) - 3.*(s1(i,ii)**4.))/&
                  (var(i,ii)**2.)

      end do
   end do

   deallocate(momentCounts);

   !!  SECOND CHECKPOINT
   call system_clock(t1,CountRate);
   print*,"Execution Time: ",(t1-t0)*1./CountRate;


   !! DATA STORING SECTION
   !! --------------------

   !! DISTRIBUTION
   open(1, file='../../../data/position-distribution/distributions/resetting_distribution_pos.txt', status='replace')
      do i=1,2*maxDistance+1
            write(1,*) counts(1,i), counts(2,i), counts(3,i), counts(4,i);
      end do
   close(1)

   !! CASES
   open(1, file='../../../data/position-distribution/information/resetting_cases_pos.txt', status='replace')
      do i=1,cases
         write(1,*) probability(i)
      end do
   close(1)

   !! VARIANCE
   open(1, file='../../../data/position-distribution/statistics/resetting_variance_pos.txt', status='replace')
      do i=1,momentsTime
         write(1,*) var(1,i), var(2,i), var(3,i), var(4,i);
      end do
   close(1)

   !! SKEWNESS
   open(1, file='../../../data/position-distribution/statistics/resetting_skewness_pos.txt', status='replace')
      do i=1,momentsTime
         write(1,*) skw(1,i), skw(2,i), skw(3,i), skw(4,i);
      end do
   close(1)

   !! KURTOSIS
   open(1, file='../../../data/position-distribution/statistics/resetting_kurtosis_pos.txt', status='replace')
      do i=1,momentsTime
         write(1,*) krt(1,i), krt(2,i), krt(3,i), krt(4,i);
      end do
   close(1)

   !! SIMULATION INFO
   open(1, file='../../../data/position-distribution/information/resetting_info_pos.txt', status='replace')
      do i=1,5
         write(1,*) infoName(i), info(i);
      end do
   close(1)

   !! RAW MOMENTS
   open(1, file='../../../data/position-distribution/statistical-moments/resetting_m1_pos.txt', status='replace')
      do i=1,momentsTime
         write(1,*) s1(1,i), s1(2,i), s1(3,i), s1(4,i);
      end do
   close(1)
   open(1, file='../../../data/position-distribution/statistical-moments/resetting_m2_pos.txt', status='replace')
      do i=1,momentsTime
         write(1,*) s2(1,i), s2(2,i), s2(3,i), s2(4,i);
      end do
   close(1)
   open(1, file='../../../data/position-distribution/statistical-moments/resetting_m3_pos.txt', status='replace')
      do i=1,momentsTime
         write(1,*) s3(1,i), s3(2,i), s3(3,i), s3(4,i);
      end do
   close(1)
   open(1, file='../../../data/position-distribution/statistical-moments/resetting_m4_pos.txt', status='replace')
      do i=1,momentsTime
         write(1,*) s4(1,i), s4(2,i), s4(3,i), s4(4,i);
      end do
   close(1)

end program
