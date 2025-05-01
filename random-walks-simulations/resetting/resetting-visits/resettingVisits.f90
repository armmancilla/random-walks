include "resettingVisitsModule.f90"
program resettingVisits
   !!  --------------------------------------------------------------------
   !!  |  SIMULATION PROGRAM FOR THE SIMPLE RANDOM WALKER FOR THE VISITS  |
   !!  --------------------------------------------------------------------
   !!
   use omp_lib
   use resettingVisitsModule
   implicit none

   !!  VARIABLES TO MEASURE EXECUTION TIME
   integer :: t0, t1, CountRate;

   !!  FIRST CHECKPOINT
   call system_clock(t0);
   
   !! INITIALIZING LISTS AS ZERO
   s1 = 0.; s2 = 0.; s3 = 0.; s4 = 0.;
   skw = 0.; krt = 0.; var = 0.;
   
   !! POSITIONS LIST
   do i=1,2*maxDistance+1
      positions(i) = i-maxDistance-1
   end do

   !! INITIALIZING RANDOM SEEDS
   call random_seed(size=seedSize);
   !$OMP PARALLEL
   !$OMP SINGLE
   allocate(seedsArray(seedSize,omp_get_num_threads()));
   !$OMP END SINGLE
   call init_random_seed(seedsArray(:,omp_get_thread_num()+1));
   !$OMP END PARALLEL

   !! INITIALIZING RANDOM WALK LISTS
   allocate(randomWalk(realizations));
   allocate(randomWalkTemp(realizations));
   randomWalk = 0; randomWalkTemp = 0;
   
   !! RANDOM WALK SIMULATION
   do i=1,cases
      probabilityTemp = probability(i);
      randomWalk = 0;
      
      !! STEPS LOOP
      do ii=1,steps

         !$OMP PARALLEL PRIVATE(tempRand,iii,k,countsTemp)
         if (ii==1) then
            call random_seed(put=seedsArray(:,omp_get_thread_num()+1));
         end if
         
         countsTemp = 0;
         !$OMP DO
         do iii=1,realizations
            k = randomWalk(iii) + maxDistance + 1;

            ! Counting the visits
            if (abs(randomWalk(iii))<=maxDistance) then
               !$OMP ATOMIC
               counts(i,k) = counts(i,k) + 1;
            end if

            !! STEP DYNAMICS
            call random_number(tempRand);
            if (tempRand <= probabilityTemp) then
               randomWalk(iii) = 0;
            else
               call random_number(tempRand);
               randomWalk(iii) = randomWalk(iii) + 2*nint(tempRand)-1;
            end if

            ! print*,"NT: ",omp_get_thread_num()," II: ",ii," III: ",iii," RN: ",tempRand," RW: ",randomWalk(iii), &
            !       " C: ",(countsTemp(k));

         end do
         !$OMP END DO
         !$OMP END PARALLEL

         !! MOMENTS CALCULATION
         if (ii<=momentsTime) then

            do j=1,2*maxDistance+1
            
               s1(i,ii) = s1(i,ii) + 1.0*counts(i,j)*positions(j)/(realizations*ii)
               s2(i,ii) = s2(i,ii) + 1.0*counts(i,j)*(positions(j)**2)/(realizations*ii)
               s3(i,ii) = s3(i,ii) + 1.0*counts(i,j)*(positions(j)**3)/(realizations*ii)
               s4(i,ii) = s4(i,ii) + 1.0*counts(i,j)*(positions(j)**4)/(realizations*ii)
            
            end do
            
            !! VARIANCE, SKEWNESS AND KURTOSIS
            var(i,ii) = ( s2(i,ii) - s1(i,ii)**2. )
            skw(i,ii) = ( s3(i,ii) - 3.*s1(i,ii)*( var(i,ii) )&
                        - s1(i,ii)**3. )/( var(i,ii)**(3./2.) )
            krt(i,ii) = ( s4(i,ii) - 4.*s1(i,ii)*s3(i,ii) &
                        + 6.*(s1(i,ii)**2.)*s2(i,ii) - 3.*(s1(i,ii)**4.))/&
                        (var(i,ii)**2.)

         end if
      end do
   end do

   deallocate(randomWalk); deallocate(randomWalkTemp);
   deallocate(seedsArray);

   !!  SECOND CHECKPOINT
   call system_clock(t1,CountRate);
   print*,"Execution Time: ",int((t1-t0)*1000./CountRate)/1000.;

   !! DATA STORING SECTION
   !! --------------------

   !! SIMULATION INFO
   open(1, file='../../../data/visits-distribution/information/resetting_info_vis.txt', status='replace')
      do i=1,5
         write(1,*) infoName(i), info(i);
      end do
   close(1)

   !! DISTRIBUTION
   open(1, file='../../../data/visits-distribution/distributions/resetting_distribution_vis.txt', status='replace')
      do i=1,2*maxDistance+1
            write(1,*) counts(1,i), counts(2,i), counts(3,i), counts(4,i);
      end do
   close(1)

   !! CASES
   open(1, file='../../../data/visits-distribution/information/resetting_cases_vis.txt', status='replace')
      do i=1,cases
         write(1,*) probability(i)
      end do
   close(1)

   !! VARIANCE
   open(1, file='../../../data/visits-distribution/statistics/resetting_variance_vis.txt', status='replace')
      do i=1,momentsTime
         write(1,*) var(1,i), var(2,i), var(3,i), var(4,i);
      end do
   close(1)

   !! SKEWNESS
   open(1, file='../../../data/visits-distribution/statistics/resetting_skewness_vis.txt', status='replace')
      do i=1,momentsTime
         write(1,*) skw(1,i), skw(2,i), skw(3,i), skw(4,i);
      end do
   close(1)

   !! KURTOSIS
   open(1, file='../../../data/visits-distribution/statistics/resetting_kurtosis_vis.txt', status='replace')
      do i=1,momentsTime
         write(1,*) krt(1,i), krt(2,i), krt(3,i), krt(4,i);
      end do
   close(1)

   !! RAW MOMENTS
   open(1, file='../../../data/visits-distribution/statistical-moments/resetting_m1_vis.txt', status='replace')
      do i=1,momentsTime
         write(1,*) s1(1,i), s1(2,i), s1(3,i), s1(4,i);
      end do
   close(1)
   open(1, file='../../../data/visits-distribution/statistical-moments/resetting_m2_vis.txt', status='replace')
      do i=1,momentsTime
         write(1,*) s2(1,i), s2(2,i), s2(3,i), s2(4,i);
      end do
   close(1)
   open(1, file='../../../data/visits-distribution/statistical-moments/resetting_m3_vis.txt', status='replace')
      do i=1,momentsTime
         write(1,*) s3(1,i), s3(2,i), s3(3,i), s3(4,i);
      end do
   close(1)
   open(1, file='../../../data/visits-distribution/statistical-moments/resetting_m4_vis.txt', status='replace')
      do i=1,momentsTime
         write(1,*) s4(1,i), s4(2,i), s4(3,i), s4(4,i);
      end do
   close(1)

end program
