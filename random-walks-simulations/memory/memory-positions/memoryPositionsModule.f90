module memoryPositionsModule
   use omp_lib
   implicit none

   !! RANDOM WALK PARAMETERS
   integer, parameter :: cases = 4;
   integer, parameter :: realizations = 100000;
   integer, parameter :: maxDistance = 400;
   integer, parameter :: steps = 10000;
   integer, dimension(2*maxDistance+1) :: positions = 0;
   integer, dimension(cases,2*maxDistance+1) :: counts = 0;
   integer :: randomWalk = 0
   integer :: choice
   real :: tempRand

   real, dimension(cases) :: probability = [0.0, 0.05, 0.1, 0.15];
   real :: probabilityTemp;
   integer, dimension(2*maxDistance+1) :: memory

   !! RAW MOMENTS. VARIANCE, SKEWNESS AND KURTOSIS
   integer, parameter :: momentsTime = 1000;
   integer, allocatable :: momentCounts(:,:);
   real, dimension(cases,momentsTime) :: s1, s2, s3, s4
   real, dimension(cases,momentsTime) :: skw, krt, var

   !! INFORMATION STORING
   integer, dimension(5) :: info = [cases,realizations,steps,maxDistance,momentsTime]
   character(len=15),dimension(5) :: infoName = [character(len=15) :: "Cases","Realizations",&
                                                "Steps","MaxDistance","Time"];

   !! INDICES
   integer :: i,ii,iii,j,z,k;

contains

   !!  INITIALIZER OF THE RANDOM SEED
   !!  
   !!  Default subroutine to use CALL RANDOM_NUMBER()
   subroutine init_random_seed()
      ! Default parameters
      integer, allocatable :: seed(:)
      integer :: n, thread_id
      integer :: t

      thread_id = omp_get_thread_num()
      ! Random seed generator
      call random_seed(size = n)
      allocate(seed(n))
      call system_clock(count=t)
      seed = t + (thread_id + 1) * 12345 + [(i,i=1,n)]
      call random_seed(put=seed)
      deallocate(seed)
   end subroutine

end module
