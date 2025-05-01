module simplePositionsModule
   use omp_lib
   implicit none

   !! RANDOM WALK PARAMETERS
   integer, parameter :: cases = 3;
   integer, parameter :: realizations = 100000;
   integer, parameter :: maxDistance = 400;
   integer, dimension(cases) :: steps = [100, 1000, 10000];
   integer, dimension(2*maxDistance+1) :: positions = 0;
   integer, dimension(cases,2*maxDistance+1) :: counts = 0;
   integer :: randomWalk = 0
   integer :: choice
   real :: tempRand

   !! RAW MOMENTS. VARIANCE, SKEWNESS AND KURTOSIS
   integer, parameter :: momentsTime = 100
   integer, allocatable :: momentCounts(:,:);
   real, dimension(cases,momentsTime) :: s1, s2, s3, s4
   real, dimension(cases,momentsTime) :: skw, krt, var

   !! INFORMATION STORING
   integer, dimension(4) :: info = [cases,realizations,maxDistance,momentsTime]
   character(len=15),dimension(4) :: infoName = [character(len=15) :: "Cases","Realizations","MaxDistance","Time"];

   !! INDICES
   integer :: i,ii,iii,j;

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
