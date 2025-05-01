module resettingVisitsModule
   use omp_lib
   implicit none

   !! RANDOM WALK PARAMETERS
   integer, parameter :: cases = 4;
   integer, parameter :: realizations = 100000;
   integer, parameter :: maxDistance = 400;
   integer, parameter :: steps = 10000;
   integer, dimension(2*maxDistance+1) :: positions = 0;
   integer, dimension(cases,2*maxDistance+1) :: counts = 0;
   integer, dimension(2*maxDistance+1) :: countsTemp = 0;
   integer, allocatable :: randomWalk(:), randomWalkTemp(:);
   integer :: choice
   real :: tempRand
   integer, allocatable :: seedsArray(:,:);
   integer :: seedSize;

   real, dimension(cases) :: probability = [0.0, 0.05, 0.1, 0.15];
   real :: probabilityTemp;

   !! RAW MOMENTS. VARIANCE, SKEWNESS AND KURTOSIS
   integer, parameter :: momentsTime = 1000;
   real, dimension(cases,momentsTime) :: s1, s2, s3, s4
   real, dimension(cases,momentsTime) :: skw, krt, var

   !! INFORMATION STORING
   integer, dimension(5) :: info = [cases,realizations,steps,maxDistance,momentsTime]
   character(len=15),dimension(5) :: infoName = [character(len=15) :: "Cases","Realizations","Steps",&
                                                "MaxDistance","Time"];

   !! INDICES
   integer :: i,ii,iii,j,k;

contains

   !!  INITIALIZER OF THE RANDOM SEED
   !!  
   !!  Default subroutine to use CALL RANDOM_NUMBER()
   subroutine init_random_seed(seedOutput)
      ! Default parameters
      integer, intent(out) :: seedOutput(:)
      integer :: thread_id,clock
      integer,dimension(8) :: t

      thread_id = omp_get_thread_num()
      call system_clock(count=clock)
      call date_and_time(values=t)
      seedOutput = clock + (thread_id + 1) * 179423311 + [(i*12345,i=1,size(seedOutput))]
   end subroutine

end module
