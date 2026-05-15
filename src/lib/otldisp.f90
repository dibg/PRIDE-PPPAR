subroutine otldisp(jd,sod,otlfil,disp)
  implicit none
  integer*4,intent(in) :: jd
  real*8,intent(in) :: sod
  character*(*),intent(in) :: otlfil
  real*8,intent(out) :: disp(3)
  !
  !! function called
  integer*4, external :: get_valid_unit
  real*8, external :: timdif
  !
  !! local
  integer*4 :: lfnotl,stat,i,nrec,il,ir,nlead
  integer*4 :: mjdr
  real*8 :: sodr,enur(3)
  !
  !! in-memory ocean-loading table, loaded once per file (process-local,
  !! not thread-safe). The original routine re-opened the file and scanned
  !! it from the start on every call; here it is read once and a hint
  !! pointer keeps the per-epoch lookup at O(1) amortized.
  integer*4, save :: ntbl = 0
  integer*4, save :: hint = 0
  integer*4, allocatable, save :: mjd_tbl(:)
  real*8,    allocatable, save :: sod_tbl(:)
  real*8,    allocatable, save :: enu_tbl(:,:)
  character*256, save :: cached_otlfil = ''
  !
  !! (re)load the table when the file argument changes
  if (trim(otlfil) .ne. trim(cached_otlfil)) then
    if (allocated(mjd_tbl)) deallocate(mjd_tbl,sod_tbl,enu_tbl)
    lfnotl = get_valid_unit(10)
    !
    !! first pass: count records
    open(lfnotl,file=trim(otlfil),access='stream',form='unformatted',status="old")
    nrec = 0
    do
      read(lfnotl,iostat=stat) mjdr,sodr,enur(1:3)
      if (stat .ne. 0) exit
      nrec = nrec + 1
    end do
    close(lfnotl)
    !
    !! second pass: store records
    allocate(mjd_tbl(nrec),sod_tbl(nrec),enu_tbl(3,nrec))
    open(lfnotl,file=trim(otlfil),access='stream',form='unformatted',status="old")
    do i=1,nrec
      read(lfnotl,iostat=stat) mjd_tbl(i),sod_tbl(i),enu_tbl(1:3,i)
      if (stat .ne. 0) exit
    end do
    close(lfnotl)
    ntbl = nrec
    hint = 0
    cached_otlfil = otlfil
  endif
  !
  if (ntbl .le. 0) then
    disp(1:3) = 0.d0
    return
  endif
  !
  !! nlead = number of leading records at or before the requested epoch.
  !! Records are time-ordered, so timdif is monotonically decreasing; the
  !! hint seed only changes how many steps the walk below takes, not nlead.
  nlead = min(max(hint,0),ntbl)
  do while (nlead .lt. ntbl)
    if (timdif(jd,sod,mjd_tbl(nlead+1),sod_tbl(nlead+1)) .ge. 0.d0) then
      nlead = nlead + 1
    else
      exit
    endif
  end do
  do while (nlead .gt. 0)
    if (timdif(jd,sod,mjd_tbl(nlead),sod_tbl(nlead)) .ge. 0.d0) exit
    nlead = nlead - 1
  end do
  hint = nlead
  !
  !! bracketing indices, reproducing the original linear-scan semantics:
  !! "left" starts at record 1, "right" is the first record past the epoch.
  il = max(nlead,1)
  ir = min(max(nlead,1)+1,ntbl)
  !
  !! return whichever bracketing record is closer in time
  if (abs(timdif(jd,sod,mjd_tbl(il),sod_tbl(il))) .ge. &
      abs(timdif(jd,sod,mjd_tbl(ir),sod_tbl(ir)))) then
    disp(1:3) = enu_tbl(1:3,ir)
  else
    disp(1:3) = enu_tbl(1:3,il)
  endif

end subroutine
