!**********************************************************************************************************************!
!                                        Modification by Mian Muhammad Faheem.                                         !
!**********************************************************************************************************************!
module smss_module

  use kinds_f90
  implicit none

  public
  save

  !Logical file units to read/print SMSS ensemble information.
  integer, parameter :: lusmss = 50

  !Vectors for holding indices for SMSS ensemble.
  integer              :: qmes
  integer, allocatable :: qmid(:), qmit(:)

  interface ReadEnsembleInput
    module procedure ReadEnsembleInput_void
  end interface ReadEnsembleInput

  interface WriteEnsembleForces
    module procedure WriteEnsembleForces_void
  end interface WriteEnsembleForces

  contains

  !******************************************************************************************************************!
  !                                           INTERFACE READENSEMBLEINPUT                                            !
  !******************************************************************************************************************!
  subroutine ReadEnsembleInput_void ()
    implicit none

    character(len=256) :: buff
    integer            :: ii, jj, code

    !Open MMS_ENSEMBLE file and locate $MMS_SHIFTS section.
    open (unit=lusmss, file='MMS_ENSEMBLE', status='OLD', action='READ', position='REWIND', iostat=code)
      if (code>0) call error (9001)
    jj = 0
    do while (.true.)
      read (lusmss,'(A)',end=100,iostat=code) buff
        if (code>0) call error (9002)
      jj = jj + 1
      buff = adjustl(buff)
      if (buff(1:11)=='$MMS_SHIFTS') exit
    end do
    100 continue
    if (jj==0) call error (9021)

    !Find the size of SMSS ensemble.
    qmes = 0
    do while (.true.)
      read (lusmss,'(A)',end=200,iostat=code) buff
        if (code>0) call error (9002)
      buff = adjustl(buff)
      if (buff(1:1)=='$') exit
      qmes = qmes + 1
    end do
    200 continue
    if (qmes==0) call error (9022)

    !Allocate memory for respective arrays.
    if (.not. allocated(qmid)) allocate(qmid(1:qmes))
    if (.not. allocated(qmit)) allocate(qmit(1:qmes))

    !Rewind the file and skip to $MMS_SHIFTS section.
    rewind (lusmss)
    do ii = 1, jj
      read (lusmss,'(A)',iostat=code) buff
        if (code>0) call error (9002)
    end do

    !Now fill data arrays.
    do ii = 1, qmes
      read (lusmss,*,iostat=code) qmit(ii), qmid(ii)
        if (code>0) call error (9002)
    end do

    close (lusmss)

  end subroutine ReadEnsembleInput_void
  !******************************************************************************************************************!

  !******************************************************************************************************************!
  !                                          INTERFACE WRITEENSEMBLEFORCES                                           !
  !******************************************************************************************************************!
  subroutine WriteEnsembleForces_void ()
    use kinds_f90
    use config_module, only: natms, ltg, fxx, fyy, fzz
    use comms_module,  only: idnode, mxnode, gsync
    implicit none

    character(len=256) :: file, cmdl
    integer            :: ii, jj, code
    real(kind=wp)      :: qmfx, qmfy, qmfz

    !Validate SMSS ensemble.
    if ((qmes<1) .or. (qmes>size(ltg))) call error (9024)
    if ((any(qmid<1)) .or. (any(qmid>size(ltg)))) call error (9023)

    !Create a unique file name for each node.
    write (file,'(I6)') idnode
    file = 'smss-tmp-forces-' // trim(adjustl(file))

    !Create/open temporary files.
    open (unit=lusmss, file=trim(file), status='REPLACE', action='WRITE', position='REWIND', iostat=code)
      if (code>0) call error (9011)

    !Identify SMSS ensemble atoms and write indices and forces (in atomic units).
    do ii = 1, qmes
      do jj = 1, natms
        if (qmid(ii)==ltg(jj)) then
          qmfx = fxx(jj) * 2.01553014e-6_wp
          qmfy = fyy(jj) * 2.01553014e-6_wp
          qmfz = fzz(jj) * 2.01553014e-6_wp
          write (lusmss,'(2I10,3F20.14)',iostat=code) qmit(ii), qmid(ii), qmfx, qmfy, qmfz
            if (code>0) call error (9012)
        end if
      end do
    end do

    close (lusmss)
    if (mxnode>1) call gsync()

    !Combine all temporary files.
    if (idnode==0) then
      do ii = 0, mxnode-1
        write (file,'(I6)') ii
        file = 'smss-tmp-forces-' // trim(adjustl(file))
        cmdl = 'cat ' // trim(file) // ' >> MMS_REPLAY'
        call system (trim(cmdl))
      end do
      cmdl = 'rm smss-tmp-forces-*'
      call system (trim(cmdl))
    end if

  end subroutine WriteEnsembleForces_void
  !******************************************************************************************************************!

end module smss_module
!**********************************************************************************************************************!
!                                                 End of Modification.                                                 !
!**********************************************************************************************************************!
