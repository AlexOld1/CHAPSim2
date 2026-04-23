module io_tools_mod
  use print_msg_mod
  use parameters_constant_mod
  use decomp_2d_io
  use udf_type_mod
  use io_files_mod
  implicit none

  !----------------------------------------------------------------------------------------------------------
  ! io parameters
  !----------------------------------------------------------------------------------------------------------
  character(*), parameter :: io_restart = "restart-io"
  character(*), parameter :: io_in2outlet = "outlet2inlet-io"

  public :: initialise_decomp_io
  !public :: generate_file_name
  public :: generate_pathfile_name

  public  :: write_one_3d_array
  public  :: read_one_3d_array

contains

!==========================================================================================================
!==========================================================================================================
  subroutine read_one_3d_array(var, keyword, idom, iter, dtmp)
    use iso_fortran_env, only: int64
    implicit none 
    integer, intent(in) :: idom
    character(*), intent(in) :: keyword
    integer, intent(in) :: iter
    type(DECOMP_INFO), intent(in) :: dtmp
    real(WP), dimension(:, :, :), intent(out) :: var( dtmp%xsz(1), &
                                                      dtmp%xsz(2), &
                                                      dtmp%xsz(3))
    character(64):: data_flname_path
    character(512) :: msg
    integer(int64) :: file_bytes, expected_bytes, elem_bytes

    call generate_pathfile_name(data_flname_path, idom, trim(keyword), dir_data, 'bin', iter)
    if(.not.file_exists(data_flname_path)) &
    call Print_error_msg("The file "//trim(data_flname_path)//" does not exist.")

    inquire(file=trim(data_flname_path), size=file_bytes)
    elem_bytes = int(storage_size(var), int64) / 8_int64
    expected_bytes = int(dtmp%xsz(1), int64) * &
                     int(dtmp%ysz(2), int64) * &
                     int(dtmp%zsz(3), int64) * elem_bytes
    if(file_bytes /= expected_bytes) then
      write(msg, '(A,A,A,I0,A,I0,A,I0,A,I0,A,I0,A)') &
        'Read size mismatch for ', trim(data_flname_path), ': file has ', &
        file_bytes, ' bytes, expected ', expected_bytes, ' bytes for global shape (', &
        dtmp%xsz(1), ',', dtmp%ysz(2), ',', dtmp%zsz(3), &
        '). Regenerate the matching restart/outlet file.'
      call Print_error_msg(trim(msg))
    end if

    if(nrank == 0) call Print_debug_inline_msg("Reading "//trim(data_flname_path))

    call decomp_2d_read_one(IPENCIL(1), var, trim(data_flname_path), &
          opt_decomp=dtmp, &
          opt_reduce_prec=.false.)

    return
  end subroutine
!==========================================================================================================
!==========================================================================================================
  subroutine write_one_3d_array(var, field_name, idom, iter, dtmp, io_mode)
    use typeconvert_mod
    implicit none 
    real(WP), contiguous, intent(in) :: var( :, :, :)
    type(DECOMP_INFO), intent(in) :: dtmp
    character(*), intent(in) :: field_name
    integer, intent(in) :: idom
    integer, intent(in) :: iter
    integer, intent(in) :: io_mode

    character(256):: field_file
    logical :: ex
    character(len=:), allocatable :: bak_file
    logical :: do_write

    call generate_pathfile_name(field_file, idom, trim(field_name), dir_data, 'bin', iter)

    do_write = .true.
    if(io_mode == IO_MODE_SKIP) then
      if (file_exists(trim(field_file))) then
        if (nrank == 0) then
          call Print_warning_msg("File "//trim(field_file)// &
                                " already exists; skip writing "//trim(field_name)// &
                                " at iteration "//trim(int2str(iter)))
        end if
        do_write = .false.
      end if
    else if (io_mode == IO_MODE_RENAME) then
      if (file_exists(trim(field_file))) then
        call rename_existing_file(trim(field_file))
      end if
    else
      ! do nothing, just overwrite if file exists
    end if

    if (do_write) then
      if(nrank == 0) call Print_debug_mid_msg("Writing "//trim(field_file))
      call decomp_2d_write_one(IPENCIL(1), var, trim(field_file), opt_decomp=dtmp)
    end if
    !
    return
  end subroutine
!==========================================================================================================
  subroutine rename_existing_file(file_w_path)
    implicit none
    character(*), intent(in) :: file_w_path
    character(256) :: bak_file
    logical :: ex
    !
    ex = file_exists(trim(file_w_path))
    if(nrank == 0 .and. ex) then
      bak_file = trim(file_w_path)//'.bak'
      if (file_exists(trim(bak_file))) then
        call execute_command_line('rm -f ' // trim(bak_file))
      end if
      call execute_command_line('mv ' // trim(file_w_path) // ' ' // trim(bak_file))
    end if
    return
  end subroutine rename_existing_file
  !==========================================================================================================

  
!==========================================================================================================
  subroutine initialise_decomp_io(dm)
    use udf_type_mod
    use decomp_2d_io
    implicit none 
    type(t_domain), intent(in) :: dm
    
!----------------------------------------------------------------------------------------------------------
! if not #ifdef ADIOS2, do nothing below.
!----------------------------------------------------------------------------------------------------------
    call decomp_2d_io_init()
!----------------------------------------------------------------------------------------------------------
! re-define the grid mesh size, considering the nskip
! based on decomp_info of dppp (default one defined)
!---------------------------------------------------------------------------------------------------------- 
    ! if(dm%visu_nskip(1) > 1 .or. dm%visu_nskip(2) > 1 .or. dm%visu_nskip(3) > 1) then
    !   call init_coarser_mesh_statV(dm%visu_nskip(1), dm%visu_nskip(2), dm%visu_nskip(3), from1=.true.)
    ! end if
    !call init_coarser_mesh_statS(dm%stat_nskip(1), dm%stat_nskip(2), dm%stat_nskip(3), is_start1)

  end subroutine 
!==========================================================================================================
  subroutine generate_pathfile_name(flname_path, dmtag, keyword, path, extension, opt_timetag, opt_flname)
    use typeconvert_mod
    implicit none 
    integer, intent(in)      :: dmtag
    character(*), intent(in) :: keyword
    character(*), intent(in) :: path
    character(*), intent(in) :: extension
    character(*), intent(inout), optional :: opt_flname
    character(*), intent(out) :: flname_path
    integer, intent(in), optional     :: opt_timetag
    character(64) :: flname

    if(present(opt_timetag)) then
      flname = "/domain"//trim(int2str(dmtag))//'_'//trim(keyword)//'_'//trim(int2str(opt_timetag))//"."//trim(extension)
    else 
      flname = "/domain"//trim(int2str(dmtag))//'_'//trim(keyword)//"."//trim(extension)
    end if
    if(present(opt_flname)) then
      opt_flname = flname
    end if
    flname_path = trim(path)//trim(flname)

    return
  end subroutine
!==========================================================================================================
  ! subroutine generate_file_name(flname, dmtag, keyword, extension, timetag)
  !   use typeconvert_mod
  !   implicit none 
  !   integer, intent(in)      :: dmtag
    
  !   character(*), intent(in) :: keyword
  !   character(*), intent(in) :: extension
  !   character(*), intent(out) :: flname
  !   integer, intent(in), optional      :: timetag

  !   if(present(timetag)) then
  !     flname = "domain"//trim(int2str(dmtag))//'_'//trim(keyword)//'_'//trim(int2str(timetag))//"."//trim(extension)
  !   else 
  !     flname = "domain"//trim(int2str(dmtag))//'_'//trim(keyword)//"."//trim(extension)
  !   end if

    
  !   return
  ! end subroutine
!==========================================================================================================
end module
