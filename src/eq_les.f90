module les_mod
    use decomp_2d
    use operations
    use precision_mod
    use udf_type_mod
    implicit none
    private :: init_les, calculate_cell_grad, calculate_stress_tensor, calculate_stress_tensor_square
    private :: calculate_WALE_tensor, calculate_WALE_invariants

    real(WP), allocatable :: dudx(:, :, :, :, :)
    real(WP), allocatable :: S(:, :, :, :, :)
    real(WP), allocatable :: Ssqr(:, :, :)
    real(WP), allocatable :: SD(:, :, :, :, :)
    real(WP), allocatable :: WALE_invariants(:, :, :)
    real(WP), allocatable :: nu_t(:, :, :)

    public :: calculate_les_wale
    public :: calculate_les_smag

    logical :: les_initialized = .false.

contains
    subroutine init_les(fl, dm)
        use parameters_constant_mod
        type(t_flow),   intent(inout) :: fl
        type(t_domain), intent(in)    :: dm

        allocate(dudx(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3), 3, 3))
        allocate(S   (dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3), 3, 3))
        allocate(Ssqr(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3)))
        allocate(SD  (dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3), 3, 3))
        allocate(WALE_invariants(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3)))
        allocate(nu_t           (dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3)))
        allocate(fl%tVisc       (dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3))) ; fl%tVisc = ZERO

        les_initialized = .true.
    end subroutine init_les

    subroutine calculate_cell_grad(fl, dm)
        use parameters_constant_mod
        use boundary_conditions_mod
        use transpose_extended_mod
        type(t_flow),   intent(in) :: fl
        type(t_domain), intent(in) :: dm

        real(WP), dimension( dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3) ) :: accc_xpencil
        real(WP), dimension( dm%dpcc%xsz(1), dm%dpcc%xsz(2), dm%dpcc%xsz(3) ) :: apcc_xpencil
        real(WP), dimension( dm%dppc%xsz(1), dm%dppc%xsz(2), dm%dppc%xsz(3) ) :: appc_xpencil
        real(WP), dimension( dm%dcpc%xsz(1), dm%dcpc%xsz(2), dm%dcpc%xsz(3) ) :: acpc_xpencil
        real(WP), dimension( dm%dccp%xsz(1), dm%dccp%xsz(2), dm%dccp%xsz(3) ) :: accp_xpencil
        real(WP), dimension( dm%dpcp%xsz(1), dm%dpcp%xsz(2), dm%dpcp%xsz(3) ) :: apcp_xpencil
        real(WP), dimension( dm%dpcc%ysz(1), dm%dpcc%ysz(2), dm%dpcc%ysz(3) ) :: apcc_ypencil
        real(WP), dimension( dm%dppc%ysz(1), dm%dppc%ysz(2), dm%dppc%ysz(3) ) :: appc_ypencil
        real(WP), dimension( dm%dccc%ysz(1), dm%dccc%ysz(2), dm%dccc%ysz(3) ) :: accc_ypencil
        real(WP), dimension( dm%dcpc%ysz(1), dm%dcpc%ysz(2), dm%dcpc%ysz(3) ) :: acpc_ypencil
        real(WP), dimension( dm%dccp%ysz(1), dm%dccp%ysz(2), dm%dccp%ysz(3) ) :: accp_ypencil
        real(WP), dimension( dm%dcpp%ysz(1), dm%dcpp%ysz(2), dm%dcpp%ysz(3) ) :: acpp_ypencil
        real(WP), dimension( dm%dpcc%zsz(1), dm%dpcc%zsz(2), dm%dpcc%zsz(3) ) :: apcc_zpencil
        real(WP), dimension( dm%dpcp%zsz(1), dm%dpcp%zsz(2), dm%dpcp%zsz(3) ) :: apcp_zpencil
        real(WP), dimension( dm%dccc%zsz(1), dm%dccc%zsz(2), dm%dccc%zsz(3) ) :: accc_zpencil
        real(WP), dimension( dm%dccp%zsz(1), dm%dccp%zsz(2), dm%dccp%zsz(3) ) :: accp_zpencil
        real(WP), dimension( dm%dcpc%zsz(1), dm%dcpc%zsz(2), dm%dcpc%zsz(3) ) :: acpc_zpencil
        real(WP), dimension( dm%dcpp%zsz(1), dm%dcpp%zsz(2), dm%dcpp%zsz(3) ) :: acpp_zpencil
        real(WP), dimension( dm%dppc%ysz(1), 4, dm%dppc%ysz(3) ) :: fbcy_p4c
        real(WP), dimension( dm%dcpp%ysz(1), 4, dm%dcpp%ysz(3) ) :: fbcy_c4p

        ! du/dx, du/dy, du/dz
        call Get_x_1der_P2C_3D(fl%qx, accc_xpencil, dm, dm%iAccuracy, dm%ibcx_qx, dm%fbcx_qx)
        dudx(:, :, :, 1, 1) = accc_xpencil(:, :, :)
        call transpose_x_to_y(fl%qx, apcc_ypencil, dm%dpcc)
        call Get_y_1der_C2P_3D(apcc_ypencil, appc_ypencil, dm, dm%iAccuracy, dm%ibcy_qx, dm%fbcy_qx)
        fbcy_p4c = MAXP
        if(dm%icase == ICASE_PIPE) then
        call axis_mirror_fbcy(appc_ypencil, IPENCIL(2), fbcy_p4c, dm%knc_sym, dm%dppc, is_odd = .true., &
                                axis_mode = AXIS_RECON_M1, assign_axis_to_var = .true., nr = 0, opt_dz = dm%h(3))
        end if
        call Get_y_midp_P2C_3D(appc_ypencil, apcc_ypencil, dm, dm%iAccuracy, dm%ibcy_qx)
        call transpose_y_to_x(apcc_ypencil, apcc_xpencil, dm%dpcc)
        call Get_x_midp_P2C_3D(apcc_xpencil, accc_xpencil, dm, dm%iAccuracy, dm%ibcx_qx)
        dudx(:, :, :, 1, 2) = accc_xpencil(:, :, :)
        call transpose_to_z_pencil(fl%qx, apcc_zpencil, dm%dpcc, IPENCIL(1))
        call Get_z_1der_C2P_3D(apcc_zpencil, apcp_zpencil, dm, dm%iAccuracy, dm%ibcz_qx, dm%fbcz_qx)
        call Get_z_midp_P2C_3D(apcp_zpencil, apcc_zpencil, dm, dm%iAccuracy, dm%ibcz_qx)
        call transpose_from_z_pencil(apcc_zpencil, apcc_xpencil, dm%dccc, IPENCIL(1))
        call Get_x_midp_P2C_3D(apcc_xpencil, accc_xpencil, dm, dm%iAccuracy, dm%ibcx_qx)
        dudx(:, :, :, 1, 3) = accc_xpencil(:, :, :)
        ! dv/dx, dv/dy, dv/dz
        call Get_x_1der_C2P_3D(fl%qy, appc_xpencil, dm, dm%iAccuracy, dm%ibcx_qy, dm%fbcx_qy)
        call Get_x_midp_P2C_3D(appc_xpencil, acpc_xpencil, dm, dm%iAccuracy, dm%ibcx_qy)
        call transpose_x_to_y(acpc_xpencil, acpc_ypencil, dm%dcpc)
        call Get_y_midp_P2C_3D(acpc_ypencil, accc_ypencil, dm, dm%iAccuracy, dm%ibcy_qy)
        call transpose_y_to_x(accc_ypencil, accc_xpencil, dm%dccc)
        dudx(:, :, :, 2, 1) = accc_xpencil(:, :, :)
        call transpose_x_to_y(fl%qy, acpc_ypencil, dm%dcpc)
        call Get_y_1der_P2C_3D(acpc_ypencil, accc_ypencil, dm, dm%iAccuracy, dm%ibcy_qy, dm%fbcy_qy)
        call transpose_y_to_x(accc_ypencil, accc_xpencil, dm%dccc)
        dudx(:, :, :, 2, 2) = accc_xpencil(:, :, :)
        call transpose_to_z_pencil(fl%qy, acpc_zpencil, dm%dcpc, IPENCIL(1))
        call Get_z_1der_C2P_3D(acpc_zpencil, acpp_zpencil, dm, dm%iAccuracy, dm%ibcz_qy, dm%fbcz_qy)
        call Get_z_midp_P2C_3D(acpp_zpencil, acpc_zpencil, dm, dm%iAccuracy, dm%ibcz_qy)
        call transpose_z_to_y(acpc_zpencil, acpc_ypencil, dm%dcpc)
        call Get_y_midp_P2C_3D(acpc_ypencil, accc_ypencil, dm, dm%iAccuracy, dm%ibcy_qy)
        call transpose_y_to_x(accc_ypencil, accc_xpencil, dm%dccc)
        dudx(:, :, :, 2, 3) = accc_xpencil(:, :, :)
        ! dw/dx, dw/dy, dw/dz
        call Get_x_1der_C2P_3D(fl%qz, apcp_xpencil, dm, dm%iAccuracy, dm%ibcx_qz, dm%fbcx_qz)
        call Get_x_midp_P2C_3D(apcp_xpencil, accp_xpencil, dm, dm%iAccuracy, dm%ibcx_qz)
        call transpose_to_z_pencil(accp_xpencil, accp_zpencil, dm%dccp, IPENCIL(1))
        call Get_z_midp_P2C_3D(accp_zpencil, accc_zpencil, dm, dm%iAccuracy, dm%ibcz_qz)
        call transpose_from_z_pencil(accc_zpencil, accc_xpencil, dm%dccc, IPENCIL(1))
        dudx(:, :, :, 3, 1) = accc_xpencil(:, :, :)
        call transpose_x_to_y(fl%qz, accp_ypencil, dm%dccp)
        call Get_y_1der_C2P_3D(accp_ypencil, acpp_ypencil, dm, dm%iAccuracy, dm%ibcy_qz, dm%fbcy_qz)
        fbcy_c4p = MAXP
        if(dm%icase == ICASE_PIPE) then
        call axis_mirror_fbcy(acpp_ypencil, IPENCIL(2), fbcy_c4p, dm%knc_sym, dm%dcpp, is_odd = .false., &
                                axis_mode = AXIS_RECON_M0_M2, assign_axis_to_var = .true., nr = 0, opt_dz = dm%h(3))
        end if
        call Get_y_midp_P2C_3D(acpp_ypencil, accp_ypencil, dm, dm%iAccuracy, dm%ibcy_qz)
        call transpose_to_z_pencil(accp_ypencil, accp_zpencil, dm%dccp, IPENCIL(2))
        call Get_z_midp_P2C_3D(accp_zpencil, accc_zpencil, dm, dm%iAccuracy, dm%ibcz_qz)
        call transpose_from_z_pencil(accc_zpencil, accc_xpencil, dm%dccc, IPENCIL(1))
        dudx(:, :, :, 3, 2) = accc_xpencil(:, :, :)
        call transpose_to_z_pencil(fl%qz, accp_zpencil, dm%dccp, IPENCIL(1))
        call Get_z_1der_P2C_3D(accp_zpencil, accc_zpencil, dm, dm%iAccuracy, dm%ibcz_qz, dm%fbcz_qz)
        call transpose_from_z_pencil(accc_zpencil, accc_xpencil, dm%dccc, IPENCIL(1))
        dudx(:, :, :, 3, 3) = accc_xpencil(:, :, :)

    end subroutine calculate_cell_grad


    subroutine calculate_stress_tensor()
        integer :: i, j
        do i = 1, 3
            do j = 1, 3
                S(:, :, :, i, j) = 0.5_WP * (dudx(:, :, :, i, j) + dudx(:, :, :, j, i))
            end do
        end do
    end subroutine calculate_stress_tensor


    subroutine calculate_stress_tensor_square()
        Ssqr(:, :, :) = S(:, :, :, 1, 1)**2 + S(:, :, :, 1, 2)**2 + S(:, :, :, 1, 3)**2 &
                        + S(:, :, :, 2, 1)**2 + S(:, :, :, 2, 2)**2 + S(:, :, :, 2, 3)**2 &
                        + S(:, :, :, 3, 1)**2 + S(:, :, :, 3, 2)**2 + S(:, :, :, 3, 3)**2
    end subroutine calculate_stress_tensor_square


    subroutine calculate_wale_tensor()
        integer :: i, j
        real(WP), dimension(size(dudx,1), size(dudx,2), size(dudx,3)) :: trace_dudx2

        trace_dudx2 = dudx(:, :, :, 1, 1)**2 + dudx(:, :, :, 2, 2)**2 + dudx(:, :, :, 3, 3)**2

        do i = 1, 3
            do j = 1, 3
                SD(:, :, :, i, j) = 0.5_WP * (dudx(:, :, :, i, 1) * dudx(:, :, :, 1, j) + &
                                               dudx(:, :, :, i, 2) * dudx(:, :, :, 2, j) + &
                                               dudx(:, :, :, i, 3) * dudx(:, :, :, 3, j) + &
                                               dudx(:, :, :, j, 1) * dudx(:, :, :, 1, i) + &
                                               dudx(:, :, :, j, 2) * dudx(:, :, :, 2, i) + &
                                               dudx(:, :, :, j, 3) * dudx(:, :, :, 3, i))
                if (i == j) then
                    SD(:, :, :, i, j) = SD(:, :, :, i, j) - (1.0_WP/3.0_WP) * trace_dudx2
                end if
            end do
        end do

    end subroutine calculate_wale_tensor


    subroutine calculate_wale_invariants()
        WALE_invariants(:, :, :) = SD(:, :, :, 1, 1)**2 + SD(:, :, :, 1, 2)**2 + SD(:, :, :, 1, 3)**2 &
                                    + SD(:, :, :, 2, 1)**2 + SD(:, :, :, 2, 2)**2 + SD(:, :, :, 2, 3)**2 &
                                    + SD(:, :, :, 3, 1)**2 + SD(:, :, :, 3, 2)**2 + SD(:, :, :, 3, 3)**2
    end subroutine calculate_wale_invariants

    subroutine calculate_eddy_viscosity_wale(dm)
        use parameters_constant_mod
        type(t_domain), intent(in) :: dm
        real(WP) :: Cw, delta

        Cw    = 0.5_WP
        delta = (dm%h(1) * dm%h(2) * dm%h(3))**(1.0_WP/3.0_WP)

        nu_t = (Cw * delta)**2 * (WALE_invariants)**(3.0_WP/2.0_WP) / &
               ((Ssqr)**(3.0_WP/2.0_WP) + (WALE_invariants)**(5.0_WP/4.0_WP))

    end subroutine calculate_eddy_viscosity_wale

    subroutine calculate_les_wale(fl, dm)
        use parameters_constant_mod
        type(t_flow),   intent(inout) :: fl
        type(t_domain), intent(in)    :: dm

        if (.not. les_initialized) call init_les(fl, dm)

        call calculate_cell_grad(fl, dm)
        call calculate_stress_tensor()
        call calculate_stress_tensor_square()
        call calculate_WALE_tensor()
        call calculate_WALE_invariants()
        call calculate_eddy_viscosity_wale(dm)
        
        if (allocated(fl%dDens)) then
            fl%tVisc = nu_t * fl%dDens  ! Density only allocated if running with variable properties
        else
            fl%tVisc = nu_t             ! For isothermal flows
        end if

    end subroutine calculate_les_wale

    subroutine calculate_les_smag(fl, dm)
        use parameters_constant_mod
        type(t_flow),   intent(inout) :: fl
        type(t_domain), intent(in)    :: dm

        ! call calculate_cell_grad(fl, dm)
        ! call calculate_stress_tensor()
        ! call calculate_stress_tensor_square()

        ! mag_S = sqrt(2 * Ssqr)
        ! nu_t_smag = (Cs * delta)**2 * mag_S
        ! visc = visc + nu_t_smag

    end subroutine calculate_les_smag
end module les_mod
