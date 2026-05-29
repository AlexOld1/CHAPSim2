module les_mod
  use decomp_2d
  use operations
  use precision_mod, only: WP
  use udf_type_mod, only: t_domain, t_flow
  implicit none
  !
  private :: calculate_velocity_gradient, calculate_strain_rate_tensor
  private :: calculate_strain_rate_magnitude_square, calculate_wale_tensor
  private :: calculate_wale_tensor_magnitude_square, calculate_eddy_viscosity_wale
  !
  public :: calculate_les_wale ! ref; https://www.cfd-online.com/Wiki/Wall-adapting_local_eddy-viscosity_(WALE)_model
  !public :: calculate_les_smag

contains
!==========================================================================================================
  !> Calculate velocity-gradient tensor du_i/dx_j at cell centres.
  !> - fl: flow variables containing staggered velocity components.
  !> - dm: domain/decomposition and boundary-condition metadata.
  !> - velocity_gradient: output tensor with component order (i velocity, j derivative).
  subroutine calculate_velocity_gradient(fl, dm, velocity_gradient)
    use boundary_conditions_mod, only: AXIS_RECON_M0_M2, AXIS_RECON_M1, axis_mirror_fbcy
    use parameters_constant_mod, only: ICASE_PIPE, IPENCIL, MAXP
    use transpose_extended_mod, only: transpose_from_z_pencil, transpose_to_z_pencil
    implicit none
    type(t_flow),   intent(in) :: fl
    type(t_domain), intent(in) :: dm
    real(WP),       intent(out) :: velocity_gradient(:, :, :, :, :)

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
    velocity_gradient(:, :, :, 1, 1) = accc_xpencil(:, :, :)
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
    velocity_gradient(:, :, :, 1, 2) = accc_xpencil(:, :, :)
    call transpose_to_z_pencil(fl%qx, apcc_zpencil, dm%dpcc, IPENCIL(1))
    call Get_z_1der_C2P_3D(apcc_zpencil, apcp_zpencil, dm, dm%iAccuracy, dm%ibcz_qx, dm%fbcz_qx)
    call Get_z_midp_P2C_3D(apcp_zpencil, apcc_zpencil, dm, dm%iAccuracy, dm%ibcz_qx)
    call transpose_from_z_pencil(apcc_zpencil, apcc_xpencil, dm%dccc, IPENCIL(1))
    call Get_x_midp_P2C_3D(apcc_xpencil, accc_xpencil, dm, dm%iAccuracy, dm%ibcx_qx)
    velocity_gradient(:, :, :, 1, 3) = accc_xpencil(:, :, :)
    ! dv/dx, dv/dy, dv/dz
    call Get_x_1der_C2P_3D(fl%qy, appc_xpencil, dm, dm%iAccuracy, dm%ibcx_qy, dm%fbcx_qy)
    call Get_x_midp_P2C_3D(appc_xpencil, acpc_xpencil, dm, dm%iAccuracy, dm%ibcx_qy)
    call transpose_x_to_y(acpc_xpencil, acpc_ypencil, dm%dcpc)
    call Get_y_midp_P2C_3D(acpc_ypencil, accc_ypencil, dm, dm%iAccuracy, dm%ibcy_qy)
    call transpose_y_to_x(accc_ypencil, accc_xpencil, dm%dccc)
    velocity_gradient(:, :, :, 2, 1) = accc_xpencil(:, :, :)
    call transpose_x_to_y(fl%qy, acpc_ypencil, dm%dcpc)
    call Get_y_1der_P2C_3D(acpc_ypencil, accc_ypencil, dm, dm%iAccuracy, dm%ibcy_qy, dm%fbcy_qy)
    call transpose_y_to_x(accc_ypencil, accc_xpencil, dm%dccc)
    velocity_gradient(:, :, :, 2, 2) = accc_xpencil(:, :, :)
    call transpose_to_z_pencil(fl%qy, acpc_zpencil, dm%dcpc, IPENCIL(1))
    call Get_z_1der_C2P_3D(acpc_zpencil, acpp_zpencil, dm, dm%iAccuracy, dm%ibcz_qy, dm%fbcz_qy)
    call Get_z_midp_P2C_3D(acpp_zpencil, acpc_zpencil, dm, dm%iAccuracy, dm%ibcz_qy)
    call transpose_z_to_y(acpc_zpencil, acpc_ypencil, dm%dcpc)
    call Get_y_midp_P2C_3D(acpc_ypencil, accc_ypencil, dm, dm%iAccuracy, dm%ibcy_qy)
    call transpose_y_to_x(accc_ypencil, accc_xpencil, dm%dccc)
    velocity_gradient(:, :, :, 2, 3) = accc_xpencil(:, :, :)
    ! dw/dx, dw/dy, dw/dz
    call Get_x_1der_C2P_3D(fl%qz, apcp_xpencil, dm, dm%iAccuracy, dm%ibcx_qz, dm%fbcx_qz)
    call Get_x_midp_P2C_3D(apcp_xpencil, accp_xpencil, dm, dm%iAccuracy, dm%ibcx_qz)
    call transpose_to_z_pencil(accp_xpencil, accp_zpencil, dm%dccp, IPENCIL(1))
    call Get_z_midp_P2C_3D(accp_zpencil, accc_zpencil, dm, dm%iAccuracy, dm%ibcz_qz)
    call transpose_from_z_pencil(accc_zpencil, accc_xpencil, dm%dccc, IPENCIL(1))
    velocity_gradient(:, :, :, 3, 1) = accc_xpencil(:, :, :)
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
    velocity_gradient(:, :, :, 3, 2) = accc_xpencil(:, :, :)
    call transpose_to_z_pencil(fl%qz, accp_zpencil, dm%dccp, IPENCIL(1))
    call Get_z_1der_P2C_3D(accp_zpencil, accc_zpencil, dm, dm%iAccuracy, dm%ibcz_qz, dm%fbcz_qz)
    call transpose_from_z_pencil(accc_zpencil, accc_xpencil, dm%dccc, IPENCIL(1))
    velocity_gradient(:, :, :, 3, 3) = accc_xpencil(:, :, :)
    return
  end subroutine calculate_velocity_gradient
!==========================================================================================================
  !> Form the resolved strain-rate tensor S_ij = 0.5*(du_i/dx_j + du_j/dx_i).
  !> - velocity_gradient: cell-centred velocity-gradient tensor.
  !> - strain_rate: output symmetric strain-rate tensor.
  subroutine calculate_strain_rate_tensor(velocity_gradient, strain_rate)
    implicit none
    real(WP), intent(in)  :: velocity_gradient(:, :, :, :, :)
    real(WP), intent(out) :: strain_rate(:, :, :, :, :)
    integer :: i, j
    do i = 1, 3
        do j = 1, 3
            strain_rate(:, :, :, i, j) = 0.5_WP * (velocity_gradient(:, :, :, i, j) + &
                                                   velocity_gradient(:, :, :, j, i))
        end do
    end do
    return
  end subroutine calculate_strain_rate_tensor

!==========================================================================================================
  !> Compute S_ij S_ij at each cell centre.
  !> - strain_rate: resolved strain-rate tensor.
  !> - strain_rate_mag2: output tensor contraction S_ij S_ij.
  subroutine calculate_strain_rate_magnitude_square(strain_rate, strain_rate_mag2)
    implicit none
    real(WP), intent(in)  :: strain_rate(:, :, :, :, :)
    real(WP), intent(out) :: strain_rate_mag2(:, :, :)

    strain_rate_mag2(:, :, :) = strain_rate(:, :, :, 1, 1)**2 + strain_rate(:, :, :, 1, 2)**2 + &
                                strain_rate(:, :, :, 1, 3)**2 + strain_rate(:, :, :, 2, 1)**2 + &
                                strain_rate(:, :, :, 2, 2)**2 + strain_rate(:, :, :, 2, 3)**2 + &
                                strain_rate(:, :, :, 3, 1)**2 + strain_rate(:, :, :, 3, 2)**2 + &
                                strain_rate(:, :, :, 3, 3)**2
    return
  end subroutine calculate_strain_rate_magnitude_square

!==========================================================================================================
  !> Compute the traceless symmetric square-gradient tensor used by the WALE model.
  !> - velocity_gradient: cell-centred velocity-gradient tensor.
  !> - wale_tensor: output S^d_ij tensor.
  subroutine calculate_wale_tensor(velocity_gradient, wale_tensor)
      implicit none
      real(WP), intent(in)  :: velocity_gradient(:, :, :, :, :)
      real(WP), intent(out) :: wale_tensor(:, :, :, :, :)
      integer :: i, j
      real(WP), dimension(size(velocity_gradient,1), size(velocity_gradient,2), size(velocity_gradient,3)) :: &
          trace_gradient_square

      trace_gradient_square = velocity_gradient(:, :, :, 1, 1) * velocity_gradient(:, :, :, 1, 1) + &
                              velocity_gradient(:, :, :, 1, 2) * velocity_gradient(:, :, :, 2, 1) + &
                              velocity_gradient(:, :, :, 1, 3) * velocity_gradient(:, :, :, 3, 1) + &
                              velocity_gradient(:, :, :, 2, 1) * velocity_gradient(:, :, :, 1, 2) + &
                              velocity_gradient(:, :, :, 2, 2) * velocity_gradient(:, :, :, 2, 2) + &
                              velocity_gradient(:, :, :, 2, 3) * velocity_gradient(:, :, :, 3, 2) + &
                              velocity_gradient(:, :, :, 3, 1) * velocity_gradient(:, :, :, 1, 3) + &
                              velocity_gradient(:, :, :, 3, 2) * velocity_gradient(:, :, :, 2, 3) + &
                              velocity_gradient(:, :, :, 3, 3) * velocity_gradient(:, :, :, 3, 3)

      do i = 1, 3
          do j = 1, 3
              wale_tensor(:, :, :, i, j) = 0.5_WP * &
                  (velocity_gradient(:, :, :, i, 1) * velocity_gradient(:, :, :, 1, j) + &
                    velocity_gradient(:, :, :, i, 2) * velocity_gradient(:, :, :, 2, j) + &
                    velocity_gradient(:, :, :, i, 3) * velocity_gradient(:, :, :, 3, j) + &
                    velocity_gradient(:, :, :, j, 1) * velocity_gradient(:, :, :, 1, i) + &
                    velocity_gradient(:, :, :, j, 2) * velocity_gradient(:, :, :, 2, i) + &
                    velocity_gradient(:, :, :, j, 3) * velocity_gradient(:, :, :, 3, i))
              if (i == j) then
                  wale_tensor(:, :, :, i, j) = wale_tensor(:, :, :, i, j) - &
                                                (1.0_WP/3.0_WP) * trace_gradient_square
              end if
          end do
      end do
    return
  end subroutine calculate_wale_tensor
!==========================================================================================================
  !> Compute S^d_ij S^d_ij at each cell centre.
  !> - wale_tensor: WALE traceless symmetric square-gradient tensor.
  !> - wale_tensor_mag2: output tensor contraction S^d_ij S^d_ij.
  subroutine calculate_wale_tensor_magnitude_square(wale_tensor, wale_tensor_mag2)
    implicit none
    real(WP), intent(in)  :: wale_tensor(:, :, :, :, :)
    real(WP), intent(out) :: wale_tensor_mag2(:, :, :)

    wale_tensor_mag2(:, :, :) = wale_tensor(:, :, :, 1, 1)**2 + wale_tensor(:, :, :, 1, 2)**2 + &
                                wale_tensor(:, :, :, 1, 3)**2 + wale_tensor(:, :, :, 2, 1)**2 + &
                                wale_tensor(:, :, :, 2, 2)**2 + wale_tensor(:, :, :, 2, 3)**2 + &
                                wale_tensor(:, :, :, 3, 1)**2 + wale_tensor(:, :, :, 3, 2)**2 + &
                                wale_tensor(:, :, :, 3, 3)**2
    return
  end subroutine calculate_wale_tensor_magnitude_square
!==========================================================================================================
  !> Evaluate the WALE kinematic eddy viscosity from strain and WALE invariants.
  !> - dm: domain descriptor used for grid spacing.
  !> - strain_rate_mag2: S_ij S_ij field.
  !> - wale_tensor_mag2: S^d_ij S^d_ij field.
  !> - eddy_visc_kinematic: output nondimensional kinematic eddy viscosity.
  subroutine calculate_eddy_viscosity_wale(dm, strain_rate_mag2, wale_tensor_mag2, eddy_visc_kinematic)
    use parameters_constant_mod, only: ZERO
    implicit none
    type(t_domain), intent(in) :: dm
    real(WP),       intent(in)  :: strain_rate_mag2(:, :, :)
    real(WP),       intent(in)  :: wale_tensor_mag2(:, :, :)
    real(WP),       intent(out) :: eddy_visc_kinematic(:, :, :)
    real(WP) :: Cw, delta
    real(WP), dimension(size(eddy_visc_kinematic,1), size(eddy_visc_kinematic,2), size(eddy_visc_kinematic,3)) :: denominator

    Cw    = 0.5_WP
    delta = (dm%h(1) * dm%h(2) * dm%h(3))**(1.0_WP/3.0_WP)

    denominator = (strain_rate_mag2)**(5.0_WP/2.0_WP) + (wale_tensor_mag2)**(5.0_WP/4.0_WP)
    where (denominator > ZERO)
        eddy_visc_kinematic = (Cw * delta)**2 * (wale_tensor_mag2)**(3.0_WP/2.0_WP) / denominator
    elsewhere
        eddy_visc_kinematic = ZERO
    end where
    return
  end subroutine calculate_eddy_viscosity_wale
!==========================================================================================================
  !> Calculate turbulent dynamic viscosity using the WALE LES model.
  !> - fl: flow state; fl%tVisc is overwritten with turbulent dynamic viscosity.
  !> - dm: domain/decomposition and boundary-condition metadata.
  subroutine calculate_les_wale(fl, dm)
    implicit none
    type(t_flow),   intent(inout) :: fl
    type(t_domain), intent(in)    :: dm
    real(WP), dimension(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3), 3, 3) :: velocity_gradient
    real(WP), dimension(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3), 3, 3) :: strain_rate
    real(WP), dimension(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3), 3, 3) :: wale_tensor
    real(WP), dimension(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3)) :: strain_rate_mag2
    real(WP), dimension(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3)) :: wale_tensor_mag2
    real(WP), dimension(dm%dccc%xsz(1), dm%dccc%xsz(2), dm%dccc%xsz(3)) :: eddy_visc_kinematic

    call calculate_velocity_gradient(fl, dm, velocity_gradient)
    call calculate_strain_rate_tensor(velocity_gradient, strain_rate)
    call calculate_strain_rate_magnitude_square(strain_rate, strain_rate_mag2)
    call calculate_wale_tensor(velocity_gradient, wale_tensor)
    call calculate_wale_tensor_magnitude_square(wale_tensor, wale_tensor_mag2)
    call calculate_eddy_viscosity_wale(dm, strain_rate_mag2, wale_tensor_mag2, eddy_visc_kinematic)

    if (allocated(fl%dDens)) then
        fl%tVisc = fl%dDens * eddy_visc_kinematic
    else
        fl%tVisc = eddy_visc_kinematic
    end if
    return
  end subroutine calculate_les_wale

end module les_mod
