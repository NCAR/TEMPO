! TEMPO neural network module designed for cloud droplet number concentration prediction
! A big thanks to David John Gagne (NCAR) for very useful discussions
! Please see https://github.com/NCAR/mlmicrophysics for a more detailed way to connect
! ML models in python to Fortran microphysics code

module module_mp_tempo_ml

  implicit none
  private
  public :: MLdata, tempo_save_or_read_ml_data, predict_number, predict_number_sub

  type MLdata
     integer :: input_size
     integer :: output_size
     integer :: node_size
     real, allocatable, dimension(:) :: transform_mean
     real, allocatable, dimension(:) :: transform_var
     real, allocatable, dimension(:,:) :: weights00
     real, allocatable, dimension(:) :: bias00
     real, allocatable, dimension(:,:) :: weights01
     real, allocatable, dimension(:) :: bias01
  end type MLdata

contains

!----------------------------------------------------------------------------------
  ! Initializes or returns neural network information.
  ! Saves nueral network data information if a neural network is not
  ! initialized. If output data is requested, the initialized and saved
  ! neural network will be used.
  subroutine tempo_save_or_read_ml_data(ml_data_in, ml_data_out)

    logical, save :: not_initialized = .true.
    type(MLdata), dimension(2), intent(in), optional :: ml_data_in
    type(MLdata), dimension(2), intent(out), optional :: ml_data_out
    type(MLdata), dimension(2), save :: tempo_ml_data_save

    if (not_initialized) then
       if (present(ml_data_in)) tempo_ml_data_save = ml_data_in
       not_initialized = .false.
    endif

    if (present(ml_data_out)) then
       ml_data_out = tempo_ml_data_save
    endif

  end subroutine tempo_save_or_read_ml_data

!------------------------------------------------------------------------------------------------------
  ! Predicts number concentration
  real function predict_number(qc, qr, qi, qs, pres, temp, w, predict_nc)

    real, intent(in) :: qc, qr, qi, qs, pres, temp, w
    logical, intent(in) :: predict_nc
    type(MLdata), dimension(2) :: get_ml_data
    type(MLdata) :: ml_data
    integer, parameter :: input_rows = 1
    real, allocatable, dimension(:,:) :: input, input_transformed
    real, allocatable, dimension(:,:) :: output00, output00Activ
    real, allocatable, dimension(:,:) :: output01, output01Activ
    real, parameter :: logMin = -6.0 ! R2
    real, parameter :: logMax = 9.3010299957 ! 2000 cm^-3
    double precision :: predictExp
    
    ! Get neural network data
    call tempo_save_or_read_ml_data(ml_data_out=get_ml_data)

    if (predict_nc) then
       ml_data = get_ml_data(1)
    else
       ml_data = get_ml_data(2)
    endif
    
    ! Allocate arrays for input data and transformation
    if (.not. allocated(input)) then
       allocate(input(input_rows, ml_data%input_size))
    endif
    if (.not. allocated(input_transformed)) then
       allocate(input_transformed(input_rows, ml_data%input_size))
    endif
    
    ! Collect input data
    input(1,1) = qc
    input(1,2) = qr
    input(1,3) = qi
    input(1,4) = qs
    input(1,5) = pres
    input(1,6) = temp
    input(1,7) = w

    ! Transform input data
    call standard_scaler_transform(mean=ml_data%transform_mean, var=ml_data%transform_var, &
         raw_data=input, transformed_data=input_transformed)

    ! Allocate arrays
    if (.not. allocated(output00)) then
       allocate(output00(ml_data%node_size, ml_data%output_size))
    endif
    if (.not. allocated(output00Activ)) then
       allocate(output00Activ(ml_data%node_size, ml_data%output_size))
    endif
    if (.not. allocated(output01)) then
       allocate(output01(ml_data%output_size, ml_data%output_size))
    endif
    if (.not. allocated(output01Activ)) then
       allocate(output01Activ(ml_data%output_size, ml_data%output_size))
    endif

    ! Reconstruct neural network
    ! First layer
    output00 = matmul(transpose(ml_data%weights00), transpose(input_transformed)) + &
         reshape(ml_data%bias00, (/size(ml_data%bias00), ml_data%output_size/))
    call relu_activation(input=output00, output=output00Activ)

    ! Second layer
    output01 = matmul(transpose(ml_data%weights01), output00Activ) + &
         reshape(ml_data%bias01, (/size(ml_data%bias01), ml_data%output_size/))
    call relu_activation(input=output01, output=output01Activ)

    ! Prediction
    predictExp = min(logMax, max(logMin, output01Activ(1,1)))
    predict_number = 10.**predictExp
    
  end function predict_number

!------------------------------------------------------------------------------------------------------
  ! Predicts number concentration
  subroutine predict_number_sub(kts, kte, qc, qr, qi, qs, pres, temp, w, predict_number, predict_nc)

    integer, intent(in) :: kts, kte
    real, dimension(kts:kte), intent(in) :: qc, qr, qi, qs, pres, temp, w
    double precision, dimension(kts:kte), intent(inout) :: predict_number
    logical, intent(in) :: predict_nc
    type(MLdata), dimension(2) :: get_ml_data
    type(MLdata) :: ml_data
    integer, parameter :: input_rows = 1
    real, allocatable, dimension(:,:) :: input, input_transformed
    real, allocatable, dimension(:,:) :: output00, output00Activ, reshaped_bias00
    real, allocatable, dimension(:,:) :: output01, output01Activ, reshaped_bias01
    real, parameter :: logMin = -6.0 ! R2
    real, parameter :: logMax = 9.3010299957 ! 2000 cm^-3
    double precision :: predictExp, bias_corr
    integer :: k

    ! Get neural network data
    call tempo_save_or_read_ml_data(ml_data_out=get_ml_data)

    if (predict_nc) then
       ml_data = get_ml_data(1)
    else
       ml_data = get_ml_data(2)
    endif

    ! Allocate arrays for input data and transformation
    if (.not. allocated(input)) then
       allocate(input(ml_data%input_size, kte))
    endif
    if (.not. allocated(input_transformed)) then
       allocate(input_transformed(ml_data%input_size, kte))
    endif

    ! Collect input data
    input(1,:) = qc
    input(2,:) = qr
    input(3,:) = qi
    input(4,:) = qs
    input(5,:) = pres
    input(6,:) = temp
    input(7,:) = w

    ! Transform input data
    call standard_scaler_transform(mean=ml_data%transform_mean, var=ml_data%transform_var, &
         raw_data=input, transformed_data=input_transformed)

    ! Allocate arrays
    if (.not. allocated(output00)) then
       allocate(output00(ml_data%node_size, kte))
    endif
    if (.not. allocated(output00Activ)) then
       allocate(output00Activ(ml_data%node_size, kte))
    endif
    if (.not. allocated(reshaped_bias00)) then
       allocate(reshaped_bias00(ml_data%node_size, kte))
    endif

    if (.not. allocated(output01)) then
       allocate(output01(ml_data%output_size, kte))
    endif
    if (.not. allocated(output01Activ)) then
       allocate(output01Activ(ml_data%output_size, kte))
    endif
    if (.not. allocated(reshaped_bias01)) then
       allocate(reshaped_bias01(ml_data%output_size, kte))
    endif

    do k = kts, kte
       reshaped_bias00(:,k) = ml_data%bias00
       reshaped_bias01(1,k) = ml_data%bias01(1)
    enddo

    ! Reconstruct neural network
    ! First layer
    output00 = matmul(ml_data%weights00, input_transformed) + reshaped_bias00
    call relu_activation(input=output00, output=output00Activ)

    ! Second layer
    output01 = matmul(ml_data%weights01, output00Activ) + reshaped_bias01
    call relu_activation(input=output01, output=output01Activ)

    do k = kts, kte
       ! Prediction
       predictExp = min(logMax, max(logMin, output01Activ(1,k)))

       ! Bias correction
       bias_corr = 1.0
       if ((predictExp >= 0.) .and. (predictExp < 3.)) then
          bias_corr = -0.2704*predictExp**5 + 1.838*predictExp**4 - 5.127*predictExp**3 + &
               8.547*predictExp**2 - 8.439*predictExp + 4.297
       endif

       predict_number(k) = bias_corr * (10.**predictExp)
    enddo

  end subroutine predict_number_sub

!------------------------------------------------------------------------------------------------------  
  ! Standard transformer
  subroutine standard_scaler_transform(mean, var, raw_data, transformed_data)
    real, dimension(:,:), intent(in) :: raw_data
    real, dimension(:), intent(in) :: mean, var
    real, intent(out) :: transformed_data(size(raw_data, 1), size(raw_data, 2))
    integer :: i
    
    do i = 1, size(raw_data, 1)
       transformed_data(i,:) = (raw_data(i,:) - mean(i)) / sqrt(var(i))
    end do
    
  end subroutine standard_scaler_transform

!----------------------------------------------------------------------------------
  ! Activation function 
  subroutine relu_activation(input, output)

    real, dimension(:,:), intent(in) :: input
    real, dimension(size(input, 1), size(input, 2)), intent(out) :: output
    integer :: i, j
    
    do i = 1, size(input, 1)
       do j = 1, size(input,2)
          output(i, j) = max(input(i,j), 0.)
       end do
    end do

  end subroutine relu_activation
!----------------------------------------------------------------------------------
  
end module module_mp_tempo_ml
