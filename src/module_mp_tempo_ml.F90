! TEMPO neural network module designed for cloud droplet number concentration prediction
! A big thanks to David John Gagne (NCAR) for very useful discussions
! Please see https://github.com/NCAR/mlmicrophysics for a more detailed way to connect
! ML models in python to Fortran microphysics code

module module_mp_tempo_ml
  use module_mp_tempo_params, only : wp, sp, dp
  
  implicit none
  ! public :: MLdata, tempo_save_or_read_ml_data, predict_number, predict_number_sub
   public 
   
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

  ! ML data
    integer, parameter :: nc_ml_input = 7
    integer, parameter :: nc_ml_nodes = 24
    integer, parameter :: nc_ml_output = 1

    integer, parameter :: nr_ml_input = 7
    integer, parameter :: nr_ml_nodes = 24
    integer, parameter :: nr_ml_output = 1

    real(wp), dimension(nc_ml_input), parameter :: &
         nc_ml_trans_mean = (/0.000191556196486247, 3.58145042772654e-05, &
         3.12611085359273e-07, 5.74303078738579e-05, 84191.2092225319, &
         279.070551773565, 0.123679354084004/)
    real(wp), dimension(nc_ml_input), parameter :: &
         nc_ml_trans_var = (/5.78143564171777e-08, 3.22834309750552e-08, &
         6.45745893307455e-11, 4.16625579383794e-08, 215694631.771185, &
         94.6576255386858, 0.384841247662964/)

    real(wp), dimension(nc_ml_input * nc_ml_nodes), parameter :: &
         nc_ml_w00 = (/-2.006957, -0.2812008, -0.339073, 1.596426, 2.395225, 1.76315, &
         -0.0626798, 1.267002, -0.02234177, -6.522605e-33, 0.4792154, -0.1253034, &
         -3.217191, -3.092887, -0.0863651, 1.071625, 0.09741028, 0.2255831, &
         -0.6929023, -0.02693799, -3.432344e-33, -0.8791879, -0.9359049, 1.083484, &
         -0.07909214, -0.0122418, -0.02815927, 0.1676407, 0.08252326, 0.6697816, &
         -0.4019359, 0.4687141, 0.001813132, -9.792186e-33, 0.0409322, 0.0113192, &
         -0.01354596, 0.00307771, -0.4635534, 0.03835761, -0.1015553, 0.7316446, &
         -0.05791711, -0.0002690362, -7.920147e-33, -0.1216918, -0.3190572, 0.09809405, &
         -0.16476, -0.03387314, 0.005422261, 0.04043967, 0.03901243, 0.07444729, &
         0.01954299, 0.06918761, 0.04823543, -8.637957e-33, 0.06371575, -0.09250915, &
         -1.109653, -1.373999, -0.2412623, -0.04482195, 0.1584691, 0.06353725, &
         0.0006248798, 0.04593191, -8.878673e-33, -0.4988684, 0.01110262, 0.04623203, &
         0.006581791, 0.03536217, -0.1890567, -0.08839592, 0.1327181, 0.03478973, &
         -0.1565902, -0.100401, -0.1179777, -8.879818e-33, -0.1383738, 0.02847495, &
         -0.005902881, 0.005615512, -0.6308192, -0.02431803, -0.141971, -0.3490018, &
         -0.9850957, -0.1449479, -8.059166e-33, -0.1186465, -1.165381, 0.069015, &
         0.003388841, -0.04041302, 0.1638467, 0.1147008, -0.04833491, -0.07755993, &
         -0.5137688, 0.04546477, 0.04101883, 6.752353e-33, 0.3541977, 0.04880851, &
         -0.00102834, -0.01280629, -0.1116254, -0.02204754, 0.07100908, 0.2354002, &
         0.07129629, 0.2489657, 8.080785e-33, -0.03449865, 0.06037927, -0.02023619, &
         0.7779589, 0.04680278, 0.7492616, 0.6545208, -1.09497, -1.176524, &
         -0.451585, 0.881124, -0.4551499, -8.708624e-33, 1.006558, -0.04979523, &
         -0.0006915367, -0.002993054, 0.01654614, -0.07141764, -0.2216591, 0.8637336, &
         0.8358089, -0.7576646, 8.186339e-33, 0.1161914, 0.7121871, -1.146734, &
         0.03925339, 0.9975697, -0.06953461, -0.07598846, -0.06418022, -0.01897495, &
         -0.03612464, -0.06703389, -0.103049, 9.364858e-33, -0.03949085, 1.005938, &
         0.0001625419, -0.01027657, 0.03823901, 0.3197246, -0.1160718, 0.04449431, &
         0.009831783, 0.6551342, -8.470687e-33, 0.1374123, -0.01782138, 0.01719949/)
    real(wp), dimension(nc_ml_nodes), parameter :: &
         nc_ml_w01 = (/-4.045869, 0.558111, -1.567351, -1.64972, 2.748608, &
         -1.909901, 0.3955558, 1.507247, 0.3599722, -0.0001173223, -0.9398569, &
         -0.6867028, -61.72853, -63.13766, 1.165811, -0.6848684, 0.1931683, &
         1.1208, 2.63087, 0.740169, -21.62499, 1.545568, 3.575141, -1.299604/)
    real(wp), dimension(nc_ml_nodes), parameter :: &
         nc_ml_b00 = (/-0.9842531, 0.3064759, -0.4500185, 1.28336, 1.384105, 0.528031, &
         0.8453538, 1.579872, 2.245679, -0.008679952, 0.4549862, -0.136581, &
         -2.576741, -2.483647, -0.2089484, 0.7607977, 1.847745, 0.7316047, &
         -0.287945, 2.227298, -2.314714, -0.2561245, -0.6993448, -0.1359731/)
    real(wp), dimension(nc_ml_output), parameter :: &
         nc_ml_b01 = (/1.572826/)



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
