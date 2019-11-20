module flow_vars
     
integer frame

integer, allocatable ::  T_Bound_type(:)
     
real(8) miu, T_infinity, h_infinity

real(8), allocatable, dimension(:,:) :: Grad_T_face, Grad_T

real(8), allocatable, dimension(:) :: T, T_old, T_face, T_node, T_Bound, Flux_bound
    
end module flow_vars

Module grid_vars

    implicit none

character, allocatable :: face_type(:)

character(len = 2) dummy_char

! dummy_char is a dummy variable to read pre header 

character(len = 40) error_message

! error_message to store the error message

character(len = 4), allocatable :: bndt(:), orf_type(:)

! bndt character array to store boundary types

integer :: nc, nv, nb, version, nf, config(4,12), n_max_bnd

integer i_c, i_nb, i_v, i_v_nb, i_v_nb2, i_v_nb3, i_v_nb4, i_v_nb5, i_v_nb6, i_v_nb7, i_v_nb8, q, i_config, i_config2, i_config3, i_bnd, i_v_bnd, i_v_bnd2, i_v_bnd3, i_v_bnd4, n_undefined, n_orf, i_f, i_f1, i_f2, ncc, n_faces, i_orf_id, n_orf_faces, n_walls

integer, allocatable :: cells(:,:), bnd(:,:), face_nb(:,:), node_nb(:,:), node_nb_bnd(:,:), node_bnd(:), face_nodes(:,:,:), cell_nb(:,:), node_con_cells(:), cell_fnb(:,:), face_cnb(:,:)

! nc is the number of cells
! nv is the number of verticies
! nb is the number of outer boundries
! version is the version of gmte
! nf is the total number of faces
! config for different face arrangments
! Cells is to read the cells data then after processing it is diallocated
! bnd is to read the boundaries data then after processing it is diallocated
! i_c      for cell id counter
! i_nb     for cell neighbour id counter
! i_v      for vertex id counter
! i_v_nb   for neighbour cell vertex id counter
! i_v_nb2  for neighbour cell vertex id second counter
! i_v_nb3  for neighbour cell vertex id third counter
! i_v_nb4  for neighbour cell vertex id fourth counter
! q        for coordinates
! i_c      for boundary id counter
! i_v_bnd  for boundary face vertex id counter
! i_v_bnd2 for secondary boundary face vertex id counter
! i_v_bnd3 for third boundary face vertex id counter
! i_v_bnd4 for fourth boundary face vertex id counter
! n_undefined number of undefined boundary faces
! n_orf number orifice groups
! i_f for face loop counter
! E_f and T_f are two vectors that are calculated such that there sum is the same as the face vector and are used to caculate the cross diffusion term
!ncc is maximum number of cells connected to a single node

real(8), dimension (3)  :: edge_1, edge_2, edge_3, edge_4, center_face
real(8), dimension(0:3) :: cross_1, cross_2

real(8), allocatable :: nodes(:,:), faces(:,:,:), face_center(:,:,:), cell_Vol(:), node_nb_wgh(:,:), cell_center(:,:), r_nb(:,:,:), r_face(:,:,:),  d_nb(:,:), d_face(:,:), face_norm(:,:,:), f_dash(:,:,:), e_centeroids(:,:,:), g_c(:,:), g_c_dash(:,:), r_F_min_f_dash(:,:,:), d_f_dash(:,:), face_center_list(:,:), cell_nb_face_num(:,:), E_f(:,:,:), T_f(:,:,:), e_face(:,:,:), node_bnd_wght(:,:)

logical undefined 

integer, allocatable  :: duplicate(:,:) , orf_face_id(:)

end module grid_vars

module solution_vars

implicit none
integer i, t_step
    
real Voltest

real(8), allocatable :: PMA(:,:), PMB(:)
 
 end module solution_vars
 
 
 Module Linear_Equation_Solvers

contains

subroutine Simple_Equation_solver(MA, MB, Phi)

    use grid_vars
    use flow_vars
    use coefficents_vars

    implicit none

    integer i_iter

    real(8) MA(0:6,nc), MB(nc), Phi(nc), Phi_old(nc), s, P_R(nc)

    P_R = 0.0

    do i_iter = 1, 100
            
            
    Phi_old =  Phi
        
    do i_c = 1, nc

        s = 0.0

        do i_f = 1, 6

        if (cell_fnb (i_f,i_c) .gt. nb) then

        s = s - MA(i_f,i_c) / MA(0,i_c) * Phi(cell_nb(i_f,i_c))
           
        end if    

        end do

        Phi(i_c) =  MB(i_c)/MA(0,i_c) + s
        
        end do
        
        P_R = abs( Phi - Phi_old )

        if ( maxval(abs(P_R/Phi)) .lt. 1e-3) exit
 
     end do  

end subroutine Simple_Equation_solver

end module Linear_Equation_Solvers


    program USD
    
    use grid_vars
    use flow_vars
    use coefficents_vars
    use solution_vars

    implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Reading Grid file**
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

    call read_gmte
    
    call process_gmte

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Initialization**
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    
    
    call initialize_flow_variables

    call green_gauss(T, T_face, Grad_T, Grad_T_face, T_infinity)
 
    call Calc_USG_correction_Vectors

    call write_vtk

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Solution**
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

    do i = 1, 20

    call green_gauss(T, T_face, Grad_T, Grad_T_face, T_infinity)
    
    call calc_Diffusion

    frame = frame + 1

    call write_vtk

    end do

    end program USD
    
    
    
 subroutine abort_job
    
    use grid_vars
    
    Print*, error_message
    
    stop

end subroutine abort_job


subroutine calc_Diffusion
  
    use grid_vars
    use flow_vars
    use solution_vars
    use Linear_Equation_Solvers

    implicit none


    PMA = 0
    PMB = 0

    do i_c = 1 ,  nc

        do i_f = 1, 6

            if ( cell_fnb (i_f,i_c) .gt. nb) then

            PMA(i_f,i_c) = - miu  * E_F (i_f,i_c,0) /  d_nb(i_f,i_c)

            PMA(0,i_c) = PMA(0,i_c) + miu * E_F (i_f,i_c,0) /  d_nb(i_f,i_c)

            PMB (i_c) = PMB (i_c) + miu * ( Grad_T_face(1,cell_fnb (i_f,i_c))*T_F (i_f,i_c,1) + Grad_T_face(2,cell_fnb (i_f,i_c))*T_F (i_f,i_c,2) + Grad_T_face(3,cell_fnb (i_f,i_c))*T_F (i_f,i_c,3) )

            else if ( cell_fnb (i_f,i_c) .lt. nb+1 .and. T_Bound_type(bnd(5,cell_fnb (i_f,i_c))) .eq. 0) then ! Dirichlet Boundary Condition

            PMA(0,i_c) = PMA(0,i_c) + miu * E_F (i_f,i_c,0) /  d_face(i_f,i_c)
            
            PMB (i_c)  = PMB (i_c) + miu * E_F (i_f,i_c,0) /  d_face(i_f,i_c) * T_face(cell_fnb (i_f,i_c)) + miu * ( Grad_T_face(1,cell_fnb (i_f,i_c))*T_F (i_f,i_c,1) + Grad_T_face(2,cell_fnb (i_f,i_c))*T_F (i_f,i_c,2) + Grad_T_face(3,cell_fnb (i_f,i_c))*T_F (i_f,i_c,3) )

            else if ( cell_fnb (i_f,i_c) .lt. nb+1 .and. T_Bound_type(bnd(5,cell_fnb (i_f,i_c))) .eq. 1) then ! Neumann Boundary Condition

            PMB (i_c)  = PMB (i_c) + Flux_bound(bnd(5,cell_fnb (i_f,i_c)))

            else if ( cell_fnb (i_f,i_c) .lt. nb+1 .and. T_Bound_type(bnd(5,cell_fnb (i_f,i_c))) .eq. 2) then ! Mixed Boundary Condition

            PMA(0,i_c) = PMA(0,i_c) + (h_infinity*faces(i_f,i_c,0)*miu*E_F(i_f,i_c,0)/d_face(i_f,i_c))/( h_infinity*faces(i_f,i_c,0) + miu*E_F(i_f,i_c,0)/d_face(i_f,i_c) )

            PMB (i_c)  = PMB (i_c) -  (h_infinity*faces(i_f,i_c,0)*miu*E_F(i_f,i_c,0)/d_face(i_f,i_c))/( h_infinity*faces(i_f,i_c,0) + miu*E_F(i_f,i_c,0)/d_face(i_f,i_c) ) * T_infinity

            PMB (i_c)  = PMB (i_c) -  ( h_infinity*faces(i_f,i_c,0)*miu * ( Grad_T_face(1,cell_fnb (i_f,i_c))*T_F (i_f,i_c,1) + Grad_T_face(2,cell_fnb (i_f,i_c))*T_F (i_f,i_c,2) + Grad_T_face(3,cell_fnb (i_f,i_c))*T_F (i_f,i_c,3) ) ) / ( h_infinity*faces(i_f,i_c,0) + miu*E_F(i_f,i_c,0)/d_face(i_f,i_c) )

            end if

        end do

    end do

call Simple_Equation_solver(PMA, PMB, T)

call calc_node_values
 
end subroutine calc_Diffusion


subroutine calc_face_area
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!
!!  Face center coordinates for each face are stored in face_center
!!
!!  Face unit normal vector is stored in faces array along with the area magnitude  faces(6,nc,(-1,0,1,2,3)) = faces(6-faces, cell_ID, (face id, area magnitude,  normal vector coordinates) )
!!
!!  Face centers are again stored in the face id based list face_center_list(face_id, 1:3 (center coordinates))
!!
!!   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
    use grid_vars
 
    implicit none
    
         
    
    integer xyz
    
    real(8) theta
    
    
    
     
    
   edge_1(1:3) = nodes(1:3,cells(config(2,i_config),i_c)) - nodes(1:3,cells(config(1,i_config),i_c))

   edge_2(1:3) = nodes(1:3,cells(config(4,i_config),i_c)) - nodes(1:3,cells(config(1,i_config),i_c))

   edge_3(1:3) = nodes(1:3,cells(config(4,i_config),i_c)) - nodes(1:3,cells(config(3,i_config),i_c))
                                                
   edge_4(1:3) = nodes(1:3,cells(config(2,i_config),i_c)) - nodes(1:3,cells(config(3,i_config),i_c))
   

    

    
    
    cross_1(1) = edge_1(2) * edge_2(3) -  edge_1(3) * edge_2(2)
    
    cross_1(2) = (edge_1(1) * edge_2(3) -  edge_1(3) * edge_2(1))*-1
  
    cross_1(3) = edge_1(1) * edge_2(2) -  edge_1(2) * edge_2(1)
    
    cross_1(0) = sqrt(cross_1(1)**2 + cross_1(2)**2 + cross_1(3)**2)

    
    
    cross_2(1) = edge_3(2) * edge_4(3) -  edge_3(3) * edge_4(2)
    
    cross_2(2) = (edge_3(1) * edge_4(3) -  edge_3(3) * edge_4(1))*-1
  
    cross_2(3) = edge_3(1) * edge_4(2) -  edge_3(2) * edge_4(1)
    
    cross_2(0) = sqrt(cross_2(1)**2 + cross_2(2)**2 + cross_2(3)**2)


   
    faces(i_config,i_c,0:3) = 0.5 * cross_1(0:3) + 0.5 * cross_2(0:3)
 
    face_norm(i_config,i_c,1:3) = cross_1(1:3) / cross_1(0)


    
    do xyz = 1, 3
    
    face_center(i_config, i_c, xyz) = ( nodes(xyz,cells(config(1,i_config),i_c)) + nodes(xyz,cells(config(2,i_config),i_c)) + nodes(xyz,cells(config(3,i_config),i_c)) + nodes(xyz,cells(config(4,i_config),i_c)) )/4
    
    end do
    
    face_center_list(cell_fnb (i_config,i_c), 1:3) = face_center(i_config, i_c, 1:3)
    
    return
    
    
end subroutine calc_face_area



subroutine calc_node_values
    
    use grid_vars
    use flow_vars
    
    implicit none

    integer bnd_neigbours
    
    T_node = 0

do i_v = 1, nv

    
    if ( node_nb_bnd(1,i_v) .ne. -1000 ) then

        do i_nb = 1, n_max_bnd

            if ( node_nb_bnd(i_nb,i_v) .eq. 0) then

            exit

            else 

  
            T_node(i_v) = T_node(i_v) + node_bnd_wght(i_nb,i_v) * T_face(node_nb_bnd(i_nb,i_v)) 

            end if


        end do

                   
    end if

     

    
    
       if (node_bnd(i_v) .eq. -1000) then
        
                do i_v_nb = 1, ncc
             
                    if(node_nb(i_v_nb,i_v) .eq. 0 ) then 
                    exit
                   else 
                     T_node(i_v) = T_node(i_v) + node_nb_wgh(i_v_nb,i_v) * T(node_nb(i_v_nb,i_v)) 
                     
                   end if
                end do
         
       end if

end do

    
end subroutine calc_node_values


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Caculating the gradient at the cell centers using the Green Gauss method**
!!
!!
!!  ** Gradient values at the faces are then calculated by interpolating cell center gradients**
!!
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 



subroutine green_gauss (Phi, Phi_face, Grad_Phi, Grad_Phi_face, Phi_infinity)

    use grid_vars
    use flow_vars
    use solution_vars

    implicit none

    integer i_iter

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Allocating and intialzing variables**
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

    real(8) Phi(nc), Phi_face(n_faces), Phi_face_dash_2(n_faces), Grad_Phi(1:3,nc), Grad_Phi_old(1:3,nc), Term_1(1:3), Term_2(1:3), Term_3(1:3), Grad_Phi_face(1:3,n_faces), Phi_infinity

    Grad_Phi = 0

    Phi_face_dash_2(1:nb) = Phi_face(1:nb)

!! Option 2


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Step 1: Assume Phi_dash at the mid way between the two centroids and caculate Phi there (stored in Phi_dash_2(n_faces))**
!!
!!  ** Depending on boundary type, boundary faces will have either a fixed Phi value, a value equal to that of the neighbouring cell, or an avaraged value.
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

    

     do i_f = nb+1, n_faces
                
         Phi_face_dash_2(i_f) = (Phi(face_cnb(1,i_f)) + Phi(face_cnb(2,i_f)) ) * 0.5

    end do



    do i_f = 1, nb

        if (T_Bound_type(bnd(5,i_f)) .eq. 1) then        ! Setting face value for Neumann Boundary conditions

            Phi_face_dash_2(i_f) = Phi(face_cnb(1,i_f))

        else if (T_Bound_type(bnd(5,i_f)) .eq. 2) then   ! Setting face value for Mixed Boundary conditions

            Phi_face_dash_2(i_f) = (Phi(face_cnb(1,i_f)) + Phi_infinity)/2

        end if

    end do  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Step 2 calculate the gradient using Phi_dash_2 **
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


    do i_c = 1, nc

        do i_f = 1, 6

            Grad_Phi(1:3,i_c) = Grad_Phi(1:3,i_c) + Phi_face_dash_2(cell_fnb(i_f,i_c)) * faces(i_f,i_c,1:3)

         end do

         Grad_Phi(1:3,i_c) = Grad_Phi(1:3,i_c) / cell_vol(i_c)

    end do    

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Step 3 Calculate the Phi face values**
!!
!!  ** Depending on boundary type, boundary faces will have either a fixed Phi value, a value equal to that of the neighbouring cell, or an avaraged value.
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!     


    do i_iter = 1, 2 ! Steps 3 to 4 will be repeated two times

     do i_f = nb+1, n_faces

         Term_1(1:3) = Grad_Phi(1:3,face_cnb(1,i_f)) + Grad_Phi(1:3,face_cnb(2,i_f))
         
         Term_2(1:3) = face_center_list(i_f, 1:3) - 0.5 * (cell_center(1:3,face_cnb(1,i_f)) + cell_center(1:3,face_cnb(2,i_f)) )
                
         Phi_face(i_f) = Phi_face_dash_2(i_f) + 0.5 *  ( Term_1(1)*Term_2(1) + Term_1(2)*Term_2(2) + Term_1(3)*Term_2(3) ) 

     end do    


    do i_f = 1, nb

        if (T_Bound_type(bnd(5,i_f)) .eq. 1) then        ! Setting face values for Neumann Boundary conditions

            Phi_face(i_f) = Phi(face_cnb(1,i_f))

        else if (T_Bound_type(bnd(5,i_f)) .eq. 2) then   ! Setting face for Mixed Boundary conditions

            Phi_face(i_f) = (Phi(face_cnb(1,i_f)) + Phi_infinity)/2

        end if

    end do  

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Step 4 Update Grad_Phi using Phi at faces this time **
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    
 Grad_Phi  = 0
 
     do i_c = 1, nc

        do i_f = 1, 6

            Grad_Phi(1:3,i_c) = Grad_Phi(1:3,i_c) + Phi_face(cell_fnb(i_f,i_c)) * faces(i_f,i_c,1:3)

         end do

         Grad_Phi(1:3,i_c) = Grad_Phi(1:3,i_c) / cell_vol(i_c)

    end do   
 
  end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Interpolating the Gradient to faces **
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

    Grad_Phi_face = 0

    do i_c = 1, nc

        do i_f = 1, 6

        if(cell_fnb (i_f,i_c) .gt. nb) then

        Term_1(1:3) = Grad_Phi(1:3,i_c) * g_c(i_f,i_c) + Grad_Phi(1:3,cell_nb(i_f,i_c)) * g_c(cell_nb_face_num(i_f,i_c), cell_nb(i_f,i_c))

        Term_2(1:3) =  ( (Phi(cell_nb(i_f,i_c)) - Phi(i_c) ) / d_nb(i_f,i_c)  +  Term_1(1)*e_centeroids(i_f,i_c,1) + Term_1(2)*e_centeroids(i_f,i_c,2) + Term_1(3)*e_centeroids(i_f,i_c,3) ) * e_centeroids(i_f,i_c,1:3)

        Grad_Phi_face(1:3,i_f) = Term_1 + Term_2

        else if(cell_fnb (i_f,i_c) .lt. nb+1 .and. T_Bound_type(bnd(5,i_f)) .eq. 0) then ! Setting face gradient values for Dirichlet Boundary conditions

        Grad_Phi_face(1:3,i_f) = Grad_Phi(1:3,i_c)/2

        else if(cell_fnb (i_f,i_c) .lt. nb+1 .and. T_Bound_type(bnd(5,i_f)) .eq. 1) then ! Setting face gradient values for Neuman Boundary conditions

        Grad_Phi_face(1:3,i_f) = Grad_Phi(1:3,i_c)/2

        else if(cell_fnb (i_f,i_c) .lt. nb+1 .and. T_Bound_type(bnd(5,i_f)) .eq. 2) then ! Setting face gradient values for Mixed Boundary conditions

        Grad_Phi_face(1:3,i_f) = Grad_Phi(1:3,i_c)/2

        end if

        end do

    end do


end subroutine green_gauss



subroutine Calc_USG_correction_Vectors

    use grid_vars
    use flow_vars
    
    implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Intializing Unstructured grid_correction Vector arrays**
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
    
    real(8) cosine_theta 
   
    allocate (E_f(1:6,nc,0:3), T_f(1:6,nc,1:3))

    E_f = 0
    T_f = 0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Over Relaxed Approach**
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!    

    do i_c = 1, nc
        
        do i_f = 1, 6

            if ( cell_fnb (i_f,i_c) .gt. nb) then

                E_F (i_f,i_c,1:3) = ((faces(i_f,i_c,0))**2) / (faces(i_f,i_c,1)*e_centeroids(i_f,i_c,1) + faces(i_f,i_c,2)*e_centeroids(i_f,i_c,2) + faces(i_f,i_c,3)*e_centeroids(i_f,i_c,3)) * e_centeroids(i_f,i_c,1:3)

                cosine_theta = ( E_F(i_f,i_c,1)*faces(i_f,i_c,1) + E_F(i_f,i_c,2)*faces(i_f,i_c,2) + E_F(i_f,i_c,3)*faces(i_f,i_c,3) ) / faces(i_f,i_c,0) / (E_F(i_f,i_c,1)**2 + E_F(i_f,i_c,2)**2 + E_F(i_f,i_c,3)**2)**0.5
            
                T_F (i_f,i_c,1:3) = (face_norm(i_f,i_c,1:3) - 1/cosine_theta * e_centeroids(i_f,i_c,1:3)) * faces(i_f,i_c,0)

             else if ( cell_fnb (i_f,i_c) .lt. nb+1) then

                E_F (i_f,i_c,1:3) = ((faces(i_f,i_c,0))**2) / (faces(i_f,i_c,1)*e_face(i_f,i_c,1) + faces(i_f,i_c,2)*e_face(i_f,i_c,2) + faces(i_f,i_c,3)*e_face(i_f,i_c,3)) * e_face(i_f,i_c,1:3)

                cosine_theta = ( E_F(i_f,i_c,1)*faces(i_f,i_c,1) + E_F(i_f,i_c,2)*faces(i_f,i_c,2) + E_F(i_f,i_c,3)*faces(i_f,i_c,3) ) / faces(i_f,i_c,0) / (E_F(i_f,i_c,1)**2 + E_F(i_f,i_c,2)**2 + E_F(i_f,i_c,3)**2)**0.5

                 T_F (i_f,i_c,1:3) = (face_norm(i_f,i_c,1:3) - 1/cosine_theta * e_face(i_f,i_c,1:3)) * faces(i_f,i_c,0)

            end if


            E_F (i_f,i_c,0) = sqrt(E_F(i_f,i_c,1)**2 + E_F(i_f,i_c,2)**2 + E_F(i_f,i_c,3)**2)

        end do

     end do


end subroutine Calc_USG_correction_Vectors



subroutine process_gmte
    
 use grid_vars   
    
 implicit none
 
 real(8)  vol_m(3,3), total_node

 allocate ( faces(6,nc,-1:3), cell_vol(nc), cell_center(3,nc), cell_nb(6,nc), cell_fnb(6,nc), face_cnb(2,(nc*6-nb)/2+nb), cell_nb_face_num(6,nc) )

 allocate (orf_face_id(nb))
 
 orf_face_id = 0
 i_orf_id    = 0
 cell_center = 0
 cell_nb     = -100
 cell_fnb    = 0
 n_undefined = 0

 config(1,1) = 1
 config(2,1) = 2
 config(3,1) = 6
 config(4,1) = 5

 config(1,2) = 7
 config(2,2) = 3
 config(3,2) = 4
 config(4,2) = 8

 config(1,3) = 7
 config(2,3) = 6
 config(3,3) = 2
 config(4,3) = 3

 config(1,4) = 1
 config(2,4) = 5
 config(3,4) = 8
 config(4,4) = 4

 config(1,5) = 1
 config(2,5) = 4
 config(3,5) = 3
 config(4,5) = 2

 config(1,6) = 7
 config(2,6) = 8
 config(3,6) = 5
 config(4,6) = 6

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!! 
!!   Using one loop over the cells the following is caculated:
!!
!!      1- Centers of all cells and storing them in cell_center(coordinate, cell_id) array 
!!
!!      2- Computing cell volumes and storing them in cell_vol(cell_id) array
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

cell_vol = 0
 
 do i_c = 1, nc
     
!! Internal loop for the cell centers.     
     
     do q = 1, 3

        do i_v = 1, 8

            cell_center(q,i_c) = cell_center(q,i_c) + nodes(q,cells(i_v,i_c))

        end do

        cell_center(q,i_c) = cell_center(q,i_c)/8

     end do
  
!! Caclulating the cell volume.

    vol_m(1,1:3) = nodes(1:3,cells(7,i_c)) - nodes(1:3,cells(1,i_c))
    vol_m(2,1:3) = nodes(1:3,cells(2,i_c)) - nodes(1:3,cells(1,i_c))
    vol_m(3,1:3) = nodes(1:3,cells(3,i_c)) - nodes(1:3,cells(6,i_c)) 
    
    cell_vol(i_c) = vol_m(1,1)*(vol_m(2,2)*vol_m(3,3)-vol_m(3,2)*vol_m(2,3)) - vol_m(2,1)*(vol_m(1,2)*vol_m(3,3)-vol_m(1,3)*vol_m(3,2)) + vol_m(3,1)*(vol_m(1,2)*vol_m(2,3)-vol_m(1,3)*vol_m(2,2))
    
    vol_m(1,1:3) = nodes(1:3,cells(7,i_c)) - nodes(1:3,cells(1,i_c))
    vol_m(2,1:3) = nodes(1:3,cells(5,i_c)) - nodes(1:3,cells(1,i_c)) 
    vol_m(3,1:3) = nodes(1:3,cells(6,i_c)) - nodes(1:3,cells(8,i_c)) 
    
    cell_vol(i_c) = cell_vol(i_c) + vol_m(1,1)*(vol_m(2,2)*vol_m(3,3)-vol_m(3,2)*vol_m(2,3)) - vol_m(2,1)*(vol_m(1,2)*vol_m(3,3)-vol_m(1,3)*vol_m(3,2)) + vol_m(3,1)*(vol_m(1,2)*vol_m(2,3)-vol_m(1,3)*vol_m(2,2))
    
    vol_m(1,1:3) = nodes(1:3,cells(7,i_c)) - nodes(1:3,cells(1,i_c))
    vol_m(2,1:3) = nodes(1:3,cells(4,i_c)) - nodes(1:3,cells(1,i_c)) 
    vol_m(3,1:3) = nodes(1:3,cells(8,i_c)) - nodes(1:3,cells(3,i_c)) 
    
    cell_vol(i_c) = (cell_vol(i_c) + vol_m(1,1)*(vol_m(2,2)*vol_m(3,3)-vol_m(3,2)*vol_m(2,3)) - vol_m(2,1)*(vol_m(1,2)*vol_m(3,3)-vol_m(1,3)*vol_m(3,2)) + vol_m(3,1)*(vol_m(1,2)*vol_m(2,3)-vol_m(1,3)*vol_m(2,2)))/6
    
 end do

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Counting number of Orifice groups and storing them in the integer n_orf **
!!
!!  ** Counting number of Walls groups and storing them in the integer n_walls **
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 
n_orf   = 0
n_walls = 0

do i_bnd = 1, nb

 if (bndt(i_bnd) .eq. "ORFS") then
       
  if (n_orf .lt. bnd(5,i_bnd)) then
           
   n_orf = bnd(5,i_bnd)
       
  end if
 
 end if

 if (bndt(i_bnd) .eq. "WALS") then
       
  if (n_walls .lt. bnd(5,i_bnd)) then
           
   n_walls = bnd(5,i_bnd)
       
  end if
 
 end if
       
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Finding the node with largest number of neighbours**
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

allocate (node_con_cells(nv))

node_con_cells = 0

do i_c = 1, nc
   
    do i_v = 1, 8
       
        node_con_cells(cells(i_v,i_c)) = node_con_cells(cells(i_v,i_c)) + 1
        
    end do

end do

ncc = maxval(node_con_cells)

print*, "Maximum number of cells connected to a single node = ", ncc


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Finding the node with largest number of Boundary neighbours**
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

node_con_cells = 0


do i_f = 1, nb
   
    do i_v = 1, 4
       
        node_con_cells(bnd(i_v,i_f)) = node_con_cells(bnd(i_v,i_f)) + 1
        
    end do

end do

n_max_bnd = maxval(node_con_cells)

print*, "Maximum number of Outer Faces connected to a single node = ", n_max_bnd



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Look for node neighbours**
!!
!!  ** node_nb(ncc,cell_id) is used to store the node's cell neigbours (the maximum number of neigbours will be equalt to ncc) 
!!
!!  ** node_nb_bnd(n_max_bnd,nv) is used to store the node's boundary neigbour ids (the maximum number of neigbours will be equalt to n_max_bnd 
!!
!!  ** node_bnd(node_id) is an integer array to store the node type (internal node, boundary node, shared between two diiferent boundary types)
!!
!!  ** node_nb_wgh(ncc,nv) stores the ids of the cell neighbours (used with internal nodes)
!!
!!  ** node_bnd_wght(n_max_bnd,nv) stores the neigbouring boundary face ids for boundary nodes
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Allocating and intializing node arrays**
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


allocate (node_nb(ncc,nv), node_nb_wgh(ncc,nv), node_bnd(nv), node_nb_bnd(n_max_bnd,nv), node_bnd_wght(n_max_bnd,nv))

node_nb       = 0
node_bnd      = -1000
node_nb_bnd   = 0
node_bnd_wght = 0
node_nb_wgh   = 0.0


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Looking for cell neighbours for internal nodes**
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


do i_c = 1, nc
    
    do i_v = 1, 8
       
        do i_v_nb = 1, ncc
           
            if ( node_nb(i_v_nb,cells(i_v,i_c)) .eq. 0) then
            
            node_nb(i_v_nb,cells(i_v,i_c)) = i_c
            
            exit
            
            end if
            
        end do
                 
    end do

end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Looking for cell neighbours for boundary nodes**
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

do i_bnd = 1, nb
    
    do i_v = 1, 4
       
        do i_v_nb = 1, n_max_bnd
           
            if ( node_nb_bnd(i_v_nb,bnd(i_v,i_bnd)) .eq. 0) then
            
            node_nb_bnd(i_v_nb,bnd(i_v,i_bnd)) = i_bnd
            
            exit
         
            end if
            
        end do
                 
    end do

end do




!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Set node type**
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!   node_bnd(i_v) = :
!! 
!!  -1000 = interior node
!!   1000 = node is shared between a wall and an orifice
!!   -1,-2 = Node is a wall
!!   1,2,.= Node is an orfice number 1,2,...
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


do i_bnd = 1, nb
    
    do i_v = 1, 4
        
      if ( node_bnd(bnd(i_v,i_bnd)) .eq. -1000 ) then
          
            if (bndt(i_bnd) .eq. 'WALS') then
      
            node_bnd(bnd(i_v,i_bnd)) = -1 * bnd(5,i_bnd)
            
            else if (bndt(i_bnd) .eq. 'ORFS') then
        
            node_bnd(bnd(i_v,i_bnd)) = bnd(5,i_bnd)
            
            end if
            
      else if ( node_bnd(bnd(i_v,i_bnd)) .lt. 0 ) then
          
            if (bndt(i_bnd) .eq. 'WALS') then
      
            node_bnd(bnd(i_v,i_bnd)) = -1000 - bnd(5,i_bnd)
            
            else if (bndt(i_bnd) .eq. 'ORFS') then
        
            node_bnd(bnd(i_v,i_bnd)) = 1000 + bnd(5,i_bnd)
            
            end if
           
      else if ( node_bnd(bnd(i_v,i_bnd)) .gt. 0 ) then
          
            if (bndt(i_bnd) .eq. 'WALS') then
      
            node_bnd(bnd(i_v,i_bnd)) = 1000 + bnd(5,i_bnd)
            
            else if (bndt(i_bnd) .eq. 'ORFS') then
                
                if (node_bnd(bnd(i_v,i_bnd)) .eq. bnd(5,i_bnd)) then
        
                 node_bnd(bnd(i_v,i_bnd)) = bnd(5,i_bnd)
                else
                 node_bnd(bnd(i_v,i_bnd)) = 1000 + bnd(5,i_bnd)
                end if
                
            
            end if
            
      end if
      
        
    end do
 
end do

 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** setting node internal neighbours weight ratios**
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 


   do i_v = 1, nv
        
        total_node = 0
       
        do i_v_nb = 1, ncc
           
            if ( node_nb(i_v_nb,i_v) .gt. 0) then
            
                node_nb_wgh(i_v_nb,i_v) =  1/sqrt((cell_center(1,node_nb(i_v_nb,i_v)) - nodes(1,i_v))**2 + (cell_center(2,node_nb(i_v_nb,i_v)) - nodes(2,i_v))**2 + (cell_center(3,node_nb(i_v_nb,i_v)) - nodes(3,i_v))**2)
                
                total_node = total_node + node_nb_wgh(i_v_nb,i_v)
                
                    if (i_v_nb .eq. ncc) then
                       
                        do i_v_nb2 = 1, i_v_nb
                    
                        node_nb_wgh(i_v_nb2,i_v) =  node_nb_wgh(i_v_nb2,i_v) / total_node
                
                        end do 
  
                    end if
                    
            
            else if ( node_nb(i_v_nb,i_v) .eq. 0) then
          
                do i_v_nb2 = 1, i_v_nb
                    
                node_nb_wgh(i_v_nb2,i_v) =  node_nb_wgh(i_v_nb2,i_v) / total_node
                
                end do
                
                exit
            end if
            
        end do
                 
   end do


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** cell connectiviety and external boundaries **
!!
!!  1- i_f is a counter for the unique internal faces
!!  2- cell_fnb(cell face(1-6), cell_id) is an array storing neighbour faces
!!  3- Each face has two cell neibours stored in face_cnb(neigbouring cell(1-2),face_id)
!!  4- cell_nb_face_num(cell face(1-6),cell_id) to store at each face the original cell is located from the neighbour's prespective
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
 
i_f = 0
 
 do i_c = 1, nc
     
     do i_config = 1, 6     !! Looping throught the six possible vertex convfigurations for the face in GMTE format
         
         undefined = 1      !! Boolean to check failure in finding the neighbour

         do i_v = 1, ncc    !! Looking for the neighbours in the set of the cells that are connected to the same node
             
             if (node_nb(i_v,cells(config(1,i_config),i_c)) .eq. 0) exit    !! check if the node has a cell neighbour or not
             
             i_nb = node_nb(i_v,cells(config(1,i_config),i_c))              
             
             if (i_c .ne. i_nb) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Four Nested loops looking for a face match in the first cell neghbour picked from the previous step **
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
                                
              do i_v_nb = 1, 8                               
     
                if ( cells(config(1,i_config),i_c) .eq. cells(i_v_nb,i_nb) ) then
            
                 do i_v_nb2 = 1, 8 
             
                  if ( cells(config(2,i_config),i_c) .eq. cells(i_v_nb2,i_nb) ) then
                 
                   do i_v_nb3 = 1, 8 
                     
                    if ( cells(config(3,i_config),i_c) .eq. cells(i_v_nb3,i_nb) ) then
                                    
                     do i_v_nb4 = 1, 8 
                                     
                      if ( cells(config(4,i_config),i_c) .eq. cells(i_v_nb4,i_nb) ) then
                                            
                       cell_nb(i_config,i_c) = i_nb            !! Storing the cell id as a neighbour in case of sucess                                      
                                                
                       if (cell_fnb (i_config,i_c) .eq. 0)  then

                        i_f = i_f + 1                          !! Incrementing the number of internal unique faces

                        cell_fnb(i_config,i_c) = nb+i_f        !! storing the face id in the cell's face neighbour array
                                                    
                        face_cnb(1,nb+i_f) = i_c               !! storing the cell id in the faces's cell neighbour array

                       end if
                                                
                       undefined = 0
                                                
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** To avoid internal face duplication another four Nested loops to store the same face id on the neghbour's **
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!                               
                                                
                      do i_config3 = 1, 6

                       do i_v_nb5 = 1, 8
     
                        if ( cells(config(1,i_config3),i_nb) .eq. cells(i_v_nb5,i_c) ) then
            
                         do i_v_nb6 = 1, 8 
             
                          if ( cells(config(2,i_config3),i_nb) .eq. cells(i_v_nb6,i_c) ) then
                 
                           do i_v_nb7 = 1, 8 
                     
                            if ( cells(config(3,i_config3),i_nb) .eq. cells(i_v_nb7,i_c) ) then
                                    
                             do i_v_nb8 = 1, 8 
                                     
                              if ( cells(config(4,i_config3),i_nb) .eq. cells(i_v_nb8,i_c) ) then
                                                                                     
                               if (cell_fnb (i_config3,i_nb) .eq. 0) then

                                cell_fnb (i_config3,i_nb) = nb+i_f              !! storing the face id at the 2nd cell entry in the face neighbour array
                                                                                         
                                face_cnb(2,nb+i_f) = i_nb                       !! storing the 2nd cell id in the faces's cell neighbour array

                               end if

                                 cell_nb_face_num(i_config,i_c) = i_config3     !! storing the face number from which the 2nd cell sees the oringinal cell

                               end if
                              end do
                             end if
                            end do
                           end if
                          end do
                         end if
                        end do
                       end do
                      end if
                     end do
                    end if
                   end do
                  end if
                 end do                 
                end if
               end do
              end if
             end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Failed to find an internal neigbour, looking for the boundary face **
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

             
        
        if (undefined .eq. 1)  then     
         
         do i_v = 1,n_max_bnd
  
          if (node_nb_bnd(i_v,cells(config(1,i_config),i_c)) .eq. 0) exit
             
           i_bnd = node_nb_bnd(i_v,cells(config(1,i_config),i_c))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Four Nested loops looking for a face match in the boundary faces **
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
           
             do i_v_bnd = 1, 4
            
              if ( cells(config(1,i_config),i_c) .eq. bnd(i_v_bnd,i_bnd) ) then
                 
               do i_v_bnd2 = 1, 4
                     
                if ( cells(config(2,i_config),i_c) .eq. bnd(i_v_bnd2,i_bnd) ) then 
                         
                 do i_v_bnd3 = 1, 4
                             
                  if ( cells(config(3,i_config),i_c) .eq. bnd(i_v_bnd3,i_bnd) ) then 
                                 
                   do i_v_bnd4 = 1, 4 
                                    
                    if ( cells(config(4,i_config),i_c) .eq. bnd(i_v_bnd4,i_bnd) ) then 
                                        
                     if (bndt(i_bnd).eq. "WALS") then

                      cell_fnb (i_config,i_c) = i_bnd       !! storing the face id of the matching boundary face

                      cell_nb(i_config,i_c) = 0             !! Adding zero to that niegbour place in the cell's neighbour array

                      undefined = 0
                                       
                      face_cnb(1,i_bnd) = i_c               !! storing the cell id in the faces's cell neighbour array

                      face_cnb(2,i_bnd) = 0                 !! Adding zero to the faces's cell second neighbour ( Zero indicates it is a wall)
                                       
                     else if (bndt(i_bnd).eq. "ORFS") then
                                       
                       cell_fnb (i_config,i_c) = i_bnd      !! storing the face id of the matching boundary face
                                       
                       undefined = 0

                       face_cnb(1,i_bnd) = i_c              !! storing the cell id in the faces's cell neighbour array

                       face_cnb(2,i_bnd) = -1               !! Adding -1 to the faces's cell second neighbour ( -1 indicates it is n orifice)
                                       
                       i_orf_id = i_orf_id + 1

                       orf_face_id(i_orf_id) = i_bnd
                                       
                       end if
                      end if
                     end do
                    end if
                   end do
                  end if
                 end do
                end if
               end do
              end do
             end if
            if (undefined) then

             n_undefined = n_undefined +1

           end if
          end do
         end do
 
 n_faces = nb+i_f               !! storing the total number of unique faces
 
 n_orf_faces = i_orf_id         !! storing the total number of oriface faces
 
 print*, "Failed to find", n_undefined, "Cell neigbours"
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Computing faces areas and directions**
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
!!
!!  Face area vectors and normal vectors are stored in 
!!
!!  faces(6,nc,(-1,0,1,2,3)) = faces(face id, area magnitude, unit normal vector coordinates)
!!  
!!  face_nb(4,i_f) = face_nb((cell_neighbour-1, face_order, cell_neigbour-2, face order),i_f)
!!
!!
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 



allocate (face_center(6,nc,3), face_nodes(4,nc,6), face_norm(6,nc,3), face_center_list(n_faces, 1:3))

do i_c = 1, nc
    
    do i_config = 1, 6
                                    
   call calc_face_area
   
   faces(i_config,i_c,-1) = cell_fnb (i_config,i_c)
                                           
   end do
end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Cell neighbours and faces distances **
!!
!!  1- d_nb(cell face, cell id) array containing the distance between cell center and neigbour centers (D_CF) (If the neighbour is a boundary it will be set to zero)
!!  2- r_nb(cell face, cell id, three components) array containing the vector between cell center and neigbour centers (D_CF) (If the neighbour is a boundary it will be set to zero) 
!!  3- f_dash is the location of the intersection between line connecting the two cell centeroids and the shared face.
!!  4- g_c is a geomertric factor used in face value calculation
!!  5- g_c_dash is a geomertric factor used in face value calculation
!!  6- r_F_min_f_dash(1:6,cell_id,components) is the just f_dash - r_nb
!!  7- d_f_dash(1:6, cell_id) is the magnitude of r_F_min_f_dash and is used to claulate g_c_dash
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

allocate (d_nb(6,nc), d_face(6,nc), r_nb(6,nc,3), r_face(6,nc,3), f_dash(1:6,nc,1:3), e_centeroids(1:6,nc,1:3), g_c(1:6,nc), g_c_dash(1:6,nc), r_F_min_f_dash(1:6,nc,1:3), d_f_dash(1:6,nc), e_face(1:6,nc,1:3) )

d_nb           = 0
d_face         = 0
r_nb           = 0
r_face         = 0
f_dash         = face_center
g_c            = 0
g_c_dash       = 0
r_F_min_f_dash = 0
e_centeroids   = 0
d_f_dash       = 0
e_face         = 0

do i_c = 1, nc
    
    do i_nb = 1, 6

     r_face(i_nb,i_c,1:3) = face_center(i_nb,i_c,1:3) - cell_center(1:3,i_c)                            !! The vector connecting cell centroids
        
     d_face(i_nb,i_c)     = sqrt(r_face(i_nb,i_c,1)**2+r_face(i_nb,i_c,2)**2+r_face(i_nb,i_c,3)**2)     !! The magnitude of the vector connecting cell centroids

     e_face(i_nb,i_c,1:3) =  r_face(i_nb,i_c,1:3)/d_face(i_nb,i_c)                                      !! The unit vector of the vector connecting cell centroids
    
    if (cell_nb(i_nb,i_c) .gt. 0) then

     r_nb(i_nb,i_c,1:3)           = cell_center(1:3,cell_nb(i_nb,i_c)) - cell_center(1:3,i_c)           !! The Vector from cell center to cell neigbour center

     d_nb(i_nb,i_c)               = sqrt(r_nb(i_nb,i_c,1)**2+r_nb(i_nb,i_c,2)**2+r_nb(i_nb,i_c,3)**2)   !! The Magnitude of the previous vector

     e_centeroids(i_nb,i_c,1:3)   = r_nb(i_nb,i_c,1:3) / d_nb(i_nb,i_c)                                 !! Unit vector of the previous vector

     f_dash(i_nb,i_c,1:3)         = ((face_center(i_nb,i_c,1) * faces(i_nb,i_c,1) + face_center(i_nb,i_c,2) * faces(i_nb,i_c,2) + face_center(i_nb,i_c,3) * faces(i_nb,i_c,3))-(cell_center(1,i_c) * faces(i_nb,i_c,1) + cell_center(2,i_c) * faces(i_nb,i_c,2) + cell_center(3,i_c) * faces(i_nb,i_c,3)) )/ (e_centeroids(i_nb,i_c,1)*faces(i_nb,i_c,1) + e_centeroids(i_nb,i_c,2)*faces(i_nb,i_c,2) + e_centeroids(i_nb,i_c,3)*faces(i_nb,i_c,3)) * e_centeroids(i_nb,i_c,1:3) + cell_center(1:3,i_c) !! Point of intersection between the face and the line connecting cell centoids

     r_F_min_f_dash(i_nb,i_c,1:3) = cell_center(1:3,cell_nb(i_nb,i_c)) - f_dash(i_nb,i_c,1:3)           !! Vector created by subtracting the r_nb - f_dash

     d_f_dash(i_nb,i_c)   = sqrt(r_F_min_f_dash(i_nb,i_c,1)**2+r_F_min_f_dash(i_nb,i_c,2)**2+r_F_min_f_dash(i_nb,i_c,3)**2) !! Magnitude of the previous vector

     g_c(i_nb,i_c)        = ( d_face(i_nb,i_c) / d_nb(i_nb,i_c) )                                       !! Geometeric factor for the two cell centers based on face center
     
     g_c_dash(i_nb,i_c)   = ( d_f_dash(i_nb,i_c) / d_nb(i_nb,i_c) )                                     !! Geometeric factor for the two cell centers based on f_dash

    end if
   end do
 end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** setting boundary nodes face neighbours weight ratios**
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

     do i_v = 1, nv
        
        total_node = 0

        if (node_bnd(i_v) .ne. -1000) then

            do i_v_nb = 1, n_max_bnd

                if ( node_nb_bnd(i_v_nb,i_v) .gt. 0) then

                node_bnd_wght(i_v_nb,i_v) =  1/sqrt((face_center_list(node_nb_bnd(i_v_nb,i_v),1) - nodes(1,i_v))**2 + (face_center_list(node_nb_bnd(i_v_nb,i_v),3) - nodes(2,i_v))**2 + (face_center_list(node_nb_bnd(i_v_nb,i_v),3) - nodes(3,i_v))**2)

                total_node = total_node + node_bnd_wght(i_v_nb,i_v)

                else if ( node_nb_bnd(i_v_nb,i_v) .eq. 0) then

                exit 

                end if

            end do

               do i_v_nb2 = 1, n_max_bnd

                if ( node_nb_bnd(i_v_nb2,i_v) .gt. 0) then
            
                node_bnd_wght(i_v_nb2,i_v) = node_bnd_wght(i_v_nb2,i_v) / total_node

               end if
            end do
        end if
     end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Printing grid data**
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

print*, "Number of cells = ", nc
print*, "Number of nodes = ", nv
print*, "Number of Boundary faces = ", nb
print*, "Number of internal faces = ", n_faces-nb
print*, "Number of Orifice groups = ", n_orf

    
end subroutine process_gmte



subroutine initialize_flow_variables
    
    use grid_vars
    use flow_vars
    use solution_vars
        
    implicit none

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Alocating and intializing cell, face, nodes, and gradient values**
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 
    allocate (T_Bound_type(n_walls))
    allocate (T(nc), T_old(nc), T_face(n_faces), T_node(nv), T_Bound(n_walls), Grad_T(1:3,nc), Grad_T_face(1:3,n_faces), Flux_bound(n_walls))
 
    T      = 20
    T_node = 0
    T_face = 20
    Grad_T = 0
    frame  = 1000


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Setting Boundary Values**
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    T_Bound(1) = 20
    T_Bound(2) = 100
!   T_Bound(3) = 10

    do i_f = 1, nb

        T_face(i_f) = T_Bound(bnd(5,i_f))

    end do


    T_Bound_type(1) = 2
    T_Bound_type(2) = 0
!   T_Bound_type(3) = 1

    Flux_bound = 0
    T_infinity = 20
    h_infinity = 2
       
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Setting nodal values**
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

 call calc_node_values


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Allocating Algebric equation cooefficents Matrices **
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

     allocate(  PMA(0:6,nc), PMB(nc))
    
       PMA = 0

       PMB = 0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Fluid Properties**
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

    miu = 237 
  
end subroutine initialize_flow_variables


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Reading Grid file**
!!
!!  **  USD can read the grid file with the gmte format (ACFlux), it can b created using ICEM CFD grid genreator**
!!
!!  ** For now it is limited to the file name "test.1.gmte.c", the grid file should be named like that and placed in the debug folder**
!!
!!  ** The GMTE file has a header followed by three main sections, the first one lists the cell data, the second one for verticies, and the last one for boundary faces, these three sections are read and stored in the grid arrays
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!   

subroutine read_gmte
 
  use grid_vars
  
  implicit none
  
  integer i

    open(unit=20, file='Cells')
    open(unit=21, file='Nodes')
    open(unit=22, file='Boundaries')
    open(unit=23, file='Boundary types')
    open(unit=10, file='test.1.gmte.c')

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Reading the file Header **
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    
    read(10,*) 
    read(10,*)
    read(10,*) dummy_char, nc, nv, nb

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Allocating the grid arrays**
!!
!!  - cells array to store cell data (vertex ids for each cell, cell region id,.....)
!!  - nodes to store the vertecies coordinates. 
!!  - bnd to store the boundary faces data (verticex ids of each boundary face, boundary group id)
!!  - bndt to store the boundary type (walls, Orifice, symmetry,....)
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  

    allocate(cells(0:16,nc), nodes(0:3,nv),bnd(0:8,nb),bndt(nb))

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Reading the cell data**
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  
    
    do i = 1, nc
        read(10,*) dummy_char, cells(0,i), cells(1,i), cells(2,i), cells(3,i), cells(4,i), cells(5,i), cells(6,i), cells(7,i), cells(8,i), cells(9,i), cells(10,i), cells(11,i), cells(12,i), cells(13,i), cells(14,i), cells(15,i), cells(16,i) 
        if (dummy_char.ne."CD") then
        Print*, "Syntax error at cell number", i
        stop
        end if
    end do
    
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Reading the vertex data**
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   
    do i = 1, nv
        read(10,*) dummy_char, nodes(0,i), nodes(1,i), nodes(2,i), nodes(3,i)
        if (dummy_char.ne."VD") then
        Print*, "Syntax error at Vertex number", i
        stop
        end if
    end do

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Reading the boundary faces data**
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

    do i = 1, nb
        read(10,*) dummy_char, bnd(0,i), bnd(1,i), bnd(2,i), bnd(3,i), bnd(4,i), bnd(5,i), bnd(6,i), bnd(7,i), bnd(8,i), bndt(i)
        if (dummy_char.ne."BD") then
        Print*, "Syntax error at Boundary number", i
        stop
        end if
    end do

    end subroutine read_gmte
    
    
    subroutine write_vtk

    use grid_vars
    use flow_vars
    
      
    implicit none
    
    character (len = 1024) ::  plot
    
    write(plot,"(A7,I4,A4)")  "ABYtest", frame, ".vtk"
    
  
   open(unit=31, file=plot)
   
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Writing file header and geometry**
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   
   
   write (31, 501)
   write (31, 502) 
   write (31, 503) 
   write (31, 504) 
   write (31, 505) nv
   
   
    do i_v= 1, nv
    write (31,601) nodes(1,i_v), nodes(2,i_v), nodes(3,i_v)
    end do
    
      
    write (31,506) nc,  nc*9
    
    
    
    do i_c = 1, nc
    write (31,602) cells(1,i_c)-1, cells(2,i_c)-1, cells(3,i_c)-1, cells(4,i_c)-1, cells(5,i_c)-1, cells(6,i_c)-1, cells(7,i_c)-1, cells(8,i_c)-1
    end do
    
    
    write (31, 507)  nc
    
    do i_c = 1, nc
    write (31,603) 12
    end do
    
    
    

 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Writing Values of Temperature**
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
   
    write (31, 508)  nv
    write (31, 530)
    write (31, 510)
    
      do i_v = 1, nv
             
          write (31,604) T_node(i_v)

      end do 
      
      
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Writing Values of Temperature Cell Gradients**
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

    write (31, 515)  nc
    write (31, 531)
    write (31, 510)

     do i_c = 1, nc
             
          write (31,604) Grad_T(1,i_c)

      end do 


    write (31, 532)
    write (31, 510)

     do i_c = 1, nc
             
          write (31,604) Grad_T(2,i_c)

      end do 
    
     write (31, 533)
     write (31, 510)

     do i_c = 1, nc
             
          write (31,604) Grad_T(3,i_c)

      end do 

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!  ** Writing Values of Temperature Cell **
!! 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

     write (31, 530)
     write (31, 510)

     do i_c = 1, nc
             
          write (31,604) T(i_c)

      end do 
   

   
501 format("# vtk DataFile Version 2.0" )   
502 format("Unstructured Grid Example")   
503 format("ASCII")       
504 format("DATASET UNSTRUCTURED_GRID")       
505 format("POINTS ", i10, " float")   
506 format("CELLS ", i10, " ", i10)  
507 format("CELL_TYPES ", i10)
508 format("POINT_DATA ", i7)
509 format("SCALARS Pressure_(Pa) float 1")
510 format("LOOKUP_TABLE default")
511 format("SCALARS u_(m/sec) float 1")
512 format("SCALARS v_(m/sec) float 1")
513 format("SCALARS w_(m/sec) float 1")    
514 format("VECTORS Velocity_(m/sec)  float")
515 format("CELL_DATA ", i7)
516 format("SCALARS Pressure_(Cell) float 1")
517 format("SCALARS u_(Cell) float 1")
518 format("SCALARS v_(Cell) float 1")
519 format("SCALARS w_(Cell) float 1")    
520 format("VECTORS Velocity_(Cell)  float")  
    
521 format("SCALARS v_(star) float 1")
522 format("SCALARS v_(corr) float 1")
523 format("SCALARS p_(corr) float 1")
    
530 format("SCALARS Temperature_(K) float 1")
531 format("SCALARS Temperature_Gradient_x float 1")
532 format("SCALARS Temperature_Gradient_y float 1")
533 format("SCALARS Temperature_Gradient_z float 1")
    
    
601 format(f10.5, " ",f10.5, " ",f10.5)
602 format("8 ", i7, " ",i7, " ",i7, " ",i7," ",i7," ",i7," ",i7," ",i7)
603 format(i2)
604 format(f20.10)
 
 end subroutine write_vtk



