!      Ghanim Ullah, John R Cressman Jr, Ernest Barreto, and Steven J Schiff,
!                        July 29, 2008,
!                        Penn State University,
!                        Email: ghanim.phy@gmail.com
!      Archived with 
!      "The Influence of Sodium and Potassium Dynamics on Excitability, Seizures, 
!      and the Stability of Persistent States: II. Network and Glia Dynamics" (2009) 
!      Journal of Computational Neuroscience, 26:171-183
!      This program couples 100 inhibitory neurons and 100 excitatory neurons where the membrane potential dynamics of 
!      these neurons is taken from Gutkin et al. model, 2001,
!      J. Computational Neuroscience, 11, 121-134. 
!      The synaptic currents here are modified from that given in Gutkin et al. model, 2001 
!      The model also includes dynamic potassium and sodium concentrations
!      that builds on the model from companion paper "The Influence of Sodium and Potassium Dynamics on Excitability, Seizures, 
!      and the Stability of Persistent States: I. Single Neuron Dynamics" by
!      John R Cressman Jr, Ghanim Ullah, Jokubas Ziburkus, Steven J Schiff , and Ernest Barreto (2009),
!      Journal of Computational Neuroscience, 26:159-170.
!      NOTE: It is important to start the simulation with steady state dynamics in order to reproduce 
!      various results for given parameters set. The single cell version 
!      of this code is also archived with Cressman et al., (2009) J. Comp. Neurosci. 26:159-170. 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 

       implicit none
       integer i, j
	   integer N, N1, dim                                   !number of cells from each cell type, dimensions of each cell, number of data points
       Parameter(N=100, dim=15, N1=10000000)
       double precision time, time_step                     !total time, time step for integration
	   double precision derivs(N*dim), x(N*dim), thresh     !thresh = spiking threshold
       double precision activity_ex, activity_in            !activity of excitatory cells, activity of inhibitory cells
       integer ee, bb      
       EXTERNAL rk4                                         !rk4 algorithm is used for integration            
       EXTERNAl SystemDerives                               !function for the dynamics of the network
       time_step = 0.01d0                                   !time step (in ms)
       thresh = 0.0                                         !spiking threshold
       time = 0.0                                           !total time (in ms)

       open(unit=1, file='raster_excitatory.dat')           !saves data for raster plot of excitatory network (huge file) 
       open(unit=2, file='raster_inhibitory.dat')           !saves data for raster plot of inhibitory network (huge file)
       open(unit=3, file='activity_excitatory.dat')         !saves data for activity of excitatory network
       open(unit=4, file='activity_inhibitory.dat')         !saves data for activity of inhibitory network
       open(unit=5, file='v_cell50_excitatory.dat')         !saves membrane potential of cell 50 in excitatory network 
       open(unit=6, file='ko_cell50_excitatory.dat')        !saves extracellular potassium for cell 50 in excitatory network
       open(unit=7, file='Nai_cell50_excitatory.dat')       !saves intracellular sodium for cell 50 in in excitatory network
	   open(unit=8,file='time.dat')                         !saves total time 

!      initial conditions
       do i = 1, N
          x((i-1)*dim+1) = -65.1034452264476897             !membrane potential of excitatory network
          x((i-1)*dim+2) = 0.645021617064742286d-01         !potassium channel activation variable n for excitatory network
          x((i-1)*dim+3) = 0.981306641820766101             !sodium channel inactivation variable h for excitatory network
          x((i-1)*dim+4) = 0.217441399671948571d-05         !variable s for temporal evolution of synaptic efficacy emanating from excitatory network 
          x((i-1)*dim+5) = 0.978772044450795857d-07         !calcium concentration for excitatory network
          x((i-1)*dim+6) = -65.1034452264476897             !membrane potential of inhibitory network
          x((i-1)*dim+7) = 0.645021617064742286d-01         !potassium channel activation variable n for excitatory network
          x((i-1)*dim+8) = 0.981306641820766101             !sodium channel inactivation variable h for inhibitory network
          x((i-1)*dim+9) = 0.217441399671948571d-05         !variable s for temporal evolution of synaptic efficacy emanating from inhibitory network   
          x((i-1)*dim+10) = 7.22                            !extracellular potassium concentration for excitatory network
          x((i-1)*dim+11) = 18.54                           !intracellular sodium concentration for excitatory network
          x((i-1)*dim+12) = 7.22                            !extracellular potassium concentration for inhibitory network
          x((i-1)*dim+13) = 18.54                           !intracellular sodium concentration for inhibitory network
          x((i-1)*dim+14) = 0.0								!variable eta, modeling the synaptic block due to the depolarization for excitatory network
          x((i-1)*dim+15) = 0.0								!variable eta, modeling the synaptic block due to the depolarization for inhibitory network  
       end do 
 
!     main program, calls the subroutines SystemDerives and rk4 to integrate and tracks the activity of the network
 
       do i = 1, N1/10                                      !begin time loop
          activity_ex = 0.0                                 !activity of the excitatory cells
          activity_in = 0.0                                 !activity of the inhibitory cells
          do j = 1, 10                                      !loops the rk4 method 10 times before recording a data point, can be changed
             time = time + time_step
             call SystemDerives(time,x,derivs)
             call rk4(x,derivs,N*dim,time,time_step,x,SystemDerives)
          end do
          do j = 1, N
             ee = 0
             bb = 0
             if(x((j-1)*dim+1) .gt. thresh)then
               ee = j
               write(1,*)time,ee
               activity_ex = activity_ex + 1.0
             end if
             if(x((j-1)*dim+6) .gt. thresh)then
               bb = j+N
               write(2,*)time,bb
               activity_in = activity_in + 1.0
             end if
          end do 
          write(3,*)activity_ex
          write(4,*)activity_in
		  write(5,*)x(49*dim+1)
		  write(6,*)x(49*dim+10)
		  write(7,*)x(49*dim+11)
		  write(8,*)time
       end do                                               !end time loop
       end
             
!      This subroutine gives the network dynamics

       SUBROUTINE SystemDerives(x,y,dydx)
       implicit none
       double precision PI, PIOVERTWO, TWOPI
       Parameter (PI = 3.14159265358979323846,PIOVERTWO = PI/2.0d0, TWOPI = 2.0d0 * PI)
       double precision ran2
       integer N, dim, ii, jj
       parameter(N=100, dim=15)

!	   Various parameters and variables used in the model are declared here

       double precision C, phi, tau_e, A, g_ca, V_ca, g_l, V_l, g_na, V_na_e,V_na_i, g_k, g_ahp,V_na
       double precision V_k, V_k_e, V_k_i, tau_i, V_ee, V_ie, V_ei, V_ii, alpha_ee, alpha_ie
       double precision alpha_ei, alpha_ii,alpha_g, g_ee, g_ie, g_ei, g_ii,g_g
       double precision dx1,dx2, specnum1, specnum2, synsum1, synsum2,gapsum
       double precision alpha_n_Ve, beta_n_Ve, alpha_n_Vi, beta_n_Vi
       double precision alpha_m_Ve, beta_m_Ve, alpha_m_Vi, beta_m_Vi
       double precision alpha_h_Ve, beta_h_Ve, alpha_h_Vi, beta_h_Vi
       double precision sigma_Ve, sigma_Vi, m_inf_Ve, m_inf_Vi
       double precision Ie_mem, Ie_syn, Ie_ext, Ie_rand,V_sp
       double precision Ii_mem, Ii_syn, Ii_ext, Ii_rand,vau
       double precision x, y(N*dim), dydx(N*dim)
       double precision thresh,see,sie,sei,sii
       double precision Ik_inf_e, tau_Ik_e, A_IK_e, iksmall_inf_e, p1_e, p2_e
       double precision gamma_Ik_e, alpha_Ik_e, c1, c2, c3, c4, min_e
       double precision beta, I_pump_e, I_diff_e, I_ki_e, I_ko_e, n_k0
       double precision Ik_inf_i, tau_Ik_i, A_IK_i, iksmall_inf_i, p1_i, p2_i
       double precision gamma_Ik_i, alpha_Ik_i, min_i
       double precision I_pump_i, I_diff_i, I_ki_i, I_ko_i, diffusion, deltax
       double precision diffusion_e, diffusion_i,epsilon,G_glia
       double precision gamma_e, gamma_i, gamma_telda, v_b, beta_telda_e, beta_telda_i
       integer seed1, seed2, seed3
	   double precision expotassium_current,exsodium_current,inpotassium_current,insodium_current
	   double precision sodium_oe, sodium_oi, current_conversion
	   double precision Kin_e, Kin_i
       double precision V_l_e, V_l_i

       beta = 7.0         
       n_k0 = 14.0				            	           !n_k0 is represented as k_o,infinity in the paper
       diffusion = 250.0			                       !diffusion coefficient D
       epsilon = 0.07/0.3                                  !rest of the variables as given in the paper
	   G_glia = 5.0/0.3
       deltax = 10.0       
	   v_b = -50.0
       gamma_telda = 0.4			                       
       vau=5
       C = 1.0
       phi = 3.0
       tau_e = 4.0
       A = 20.0
       g_ca = 0.1   
       V_ca = 120.0
       g_l = 0.05
       V_l = -65.0
       g_na = 100.0 
       V_na = 55.0
       g_k = 40.0
       g_ahp = 0.01
       V_k = -80.0
       tau_i = 8.0      
       V_ee = 0.0
       V_ie = -80.0
       V_ei = 0.0
       V_ii = -80.0
       V_sp=40
       specnum1 = sqrt(100.0/PI)
       specnum2 = sqrt(30.0/PI)
       thresh=0.0
       alpha_ee = 0.12 
	   alpha_ie = 0.06     
       alpha_ei = 0.2     
       alpha_ii = 0.02

       alpha_g=0.0					                      !coupling term for gap junctions NOT INCLUDED IN THE MODEL
	
       do  ii=1, N                                        !Loop for 100 excitatory and 100 inhibitory cells
      
		 seed1 = (-1)*int(ii*(x*100))
         seed2=(-1)*int((N-ii)*(x*100))

!        rates for various ions gatting variables, as in paper. The subscripts e and i respectively 
!        represent variables for excitatory and inhibitory network.

         alpha_n_Ve = 0.01 * (y(((ii-1)*dim)+1)+34.0)/( 1.0 - exp(-0.1 * (y(((ii-1)*dim)+1)+34.0)) )
         beta_n_Ve = 0.125 * exp(-1.0*(y(((ii-1)*dim)+1)+44.0)/80.0)
         alpha_n_Vi = 0.01 * (y(((ii-1)*dim)+6)+34.0)/( 1.0 - exp(-0.1 * (y(((ii-1)*dim)+6)+34.0)) )
         beta_n_Vi = 0.125 * exp(-1.0*(y(((ii-1)*dim)+6)+44.0)/80.0)
         alpha_m_Ve = 0.1 * (y(((ii-1)*dim)+1)+30.0)/( 1.0 - exp(-0.1 * (y(((ii-1)*dim)+1)+30.0)) )
         beta_m_Ve = 4.0 * exp(-1.0*(y(((ii-1)*dim)+1)+55.0)/18.0)
         alpha_m_Vi = 0.1 * (y(((ii-1)*dim)+6)+30.0)/( 1.0 - exp(-0.1 * (y(((ii-1)*dim)+6)+30.0)) )
         beta_m_Vi = 4.0 * exp(-1.0*(y(((ii-1)*dim)+6)+55.0)/18.0)
         alpha_h_Ve = 0.07 * exp(-1.0*(y(((ii-1)*dim)+1)+44.0)/20.0)
         beta_h_Ve = 1.0/( 1.0 + exp(-0.1 * (y(((ii-1)*dim)+1)+14.0)) )
         alpha_h_Vi = 0.07 * exp(-1.0*(y(((ii-1)*dim)+6)+44.0)/20.0)
         beta_h_Vi = 1.0/( 1.0 + exp(-0.1 * (y(((ii-1)*dim)+6)+14.0)) )
         sigma_Ve = 1.0/( 1.0 + exp(-1.0*(y(((ii-1)*dim)+1)+20.0)/4.0) )
         sigma_Vi = 1.0/( 1.0 + exp(-1.0*(y(((ii-1)*dim)+6)+20.0)/4.0) )
         m_inf_Ve = alpha_m_Ve/(alpha_m_Ve + beta_m_Ve)
         m_inf_Vi = alpha_m_Vi/(alpha_m_Vi + beta_m_Vi)

!        Potassium and Sodium dynamics, I_diff and I_glia from the paper are combined into I_diff here

         I_pump_e = (1.25/(1.0+exp((25.0-y(((ii-1)*dim)+11))/3.0)))*(1.0/(1.0+exp(8.0-y(((ii-1)*dim)+10)))) 
         I_diff_e = epsilon*(y(((ii-1)*dim)+10)-n_k0)+G_glia/(1.0d0 + exp((18.0-y(((ii-1)*dim)+10))/2.5d0))
         I_pump_i = (1.25/(1.0+exp((25.0-y(((ii-1)*dim)+13))/3.0)))*(1.0/(1.0+exp(8.0-y(((ii-1)*dim)+12)))) 
         I_diff_i = epsilon*(y(((ii-1)*dim)+12)-n_k0)+G_glia/(1.0d0 + exp((18.0-y(((ii-1)*dim)+12))/2.5d0))
		 Kin_e=140.0+(18.0-y(((ii-1)*dim)+11))    
		 Kin_i=140.0+(18.0-y(((ii-1)*dim)+13))         
         sodium_oe = 144.0 - beta * (y(((ii-1)*dim)+11) - 18.0)  
		 sodium_oi = 144.0 - beta * (y(((ii-1)*dim)+13) - 18.0)

!        reversal potentials are updated based on instantaneous ion concentrations

         V_k_e = 26.64 * log(y(((ii-1)*dim)+10)/Kin_e) 
         V_k_i = 26.64 * log(y(((ii-1)*dim)+12)/Kin_i) 
         V_na_e = 26.64 * log(sodium_oe/y(((ii-1)*dim)+11)) 
         V_na_i = 26.64 * log(sodium_oi/y(((ii-1)*dim)+13)) 
         V_l_e = 26.64 * log ( (y(((ii-1)*dim)+10) + 0.065*sodium_oe + 0.6*6.0) / (Kin_e + 0.065*y(((ii-1)*dim)+11) + 0.6*130.0))
         V_l_i = 26.64 * log ( (y(((ii-1)*dim)+12) + 0.065*sodium_oi + 0.6*6.0) / (Kin_i + 0.065*y(((ii-1)*dim)+13) + 0.6*130.0))

!        diffusion of potassium from cell's extracellular space to the nearest neighbours         

         if(ii .eq. 1)then
            diffusion_e = diffusion*(y(((ii)*dim)+10)+y(((ii-1)*dim)+12)-2.0*y(((ii-1)*dim)+10))/(deltax*deltax)
            diffusion_i = diffusion*(y(((ii)*dim)+12)+y(((ii-1)*dim)+10)-2.0*y(((ii-1)*dim)+12))/(deltax*deltax)
         else if(ii .eq. N)then
            diffusion_e = diffusion*(y(((ii-2)*dim)+10)+y(((ii-1)*dim)+12)-2.0*y(((ii-1)*dim)+10))/(deltax*deltax)
            diffusion_i = diffusion*(y(((ii-2)*dim)+12)+y(((ii-1)*dim)+10)-2.0*y(((ii-1)*dim)+12))/(deltax*deltax)
         else
            diffusion_e = diffusion*(y(((ii-2)*dim)+10)+y(((ii)*dim)+10)+y(((ii-1)*dim)+12)-3.0*y(((ii-1)*dim)+10))/(deltax*deltax)
            diffusion_i = diffusion*(y(((ii-2)*dim)+12)+y(((ii)*dim)+12)+y(((ii-1)*dim)+10)-3.0*y(((ii-1)*dim)+12))/(deltax*deltax)
         end if

!        membrane currents for excitatory cells

         Ie_mem = -g_l * (y(((ii-1)*dim)+1) - V_l_e)            &
		 &       - g_na * (m_inf_Ve * m_inf_Ve * m_inf_Ve) * y(((ii-1)*dim)+3) * (y(((ii-1)*dim)+1)-V_na_e)        &
		 &       - (g_k * (y(((ii-1)*dim)+2) * y(((ii-1)*dim)+2) * y(((ii-1)*dim)+2) * y(((ii-1)*dim)+2))          &
		 &               + g_ahp * y(((ii-1)*dim)+5)/(1.0+y(((ii-1)*dim)+5)) )*(y(((ii-1)*dim)+1)-V_k_e)         
		 expotassium_current = (g_k * (y(((ii-1)*dim)+2) * y(((ii-1)*dim)+2) * y(((ii-1)*dim)+2) * y(((ii-1)*dim)+2))    &
		 &               + g_ahp * y(((ii-1)*dim)+5)/(1.0+y(((ii-1)*dim)+5)) )*(y(((ii-1)*dim)+1)-V_k_e)    
		 exsodium_current = g_na * (m_inf_Ve * m_inf_Ve * m_inf_Ve) * y(((ii-1)*dim)+3) * (y(((ii-1)*dim)+1)-V_na_e)   
		 
!******* synaptic current for excitatory cells; if uncoupled, use Ie_syn = 0.0; ***************/

         thresh=0
		 synsum1 = 0.0
         synsum2 = 0.0
         do jj=1, N
!                if((x .gt. 400.0).and.(x .lt. 430.0))THEN                   !use this loop for figure 4
!                    alpha_ee = 0.25                           
!                else
!                    alpha_ee = 0.2165
!                end if
    			 see=0
	      		 sie=0 
			     dx1=float(ii-jj)/float(N)
				 dx2=float(ii+(N-jj))/float(N)
                 g_ee = alpha_ee * specnum1 * (exp(-100.0*(dx1*dx1))+exp(-100.0*(dx2*dx2)))  
                 g_ie = alpha_ie * specnum2 * (exp(-30.0*(dx1*dx1))+exp(-30.0*(dx2*dx2))) 
			     if(ii .eq. jj)g_ee=0
			     if(y(((jj-1)*dim)+1) .gt. thresh)see=1
				 if(y(((jj-1)*dim)+6) .gt. thresh)sie=1
				 if(y(((jj-1)*dim)+14) .gt. 5.0)THEN
				   beta_telda_e = y(((jj-1)*dim)+14)
				 else
				   beta_telda_e = 0.0
				 end if  
				 if(y(((jj-1)*dim)+15) .gt. 5.0)THEN
				   beta_telda_i = y(((jj-1)*dim)+15)
				 else
				   beta_telda_i = 0.0
				 end if  
                 synsum1 = synsum1 + g_ee * y(((jj-1)*dim)+4) * exp(-beta_telda_e/vau)
                 synsum2 = synsum2 + g_ie * y(((jj-1)*dim)+9) * exp(-beta_telda_i/vau)
		 end do		 
         Ie_syn = - (y(((ii-1)*dim)+1) - V_ee) * synsum1 / float(N)         &
         &        - (y(((ii-1)*dim)+1) - V_ie) * synsum2 / float(N)

!!!**********END synaptic currents for excitatory cells**************!!!!!!!!

!		 if ( ((x .gt. 12.0) .and. (x .lt. 32.0)).and.((ii .ge. 21).and.(ii .le. 79)))then                  !use this loop for figure 2,3, 4
!             Ie_ext = 1.5*exp(-60.0*( (float(ii)-float(N)/2.0)/float(N) )*( (float(ii)-float(N)/2.0)/float(N) ) )
!        else if(((x .gt. 400.0).and.(x .le. 420.0)).and.((ii .ge. 1).and.(ii .le. 100)))then
!             Ie_ext = 0.85d0        !1.49,.99        !2nd stim 0.5d0
!        else
		      Ie_ext = 0.0
!		 end if

         Ie_rand = 20*(0.5-ran2(seed1))

!         /* currents into the inhibitory neurons */

         Ii_mem = - g_l * (y(((ii-1)*dim)+6) - V_l_i)           &
         &        - g_na * (m_inf_Vi * m_inf_Vi * m_inf_Vi) * y(((ii-1)*dim)+8) * (y(((ii-1)*dim)+6)-V_na_i)       &
         &        - g_k * (y(((ii-1)*dim)+7) * y(((ii-1)*dim)+7) * y(((ii-1)*dim)+7) * y(((ii-1)*dim)+7))          &
         &          * (y(((ii-1)*dim)+6)-V_k_i)
	     inpotassium_current = g_k * (y(((ii-1)*dim)+7) * y(((ii-1)*dim)+7) * y(((ii-1)*dim)+7) * y(((ii-1)*dim)+7))    &
         &                     * (y(((ii-1)*dim)+6)-V_k_i)  
		 insodium_current = g_na * (m_inf_Vi * m_inf_Vi * m_inf_Vi) * y(((ii-1)*dim)+8) * (y(((ii-1)*dim)+6)-V_na_i)  
         
!******** synaptic current into inhibitory neurons; if uncoupled, use Ii_syn = 0.0; *********************/

         thresh=0
		 synsum1 = 0.0
         synsum2 = 0.0
		 gapsum=0.0
		 do jj = 1, N
			     sei=0
			     sii=0
			     dx1=float(ii-jj)/float(N)
				 dx2=float(ii+(N-jj))/float(N)
				 
				 if(y(((jj-1)*dim)+1) .gt. thresh) sei=1
				 if(y(((jj-1)*dim)+6) .gt. thresh) sii=1
                 g_ei = alpha_ei * specnum2 * (exp( -30.0*(dx1*dx1) )+exp( -30.0*(dx2*dx2) )) 
                 g_ii = alpha_ii * specnum2 * (exp( -30.0*(dx1*dx1) )+exp( -30.0*(dx2*dx2) ))
!				 g_g=alpha_g*specnum2 * (exp( -30.0*(dx1*dx1) )+exp( -30.0*(dx2*dx2) ))
				 if(jj .eq. ii)g_ii=0
				 if(y(((jj-1)*dim)+14) .gt. 5.0)THEN
				   beta_telda_e = y(((jj-1)*dim)+14)
				 else
				   beta_telda_e = 0.0
				 end if  
				 if(y(((jj-1)*dim)+15) .gt. 5.0)THEN
				   beta_telda_i = y(((jj-1)*dim)+15)
				 else
				   beta_telda_i = 0.0
				 end if  
                 synsum1 = synsum1 + g_ei * y(((jj-1)*dim)+4) * exp(-beta_telda_e/vau)
                 synsum2 = synsum2 + g_ii * y(((jj-1)*dim)+9) * exp(-beta_telda_i/vau)
!                gapsum=gapsum+g_g*(y(((jj-1)*dim)+6)-y(((ii-1)*dim)+6))
		 end do
         Ii_syn = - (y(((ii-1)*dim)+6) - V_ei) * synsum1 / float(N)       &
         &        - (y(((ii-1)*dim)+6) - V_ii) * synsum2 / float(N)         ! &
!		 &		 +gapsum/float(N)

! *********END synaptic currents into inhibitory cells*****************!!!!!!!!!!!

		 Ii_ext = 0.5                  
		 Ii_rand =20.0*(0.5-ran2(seed2)) 

!!!!!!!!!variables modeling the interplay between neurons (see paper)
        
		 gamma_e = 0.0
         gamma_i = 0.0
         if((y(((ii-1)*dim)+1) .gt. -30.0).and.(y(((ii-1)*dim)+1) .lt. -10.0))THEN
            gamma_e = 0.4
         end if   
         if((y(((ii-1)*dim)+6) .gt. -30.0).and.(y(((ii-1)*dim)+6) .lt. -10.0))THEN
            gamma_i = 0.4
         end if         
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!         

!         /* excitatory neurons */

         dydx(((ii-1)*dim)+1) = (1/C) * (Ie_mem + Ie_ext + Ie_rand + Ie_syn)
         dydx(((ii-1)*dim)+2) = phi * ( alpha_n_Ve * (1.0-y(((ii-1)*dim)+2)) - beta_n_Ve * y(((ii-1)*dim)+2) )
         dydx(((ii-1)*dim)+3) = phi * ( alpha_h_Ve * (1.0-y(((ii-1)*dim)+3)) - beta_h_Ve * y(((ii-1)*dim)+3) )
         dydx(((ii-1)*dim)+4) = (1/tau_e) * (A * sigma_Ve * (1.0 - y(((ii-1)*dim)+4)) - y(((ii-1)*dim)+4))
         dydx(((ii-1)*dim)+5) = -0.002 * g_ca * (y(((ii-1)*dim)+1) - V_ca)/( 1.0 + exp( -1.0*(y(((ii-1)*dim)+1)+25.0)/2.5 ) ) - y(((ii-1)*dim)+5)/80.0

!         /* inhibitory neurons */

         dydx(((ii-1)*dim)+6) = (1/C) * (Ii_mem + Ii_ext + Ii_rand + Ii_syn)
         dydx(((ii-1)*dim)+7) = phi * ( alpha_n_Vi * (1.0-y(((ii-1)*dim)+7)) - beta_n_Vi * y(((ii-1)*dim)+7) )
         dydx(((ii-1)*dim)+8) = phi * ( alpha_h_Vi * (1.0-y(((ii-1)*dim)+8)) - beta_h_Vi * y(((ii-1)*dim)+8) )
         dydx(((ii-1)*dim)+9) = (1/tau_i) * (A * sigma_Vi * (1.0 - y(((ii-1)*dim)+9)) - y(((ii-1)*dim)+9))

!        Extracellular Potassium, Intracellular Sodium equations
!        Excitatory Neuron

         dydx(((ii-1)*dim)+10) = 0.001*(expotassium_current/3.0 - 7.0* 2.0*I_pump_e - I_diff_e + diffusion_e)        !factor 0.001 converts from seconds to milliseconds  
         dydx(((ii-1)*dim)+11) = 0.001*(-1.0*exsodium_current/(7.0*3.0) - 3.0*I_pump_e)  

!        Extracellular Potassium, Intracellular Sodium equations         
!        Inhibitory Neuron         

         dydx(((ii-1)*dim)+12) = 0.001*(inpotassium_current/3.0 - 7.0* 2.0*I_pump_i - I_diff_i + diffusion_i)
         dydx(((ii-1)*dim)+13) = 0.001*(-1.0*insodium_current/(7.0*3.0) - 3.0*I_pump_i)

!        The depolarization variable
!        Excitatory Neuron

         dydx(((ii-1)*dim)+14) = gamma_e * (y(((ii-1)*dim)+1) - v_b) - gamma_telda * y(((ii-1)*dim)+14)

!        The depolarization variable         
!        Inhibitory Neuron

         dydx(((ii-1)*dim)+15) = gamma_i * (y(((ii-1)*dim)+6) - v_b) - gamma_telda * y(((ii-1)*dim)+15)

       end do                             !end LOOP for 100 excitatory and 100 inhibitory cells
       end

       SUBROUTINE rk4(y,dydx,n,x,h,yout,derivs) 
       INTEGER n
       double precision h,x,dydx(n),y(n),yout(n) 
       EXTERNAL derivs 
       INTEGER i 
       double precision h6,hh,xh,dym(n),dyt(n),yt(n) 
       hh=h*0.5 
       h6=h/6. 
       xh=x+hh

       do i = 1, n     ! First step. 
          yt(i)=y(i)+hh*dydx(i) 
       end do  
       call derivs(xh,yt,dyt)    !  Secondstep. 
       do i = 1, n 
          yt(i)=y(i)+hh*dyt(i) 
       end do  
       call derivs(xh,yt,dym)    !  Thirdstep. 
       do i = 1, n 
          yt(i)=y(i)+h*dym(i) 
          dym(i)=dyt(i)+dym(i) 
       end do  
       call derivs(x+h,yt,dyt)   !  Fourthstep. 
       do i = 1, n          !  Accumulate increments with proper weights. 
          yout(i)=y(i)+h6*(dydx(i)+dyt(i)+2.*dym(i)) 
       end do 
       END 

       FUNCTION ran2(idum)
       INTEGER idum,IM1,IM2,IMM1,IA1,IA2,IQ1,IQ2,IR1,IR2,NTAB,NDIV
       REAL ran2,AM,EPS,RNMX
       PARAMETER (IM1=2147483563,IM2=2147483399,AM=1./IM1,      & 
       &     IMM1=IM1-1,IA1=40014,IA2=40692,IQ1=53668,          &
       &     IQ2=52774,IR1=12211,IR2=3791,NTAB=32,              &
       &     NDIV=1+IMM1/NTAB,EPS=1.2e-7,RNMX=1.-EPS)
       INTEGER idum2,j,k,iv(NTAB),iy
       SAVE iv,iy,idum2
       DATA idum2/123456789/, iv/NTAB*0/, iy/0/
          if (idum.le.0)then
             idum=max(-idum,1) 
             idum2=idum
             do  j=NTAB+8,1,-1 
                k=idum/IQ1
                idum=IA1*(idum-k*IQ1)-k*IR1
                   if (idum.lt.0) idum=idum+IM1
                     if (j .le. NTAB) iv(j)=idum
             end do 
             iy=iv(1)
          end if
          k=idum/IQ1
          idum=IA1*(idum-k*IQ1)-k*IR1 
          if (idum.lt.0) idum=idum+IM1
          k=idum2/IQ2
          idum2=IA2*(idum2-k*IQ2)-k*IR2
          if (idum2.lt.0) idum2=idum2+IM2
          j=1+iy/NDIV
          iy=iv(j)-idum2 
          iv(j)=idum
          if(iy.lt.1)iy=iy+IMM1
          ran2=min(AM*iy,RNMX) 
       END

